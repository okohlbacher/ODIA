// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/LibraryGenerator.h>

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/FORMAT/FASTAFile.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <unordered_map>

using namespace OpenMS;

namespace ODIA
{
  namespace
  {
    /// DIA-NN's substitution table: residue -> replacement.
    ///
    /// Fixed, not random, and reproduced here so the decoys match what the
    /// reference engine would produce.
    constexpr std::string_view MUTATE_FROM = "GAVLIFMPWSCTYHKRQEND";
    constexpr std::string_view MUTATE_TO   = "LLLVVLLLLTSSSSLLNDQE";

    char mutateResidue(char aa)
    {
      const auto pos = MUTATE_FROM.find(aa);
      return pos == std::string_view::npos ? aa : MUTATE_TO[pos];
    }

    /// One fragment before it is written into the library arrays.
    struct Fragment
    {
      double mz;
      FragmentType type;
      std::uint8_t ordinal;
      std::int8_t charge;
    };

    void enumerateFragments(const AASequence& peptide, const DigestParams& p,
                            int precursor_charge, std::vector<Fragment>& out)
    {
      out.clear();
      const Size n = peptide.size();
      if (n < 2) { return; }

      const int max_z = std::min(p.max_fragment_charge, std::max(1, precursor_charge - 1));
      for (int z = 1; z <= max_z; ++z)
      {
        for (Size i = 1; i < n; ++i)
        {
          const double b = peptide.getPrefix(i).getMZ(z, Residue::BIon);
          if (b >= p.fragment_mz_min && b <= p.fragment_mz_max)
          {
            out.push_back({b, FragmentType::B, static_cast<std::uint8_t>(i),
                           static_cast<std::int8_t>(z)});
          }
          const double y = peptide.getSuffix(i).getMZ(z, Residue::YIon);
          if (y >= p.fragment_mz_min && y <= p.fragment_mz_max)
          {
            out.push_back({y, FragmentType::Y, static_cast<std::uint8_t>(i),
                           static_cast<std::int8_t>(z)});
          }
        }
      }

      // Without predicted intensities there is no basis for ranking, so the cap
      // is applied by descending m/z: high-m/z fragments carry more sequence and
      // sit in a less crowded part of the spectrum. Prediction replaces this.
      std::sort(out.begin(), out.end(), [](const Fragment& a, const Fragment& b) {
        if (a.mz != b.mz) { return a.mz > b.mz; }
        if (a.type != b.type) { return a.type < b.type; }
        if (a.ordinal != b.ordinal) { return a.ordinal < b.ordinal; }
        return a.charge < b.charge;
      });
      if (out.size() > p.max_fragments) { out.resize(p.max_fragments); }
    }
  } // namespace

  DecoyMethod parseDecoyMethod(const std::string& s)
  {
    if (s == "none") { return DecoyMethod::None; }
    if (s == "pseudo_reverse") { return DecoyMethod::PseudoReverse; }
    return DecoyMethod::Mutate;
  }

  LibraryGenerator::Stats LibraryGenerator::generate(const std::string& fasta_file,
                                                     const DigestParams& params,
                                                     Library& library)
  {
    Stats stats;

    std::vector<FASTAFile::FASTAEntry> entries;
    FASTAFile().load(fasta_file, entries);
    stats.proteins = entries.size();

    ProteaseDigestion digestion;
    digestion.setEnzyme(params.enzyme);
    digestion.setMissedCleavages(params.missed_cleavages);

    // Peptides are deduplicated across proteins; a peptide seen in several
    // proteins accumulates them, which is what protein inference later needs.
    std::map<std::string, std::string> peptide_to_proteins;
    std::vector<AASequence> peptides;
    for (const auto& entry : entries)
    {
      AASequence protein;
      try { protein = AASequence::fromString(entry.sequence); }
      catch (const std::exception&) { continue; }  // non-standard residues

      peptides.clear();
      digestion.digest(protein, peptides, params.min_length, params.max_length);
      for (const auto& pep : peptides)
      {
        auto& proteins = peptide_to_proteins[pep.toUnmodifiedString()];
        if (proteins.empty()) { proteins = entry.identifier; }
        else if (proteins.find(entry.identifier) == std::string::npos)
        {
          proteins += ';';
          proteins += entry.identifier;
        }
      }
    }
    stats.peptides = peptide_to_proteins.size();

    ModifiedPeptideGenerator::MapToResidueType fixed_map =
      ModifiedPeptideGenerator::getModifications(
        std::vector<std::string>(params.fixed_modifications.begin(),
                                 params.fixed_modifications.end()));
    ModifiedPeptideGenerator::MapToResidueType variable_map =
      ModifiedPeptideGenerator::getModifications(
        std::vector<std::string>(params.variable_modifications.begin(),
                                 params.variable_modifications.end()));

    library.reserve(peptide_to_proteins.size() * params.charges.size(),
                    peptide_to_proteins.size() * params.charges.size() * params.max_fragments);

    std::vector<Fragment> fragments;
    std::vector<AASequence> forms;
    for (const auto& [sequence, proteins] : peptide_to_proteins)
    {
      AASequence base;
      try { base = AASequence::fromString(sequence); }
      catch (const std::exception&) { continue; }

      ModifiedPeptideGenerator::applyFixedModifications(fixed_map, base);

      forms.clear();
      forms.push_back(base);
      if (!params.variable_modifications.empty())
      {
        ModifiedPeptideGenerator::applyVariableModifications(
          variable_map, base, params.max_variable_modifications, forms, true);
      }

      const std::uint32_t protein_handle = library.strings().intern(proteins);

      for (const auto& form : forms)
      {
        const std::uint32_t sequence_handle =
          library.strings().intern(form.toString());

        for (int z : params.charges)
        {
          const double precursor_mz = form.getMZ(z);
          if (precursor_mz < params.precursor_mz_min || precursor_mz > params.precursor_mz_max)
          {
            ++stats.dropped_precursor_mz;
            continue;
          }

          enumerateFragments(form, params, z, fragments);
          if (fragments.size() < params.min_fragments)
          {
            ++stats.dropped_too_few_fragments;
            continue;
          }

          auto& p = library.precursors();
          auto& t = library.transitions();
          p.mz.push_back(toFixed(precursor_mz));
          p.irt.push_back(std::nanf(""));   // filled in by prediction
          p.im.push_back(std::nanf(""));
          p.charge.push_back(static_cast<std::uint8_t>(z));
          p.decoy.push_back(0);
          p.modified_sequence.push_back(sequence_handle);
          p.protein_group.push_back(protein_handle);
          p.transition_begin.push_back(static_cast<std::uint32_t>(t.product_mz.size()));
          p.transition_count.push_back(static_cast<std::uint32_t>(fragments.size()));

          for (const auto& f : fragments)
          {
            t.product_mz.push_back(toFixed(f.mz));
            t.library_intensity.push_back(1.0f);  // placeholder until prediction
            t.type.push_back(f.type);
            t.ordinal.push_back(f.ordinal);
            t.charge.push_back(f.charge);
            t.loss.push_back(LossType::None);
          }
        }
      }
    }

    stats.precursors = library.precursorCount();
    stats.transitions = library.transitionCount();
    return stats;
  }

  std::size_t LibraryGenerator::appendDecoys(Library& library, DecoyMethod method)
  {
    if (method == DecoyMethod::None) { return 0; }

    const std::size_t n_targets = library.precursorCount();
    std::size_t made = 0;

    for (std::size_t i = 0; i < n_targets; ++i)
    {
      if (library.precursors().decoy[i]) { continue; }

      const std::string sequence(library.strings().get(library.precursors().modified_sequence[i]));

      // Tokenise into residues, each carrying its own modification suffix, so
      // mutation and reversal preserve modifications instead of dropping them.
      struct Token { std::string text; char residue; bool modified; };
      std::vector<Token> tokens;
      for (std::size_t k = 0; k < sequence.size(); ++k)
      {
        if (sequence[k] == '(' || sequence[k] == '[') { continue; }  // leading N-term mod
        Token tok{std::string(1, sequence[k]), sequence[k], false};
        while (k + 1 < sequence.size() && (sequence[k + 1] == '(' || sequence[k + 1] == '['))
        {
          const char close = sequence[k + 1] == '(' ? ')' : ']';
          std::size_t j = k + 1;
          while (j < sequence.size() && sequence[j] != close) { tok.text.push_back(sequence[j++]); }
          if (j < sequence.size()) { tok.text.push_back(sequence[j]); }
          tok.modified = true;
          k = j;
        }
        tokens.push_back(std::move(tok));
      }
      if (tokens.size() < 4) { continue; }

      std::string decoy_sequence;

      if (method == DecoyMethod::Mutate)
      {
        // DIA-NN picks positions near each terminus that are not modified;
        // mutating a modified residue would silently invalidate the modification.
        std::size_t n_pos = tokens.size();
        for (std::size_t k = 1; k + 1 < tokens.size(); ++k)
        {
          if (!tokens[k].modified) { n_pos = k; break; }
        }
        std::size_t c_pos = tokens.size();
        for (std::size_t k = tokens.size() - 2; k >= 1; --k)
        {
          if (!tokens[k].modified && k != n_pos) { c_pos = k; break; }
          if (k == 1) { break; }
        }
        if (n_pos == tokens.size() || c_pos == tokens.size()) { continue; }

        tokens[n_pos].text = std::string(1, mutateResidue(tokens[n_pos].residue));
        tokens[c_pos].text = std::string(1, mutateResidue(tokens[c_pos].residue));
        for (const auto& tok : tokens) { decoy_sequence += tok.text; }
      }
      else // PseudoReverse: reverse everything but the C-terminal residue
      {
        std::vector<Token> reordered(tokens.begin(), tokens.end() - 1);
        std::reverse(reordered.begin(), reordered.end());
        reordered.push_back(tokens.back());
        for (const auto& tok : reordered) { decoy_sequence += tok.text; }
      }

      // Recompute every fragment from the decoy sequence rather than shifting
      // the target's.
      //
      // Shifting b ions by the N-terminal mass difference and y ions by the
      // C-terminal one is wrong whenever an ion spans *both* mutated residues,
      // which the long b and y ions always do -- and those are exactly the ions
      // enumerateFragments keeps, since it caps by descending m/z. Measured on
      // the human proteome before this change: 15.3% of decoy fragment m/z did
      // not correspond to the decoy sequence stored beside them, worst case
      // 76 Th out. The failure is worse than it sounds because it is asymmetric:
      // recomputing from the sequence disagrees on ~15% of decoy transitions and
      // 0% of target ones, and a criterion that behaves differently for the two
      // classes is the anti-conservative FDR mode D7 exists to avoid.
      AASequence decoy;
      try { decoy = AASequence::fromString(decoy_sequence); }
      catch (const std::exception&) { continue; }

      auto& p = library.precursors();
      auto& t = library.transitions();

      const std::uint32_t begin = p.transition_begin[i];
      const std::uint32_t count = p.transition_count[i];
      const std::uint32_t new_begin = static_cast<std::uint32_t>(t.product_mz.size());

      for (std::uint32_t k = 0; k < count; ++k)
      {
        const std::uint32_t s = begin + k;
        const auto charge = t.charge[s] == 0 ? std::int8_t{1} : t.charge[s];
        const Size ord = t.ordinal[s];
        if (ord == 0 || ord >= decoy.size()) { continue; }

        double mz = 0.0;
        try
        {
          switch (t.type[s])
          {
            case FragmentType::A: mz = decoy.getPrefix(ord).getMZ(charge, Residue::AIon); break;
            case FragmentType::B: mz = decoy.getPrefix(ord).getMZ(charge, Residue::BIon); break;
            case FragmentType::C: mz = decoy.getPrefix(ord).getMZ(charge, Residue::CIon); break;
            case FragmentType::X: mz = decoy.getSuffix(ord).getMZ(charge, Residue::XIon); break;
            case FragmentType::Y: mz = decoy.getSuffix(ord).getMZ(charge, Residue::YIon); break;
            case FragmentType::Z: mz = decoy.getSuffix(ord).getMZ(charge, Residue::ZIon); break;
            default:
              // Precursor and unrecognised types cannot be recomputed. Copying
              // the target's value would give the decoy a transition identical
              // to its target's, so drop it instead.
              continue;
          }
        }
        catch (const std::exception&) { continue; }

        t.product_mz.push_back(toFixed(mz));
        t.library_intensity.push_back(t.library_intensity[s]);
        t.type.push_back(t.type[s]);
        t.ordinal.push_back(t.ordinal[s]);
        t.charge.push_back(t.charge[s]);
        t.loss.push_back(t.loss[s]);
      }

      const auto written = static_cast<std::uint32_t>(t.product_mz.size() - new_begin);
      if (written == 0) { continue; }

      p.mz.push_back(p.mz[i]);
      p.irt.push_back(p.irt[i]);
      p.im.push_back(p.im[i]);
      p.charge.push_back(p.charge[i]);
      p.decoy.push_back(1);
      p.modified_sequence.push_back(library.strings().intern(decoy_sequence));
      p.protein_group.push_back(p.protein_group[i]);
      p.transition_begin.push_back(new_begin);
      p.transition_count.push_back(written);
      ++made;
    }
    return made;
  }

} // namespace ODIA
