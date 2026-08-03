// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <odia/Library.h>

#include <cstddef>
#include <string>
#include <vector>

namespace ODIA
{

  /// How decoy precursors are constructed.
  ///
  /// Kept selectable rather than fixed. DIA-NN's mutation scheme is the default
  /// (D7), but its own documentation identifies it as the weak point of the FDR
  /// model -- the decoy keeps the target's precursor m/z, iRT and library
  /// intensity pattern, so it is searched in the same isolation window over the
  /// same RT range with the same expected shape, and only fragment masses
  /// differ. Phase 3 needs to be able to measure that against the alternative
  /// under entrapment rather than inherit it.
  enum class DecoyMethod
  {
    None,
    Mutate,        ///< DIA-NN: substitute one residue near each terminus
    PseudoReverse  ///< reverse all but the C-terminal residue
  };

  DecoyMethod parseDecoyMethod(const std::string& s);

  struct DigestParams
  {
    /// "Trypsin/P" rather than "Trypsin": cleave after K and R even when the
    /// next residue is proline.
    ///
    /// The no-cut-before-proline rule is a partial one, and search engines
    /// have largely abandoned it -- DIA-NN's default is `--cut K*,R*`, which
    /// cuts regardless. Measured against DIA-NN's library on the human
    /// proteome, keeping "Trypsin" cost 169,044 peptides, and 48.6% of those
    /// began with a proline while none of ours did. That was the single
    /// largest source of the coverage gap, well ahead of methionine excision.
    std::string enzyme = "Trypsin/P";
    std::size_t missed_cleavages = 1;
    std::size_t min_length = 7;
    std::size_t max_length = 30;

    std::vector<int> charges{2, 3};

    /// UniMod-style names as OpenMS knows them, e.g. "Carbamidomethyl (C)".
    std::vector<std::string> fixed_modifications{"Carbamidomethyl (C)"};
    std::vector<std::string> variable_modifications{};
    std::size_t max_variable_modifications = 1;

    /// Digest the sequence again with the initiator methionine removed.
    ///
    /// Not an option in practice, though it is one here. The initiator Met is
    /// cleaved co-translationally from a large fraction of proteins, so the
    /// observed N-terminal peptide usually begins one residue in. Leaving this
    /// off cost 272,565 precursors against DIA-NN on the human proteome --
    /// 7.1% of the library, every one of them inside the m/z window and evenly
    /// split across charge, which is what identified it: a filter would not be
    /// charge-neutral.
    bool n_terminal_methionine_excision = true;

    double precursor_mz_min = 350.0;
    double precursor_mz_max = 1200.0;

    double fragment_mz_min = 200.0;
    double fragment_mz_max = 1800.0;
    int max_fragment_charge = 2;

    /// A precursor with fewer usable fragments than this is not searchable.
    std::size_t min_fragments = 4;
    /// Cap per precursor, taken in descending fragment m/z.
    std::size_t max_fragments = 12;

    DecoyMethod decoy_method = DecoyMethod::Mutate;
  };

  /// Builds an assay library from protein sequences.
  ///
  /// Produces the peptide, precursor and fragment skeleton with placeholder
  /// intensities and retention times; those are filled in by prediction.
  class LibraryGenerator
  {
  public:
    struct Stats
    {
      std::size_t proteins = 0;
      std::size_t peptides = 0;
      std::size_t precursors = 0;
      std::size_t transitions = 0;
      std::size_t dropped_too_few_fragments = 0;
      std::size_t dropped_precursor_mz = 0;
    };

    static Stats generate(const std::string& fasta_file,
                          const DigestParams& params,
                          Library& library);

    /// Fill in predicted retention times, replacing the placeholders.
    ///
    /// Predicts once per distinct modified sequence rather than once per
    /// precursor: the RT model takes no charge input, so the charge states of
    /// one peptide would otherwise be predicted identically several times over.
    /// On the human proteome that is 1.25 M predictions instead of 4.0 M.
    ///
    /// @returns the number of precursors whose iRT could not be predicted;
    ///          theirs are left NaN rather than given a made-up value.
    static std::size_t predictRetentionTimes(Library& library,
                                             const std::string& rt_model_path,
                                             bool prefer_gpu = true);

    /// Replace placeholder intensities with predicted ones, and re-choose the
    /// fragments now that there is a basis for ranking.
    ///
    /// This necessarily rebuilds the transition arrays rather than editing
    /// them: generate() caps fragments by descending m/z because it has no
    /// intensities, and the top twelve by m/z are not the top twelve by
    /// intensity, so re-ranking what survived that cap would not give the same
    /// answer as ranking before it.
    ///
    /// Predicts in blocks. A whole-proteome call would materialise ~4 M
    /// spectra -- about 2.2 GiB of payload in as many allocations -- all live
    /// at once, on top of the library being built.
    ///
    /// Call before appendDecoys: a decoy copies its target's intensity
    /// pattern, so predicting afterwards would leave every decoy with the
    /// placeholder.
    ///
    /// @returns the number of precursors left with placeholder intensities
    ///          because prediction failed for them.
    static std::size_t predictFragmentIntensities(Library& library,
                                                  const std::string& ms2_model_path,
                                                  const DigestParams& params,
                                                  float nce = 30.0f,
                                                  const std::string& instrument = "QE",
                                                  bool prefer_gpu = true);

    /// The line mapping the RT model's raw output onto the iRT scale.
    struct IrtCalibration
    {
      double slope = 1.0;
      double intercept = 0.0;
      std::size_t peptides = 0;      ///< standards that could be predicted
      double max_abs_error = 0.0;    ///< worst standard, in iRT units

      double apply(double raw) const { return slope * raw + intercept; }
    };

    /// Fit the raw-to-iRT line from the Biognosys standards, using @p model.
    ///
    /// Recomputed rather than pinned: the line depends on the model checkpoint,
    /// and a hardcoded one would rot silently on a model swap -- the outputs
    /// would still look like iRT and simply be wrong.
    static IrtCalibration fitIrtCalibration(const std::string& rt_model_path,
                                            const std::string& standards_file,
                                            bool prefer_gpu = true);

    /// Rescale a library's retention times onto the iRT scale in place.
    ///
    /// This does NOT make the predictions more accurate, and must not be
    /// described as though it does. It is monotone, so any consumer that fits
    /// its own retention-time calibration -- DIA-NN does -- is unaffected:
    /// measured on the human proteome, applying it left DIA-NN's search window
    /// unchanged to the last digit (2.18905 both ways). It is applied because a
    /// column named iRT holding a 0..1 training-gradient coordinate is a
    /// silent-units error waiting for a consumer that applies a tolerance in
    /// iRT units without calibrating first.
    static void applyIrtCalibration(Library& library, const IrtCalibration& calibration);

    /// Fill in predicted collision cross-sections, in square angstroms.
    ///
    /// Per precursor, not per peptide: the CCS model takes charge, and a
    /// peptide's charge states have genuinely different cross-sections.
    ///
    /// This does NOT populate the ion-mobility column. CCS and 1/K0 are
    /// different quantities, related through the drift gas and the
    /// instrument's calibration, and that conversion belongs downstream where
    /// the instrument is known.
    ///
    /// @returns the number of precursors left without a value; theirs are NaN.
    static std::size_t predictCollisionCrossSections(Library& library,
                                                    const std::string& ccs_model_path,
                                                    bool prefer_gpu = true);

    /// Append a decoy for every target currently in @p library.
    ///
    /// Decoys keep the target's precursor m/z, iRT and intensity pattern; only
    /// fragment m/z values move. See DecoyMethod.
    /// @param skipped_out receives the number of targets for which no decoy
    ///        could be built. Reported rather than dropped silently: 1.3% of the
    ///        DIA-NN fixture is affected, all of them N-terminally modified.
    /// @param min_fragments the same bar the targets had to clear. Applying it
    ///        to one class only is an anti-conservative FDR (D7 rule 2).
    /// Idempotent: targets that already have a decoy are skipped.
    static std::size_t appendDecoys(Library& library, DecoyMethod method,
                                    std::size_t* skipped_out = nullptr,
                                    std::size_t min_fragments = 0);
  };

} // namespace ODIA
