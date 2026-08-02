// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause
//
// GENERATED from data/peptdeep_mod_elements.txt by
// scripts/generate_peptdeep_elements.py -- do not edit by hand.

#pragma once

#include <array>
#include <string_view>

namespace ODIA
{

  /// AlphaPeptDeep's mod_elements. The index is the feature position.
  inline constexpr std::array<std::string_view, 109> PEPTDEEP_MOD_ELEMENTS{
    "C", "H", "N", "O", "P", "S", "B", "F",
    "I", "K", "U", "V", "W", "X", "Y", "Ac",
    "Ag", "Al", "Am", "Ar", "As", "At", "Au", "Ba",
    "Be", "Bi", "Bk", "Br", "Ca", "Cd", "Ce", "Cf",
    "Cl", "Cm", "Co", "Cr", "Cs", "Cu", "Dy", "Er",
    "Es", "Eu", "Fe", "Fm", "Fr", "Ga", "Gd", "Ge",
    "He", "Hf", "Hg", "Ho", "In", "Ir", "Kr", "La",
    "Li", "Lr", "Lu", "Md", "Mg", "Mn", "Mo", "Na",
    "Nb", "Nd", "Ne", "Ni", "No", "Np", "Os", "Pa",
    "Pb", "Pd", "Pm", "Po", "Pr", "Pt", "Pu", "Ra",
    "Rb", "Re", "Rh", "Rn", "Ru", "Sb", "Sc", "Se",
    "Si", "Sm", "Sn", "Sr", "Ta", "Tb", "Tc", "Te",
    "Th", "Ti", "Tl", "Tm", "Xe", "Yb", "Zn", "Zr",
    "2H", "13C", "15N", "18O", "?",
  };

  /// Index of an element, or the '?' catch-all slot if unknown.
  inline constexpr std::size_t peptDeepElementIndex(std::string_view symbol)
  {
    for (std::size_t i = 0; i < PEPTDEEP_MOD_ELEMENTS.size(); ++i)
    {
      if (PEPTDEEP_MOD_ELEMENTS[i] == symbol) { return i; }
    }
    return PEPTDEEP_MOD_ELEMENTS.size() - 1;
  }

} // namespace ODIA
