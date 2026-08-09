// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <odia/scoring/lda.h>

#include <cstdint>
#include <string>
#include <vector>

namespace ODIA
{
  namespace Scoring
  {

    /// Percolator as a drop-in scoring engine, via OpenMS 3.6's in-process API.
    ///
    /// Signature-compatible with `scoreSemiSupervisedLDA` on purpose: the two
    /// take the same (features, labels, group) and return the same
    /// `ScoredGroups`, so an engine is a one-line substitution at the call site
    /// and the rest of the pipeline cannot tell them apart. That IS the
    /// plug-and-play interface -- there is no benefit in a class hierarchy over
    /// a common function signature, and a hierarchy would make the engines
    /// harder to A/B rather than easier.
    ///
    /// WHY PERCOLATOR AND NOT MOKAPOT. mokapot (Fondrie & Noble) is a Python
    /// reimplementation of Percolator's semi-supervised SVM with cross-
    /// validation; OpenMS 3.6 ships Percolator itself, in process, with a
    /// domain-agnostic `RescoreInput`. For a C++ tool that is strictly better
    /// than mokapot: the same algorithm family, no Python dependency, no
    /// process boundary, no on-disk interchange format to keep in sync. There
    /// is no mokapot implementation in this OpenMS install -- checked.
    ///
    /// WHAT DIFFERS FROM OUR NATIVE SCORER. Both are semi-supervised
    /// target-decoy learners of the mProphet lineage, but:
    ///
    ///   * Percolator cross-validates by construction, training on folds and
    ///     scoring the held-out one, which is the standard defence against the
    ///     discriminant memorising its own training positives. Our LDA/GBT path
    ///     does its own fold handling and has never been compared against an
    ///     independent implementation of the same idea.
    ///   * It is a linear SVM, so it cannot represent an interaction that our
    ///     GBT can. On features that are individually weak but jointly
    ///     informative it should LOSE; on a small, noisy positive class it may
    ///     win precisely because it cannot overfit as hard.
    ///
    /// That second point is the reason to try it here. ODIA's open failure is a
    /// proteome-scale library where ~1.5% of targets are present, the positive
    /// class is ~98.5% noise, and the semi-supervised loop never ignites.
    ///
    /// `cv_group_keys` carries the precursor grouping, so the several candidate
    /// peak groups of one precursor cannot be split across folds -- which would
    /// leak, since they share a chromatogram.
    /// @param model_out  if non-empty, TRAIN on this data and save the model here
    /// @param model_in   if non-empty, do not train: load this model and apply it
    ///
    /// STATIC MODELLING, which is the point of the two paths.
    ///
    /// A semi-supervised scorer needs confident positives to bootstrap, and at a
    /// 1.5% true-positive rate it does not have them: GBT, LDA and Percolator
    /// all report zero on `v6_50k` while the SAME discriminant, ranked, puts 232
    /// of DIA-NN's 738 in its top 738 (31.4% precision against a 1.5% base
    /// rate). The signal is there and the q-values cannot certify it.
    ///
    /// The Noble/Keich/Kall remedy is to stop asking the hard run to train
    /// itself: fit the model where the loop DOES ignite, freeze it, and apply
    /// it. That is not a trick -- a discriminant over these sub-scores is a
    /// property of the instrument and the scoring code, not of which peptides
    /// happen to be in one sample.
    ScoredGroups scorePercolator(const std::vector<std::vector<double>>& features,
                                 const std::vector<int>& labels,
                                 const std::vector<long long>& group,
                                 const LDAParams& params,
                                 const std::vector<std::string>& feature_names,
                                 std::string* diagnostic = nullptr,
                                 const std::string& model_out = "",
                                 const std::string& model_in = "");

  } // namespace Scoring
} // namespace ODIA
