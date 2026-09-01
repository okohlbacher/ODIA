# 82 — Alkylation states: the audit, the defect, and the corrected baseline

2026-09-01. Trigger: the trace-forensics round found the S08 DIA-NN reference
had searched cysteine UNMODIFIED. This document is the project-wide audit
(verdicts BY MASS, never by id strings: 400 sampled cysteine precursors per
library tested against unmodified vs +57.02146/Cys), the root cause, and the
required parameters going forward.

## The audit

| library | used by | verdict |
|---|---|---|
| dn_pred_cam.parquet | samelib arms + corrected DIA-NN ref | **CAM** (400/400, 0 free) |
| mix10k_acq_lib.parquet | mix10k fixture arms | **CAM** (400/400) |
| human_v2.parquet | ODIA own-lib e2e (S08) | **CAM** (400/400) |
| astral_mix_lib.parquet | Astral arms (ODIA) | **CAM** (400/400) |
| corpus_lib.parquet | trace-model corpus era | **CAM** (400/400) |
| s08_targets_osw.tsv | OpenSWATH S08 | **CAM** (2,580 Carbamidomethyl annotations) |
| astral_osw.tsv | OpenSWATH Astral | **CAM-FREE — STALE**: 0 annotations in 130,692 rows; Mag-Net-era artifact (that sample had free thiols, doc/30); NOT comparable to the current neat sample |
| frozen DIA-NN S08 reference (33,330) | the benchmark denominator | **CAM-FREE — THE DEFECT** (400/400 cysteine matches at unmodified mass) |

Every ODIA-side arm ever measured on S08/mix/Astral searched the correct
chemistry. The defect is confined to the frozen S08 DIA-NN reference (and the
stale Mag-Net-era OSW Astral library, already superseded by policy).

## Where it happened, and why

Phase: the REFERENCE'S LIBRARY GENERATION, before this project's second
attempt began. `shared/ref/diann/how_s08_was_run.sh` searched
`full_lib.predicted.speclib`, a library predicted WITHOUT `--unimod4`.
DIA-NN's COMMAND LINE does not apply carbamidomethylation by default (the GUI
pre-checks it; the CLI does not), so the omission was silent. It was even
RECORDED — shared/ref/diann/README.md's own table says "searched CAM-free" —
but the number 33,330 was consumed everywhere as the reference regardless.
The Astral reference (re-cut 2026-08-17) used the correct recipe and is
unaffected.

Measured cost of the defect: the sample IS carbamidomethylated, so the
reference could not correctly match any of the ~24.8%-cysteine search space —
its accepted set is 1.9% cysteine (646, all at wrong chemistry) vs ODIA's
9.1%. All "recall of DIA-NN" numbers in docs 77-81 are measured against this
depressed, cysteine-blind denominator (their internal comparisons remain
valid — every arm shared the same reference — but the absolute recall levels
and the "36% missed" framing carry the caveat). 60.1% of the P1 port's
"complementary" IDs were simply outside the reference's search space.

## Required parameters (the record)

- **DIA-NN library generation**: `--fasta-search --predictor --gen-spec-lib
  --unimod4` (the Astral reference's recorded recipe). NEVER omit --unimod4
  on an alkylated sample; the CLI will not add it for you.
- **DIA-NN search**: the recorded flags (`--f <run.d> --lib <CAM lib>
  --fasta bench_reviewed_entrap.fasta --reannotate --threads 48
  --qvalue 0.01`); the .d folder, never the mzML (no MS2 spectra via mzML).
  DIA-NN 2.0 accepts the library as parquet.
- **ODIA / DIALibraryGenerator**: CAM is the shipped default — every audited
  ODIA library carries it. No action.
- **OpenSWATH**: S08 conversion (s08_targets_osw.tsv) is CAM-correct;
  astral_osw.tsv must be regenerated from a CAM source before any new Astral
  OSW number is quoted.
- **Samples**: S08/S30 (AGXT liver diaPASEF) and Astral NEAT are
  carbamidomethylated preparations. The retired Mag-Net Astral file was the
  free-thiol exception and is gone from the benchmark (2026-08-17 policy).

## The corrected baseline

The corrected S08 DIA-NN reference is being computed with the identical
recorded command and the CAM library (dn_pred_cam.parquet — the same library
ODIA's samelib arms search, making the corrected comparison same-library AND
same-chemistry). The corrected four-tool table lands in this document's
companion report when the run and the ODIA own-library P1 arm finish.
Until then: **the number 33,330 is SUPERSEDED for every purpose except
reproducing historical comparisons.**
