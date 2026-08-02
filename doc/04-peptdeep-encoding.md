# PeptDeep input encoding, as ODIA must implement it

**Why this document exists.** OpenMS's PeptDeep bindings reject modified
peptides outright (`PeptDeepUtils.h:48` throws *"Modified peptides are not
currently supported in this engine"*, and `generateUnmodifiedModXTensor` at
`:59-61` returns all zeros). A library that cannot carry carbamidomethyl-C is
not usable, so ODIA has to build the `mod_x` tensor itself. Getting it wrong
mostly does not fail loudly — it produces plausible, quietly wrong predictions.
Every quantity below therefore carries its source, and where an error is silent
its measured cost is given.

**Sources.** AlphaPeptDeep upstream, since the models are its
`pretrained_models_v3/generic` weights converted by OpenMS's
`tools/scripts/export_peptdeep_models_to_onnx.py`. Local copies of
`peptdeep/model/featurize.py` and `peptdeep/settings.py` under
`/scratch/kohlbach/odia/` are byte-identical to `MannLabs/alphapeptdeep@main`.

---

## 1. The element list is 109 entries, and the order is the encoding

`mod_elements` in `model_const.yaml` (lines 3–113) has exactly **109** entries.
Order is load-bearing: `settings.py:27` builds
`mod_elem_to_idx = dict(zip(mod_elements, range(mod_feature_size)))`, which is
positional, so the yaml order *is* the feature order. Index 0 is `C`; index 108
is `?`.

```
C H N O P S            <- first six, the yaml comments "do not change them"
B F I K U V W X Y Ac Ag Al Am Ar As At Au Ba Be Bi Bk Br ...
... Zr 2H 13C 15N 18O ?
```

The independent confirmation of 109 is **the shipped models themselves** —
`peptdeep_{rt,ms2,ccs}_dynamic.onnx` declare `mod_x` as
`float32 [batch_size, seq_length, 109]`. OpenMS's `PEPTDEEP_MOD_ELEMENTS = 109`
is *not* independent evidence: it and the exporter's
`MOD_HIDDEN = len(model_const["mod_elements"])` both trace to this same yaml.

## 2. A modification becomes its elemental composition — signed, unnormalised

From `_parse_mod_formula` (`settings.py:30-43`):

```python
feature = np.zeros(109)
elems = formula.strip(")").split(")")     # "H(3)C(2)N(1)O(1)" -> ["H(3","C(2","N(1","O(1"]
for elem in elems:
    chem, num = elem.split("(")
    if chem in mod_elem_to_idx: feature[mod_elem_to_idx[chem]] = int(num)
    else:                       feature[-1] += int(num)     # '?' bucket
```

**Counts are signed, and negatives are common — not an isotope-label curiosity.**
600 of 2859 rows in alphabase's `modification.tsv` carry a negative count,
including everyday modifications:

```
Deamidated@N          H(-1)N(-1)O(1)
Dehydrated@S          H(-2)O(-1)
Amidated@Any_C-term   H(1)N(1)O(-1)
Label:13C(6)@K        C(-6)13C(6)
```

Storing counts in an unsigned type, or taking `abs()` when reading an OpenMS
`EmpiricalFormula`, is silently wrong. Measured on the RT model for
`SAMPLENPEPTIDEK` with Deamidated on N7: signed **0.5231**, clamped-to-zero
0.5190, `abs()` 0.5121 (unmodified 0.5293). Wrong, with no error.

**There is no scaling or normalisation of the mod feature.** The whole path was
audited — `_parse_mod_formula` → `MOD_TO_FEATURE` (`settings.py:49-54`) →
`get_batch_mod_feature` (`featurize.py:50-79`) →
`ModelInterface._get_mod_features` (`model_interface.py:930-937`) → a plain
tensor cast. Dividing a carbamidomethyl vector by 10 moves RT from 0.5569 to
0.6089.

## 3. Placement is by site, on an (nAA + 2) axis

From `parse_mod_feature` (`featurize.py:43,46`): the tensor is
`(nAA + 2, 109)` per peptide and modifications are **added** at their site.

| site | index | meaning |
|---|---|---|
| `0` | `0` | peptide N-terminal modification |
| `1 … nAA` | same | residue modification (1-based) |
| `-1` | `nAA + 1` | peptide C-terminal modification |

`+=`, not `=`, so co-located modifications accumulate. Placement is
load-bearing: acetyl at row 0 gives RT 0.7551, the same vector at row 1 gives
0.6910.

## 4. Inputs — shapes, dtypes and scale factors

Confirmed against the shipped `.onnx` files, not only against the exporter.

| input | shape | dtype | scale | models |
|---|---|---|---|---|
| `input_sequences` / `aa_indices` | `[batch, seq_len]` | int64 | — | all three |
| `mod_x` | `[batch, seq_len, 109]` | float32 | none | all three |
| `charges` | `[batch, 1]` | float32 | **× 0.1** | MS2, CCS only |
| `nce` | `[batch, 1]` | float32 | **× 0.01** | MS2 only |
| `instrument_indices` | **`[batch]`** | int64 | — | MS2 only |

Three traps here, in descending order of damage:

- **NCE must be multiplied by 0.01.** `ms2.py:423-424`, applied at `:494`;
  OpenMS agrees (`PeptDeepInput.cpp:19` `NCE_SCALE = 0.01f`). Passing raw NCE is
  not a degradation — for `ELVISLIVESK`, charge 2, NCE 30, the cosine similarity
  between the correct spectrum and the unscaled one is **0.0028**. Unrelated
  output, no error.
- **Charge must be multiplied by 0.1.** `ms2.py:423` and `ccs.py:129`
  `charge_factor = 0.1`; OpenMS `PeptDeepInput.cpp:18` `CHARGE_SCALE = 0.1f`.
  The two agree. Cosine between `charges=0.2` and `charges=2.0` is 0.6377.
- **`instrument_indices` is rank 1**, unlike every other meta input.
  `export_peptdeep_models_to_onnx.py:61` is `torch.zeros((batch_size,))`, and
  `PeptDeepMS2Inference.cpp:137` carries the comment `// <--- KEEP AS 1D`
  against `// <--- REVERT TO 2D` for charges and nce. Feeding `[batch, 1]` is
  rejected: *"Invalid rank for input: instrument_indices Got: 2 Expected: 1"*.

Instrument mapping (`featurize.py:135-147`, `model_const.yaml`): `QE=0`,
`Lumos=1`, `timsTOF=2`, `SciexTOF=3`, `ThermoTOF=4`; anything unknown maps to
`max_instrument_num - 1 = 7`.

`seq_len` is `nAA + 2`, the terminal tokens being index 0 and `nAA+1`.
`aa_embedding_size: 27` corroborates A→1…Z→26 plus the 0 pad, which is one-hot
encoded at `building_block.py:29`.

## 5. Every peptide in a batch must have the same length

`featurize.py:57` states it: *"All sequence lengths must be the same, meaning
that nAA values must be equal."*

This is not a convention — **trailing padding is not inert**. Index 0 is one-hot
encoded rather than a zero embedding, and none of `SeqCNN` (`:114-139`),
`SeqLSTM` (bidirectional, `:335-364`) or `SeqAttentionSum` (`:399-413`) applies
a padding mask. The same peptide, `ELVISLIVESK`, with only the padding changed:

| `seq_len` | predicted RT |
|---|---|
| 13 (correct) | 0.8196 |
| 15 | 0.7022 |
| 20 | 0.4737 |
| 30 | 0.2714 |

So ODIA must group peptides by encoded length and run one batch per length, as
OpenMS does (`PeptDeepRTInference.cpp:56-64`, likewise CCS and MS2). Note that
because OpenMS always groups, its own padding branch (`PeptDeepInput.cpp:86-90`)
is dead code and has never been exercised against anything.

## 6. Outputs

| model | output shape | meaning |
|---|---|---|
| RT | `[batch, 1]` | normalised iRT |
| CCS | `[batch, 1]` | CCS |
| MS2 | `[batch, seq_len - 3, 8]` | `nAA - 1` fragment positions × 8 channels |

MS2 returns `out_x[:, 3:, :]` (`ms2.py:270`); for an 11-mer that is `(1, 10, 8)`.
Channel order comes from alphabase's `sort_charged_frag_types`, which is
`sorted(no_loss) + sorted(loss)` over `frag_types: [b, y, b_modloss, y_modloss]`
× `max_frag_charge: 2`:

```
0 b_z1   1 b_z2   2 y_z1   3 y_z2
4 b_modloss_z1   5 b_modloss_z2   6 y_modloss_z1   7 y_modloss_z2
```

## 7. Getting compositions without adding a dependency

AlphaPeptDeep reads compositions from alphabase's modification table. ODIA does
not need it: OpenMS's `ModificationsDB` already carries UniMod, and
`ResidueModification::getDiffFormula()` gives an `EmpiricalFormula` that
decomposes to (element symbol, signed count).

**The one translation needed is isotopes.** OpenMS builds isotope symbols as
`"(" + mass_number + ")" + symbol` (`ElementDB.cpp:640-641`;
`UnimodXMLHandler.cpp:142-148` converts UniMod's `13C` to `(13)C`), whereas the
element list spells them `2H`, `13C`, `15N`, `18O`. The complete set of element
symbols appearing in the shipped `unimod.xml` is

```
13C 15N 18O 2H Ag Al As Au B Br C Ca Cd Cl Co Cr Cu F Fe H Hg I K Li Mg Mn Mo
N Na Ni O P Pd Pt Ru S Se Si Zn
```

all 40 of which are in the 109-element list, so those four isotope renames are
the *only* mapping required. Without them, isotope-labelled modifications land
silently in the `?` bucket.

For completeness: across all 2859 rows of the shipped table, no composition
repeats an element token and none uses a symbol outside the list — so
`_parse_mod_formula`'s last-write-wins behaviour for known elements, and the `?`
bucket itself, are unreachable from UniMod data. They become reachable only
through a faulty symbol translation, which is the hazard above.

## 8. What this can and cannot be validated against

- **Unmodified peptides: ODIA can be checked against OpenMS.** Both should
  produce identical tensors and identical predictions from the same model. This
  is worth doing, but it is agreement with **a sibling reading of these same
  sources**, not with AlphaPeptDeep — every constant OpenMS uses was copied from
  the files above, so a shared misreading passes cleanly. It also cannot cover
  the padding path, which OpenMS never exercises (§5).
- **Modified peptides: not validatable here at all.** Nothing on this machine
  can say what AlphaPeptDeep would predict for a modified peptide.

**Therefore: no library generated with modifications should be trusted until the
modified path is checked against AlphaPeptDeep reference values.** That is a
blocking backlog item. A silently wrong intensity pattern degrades
identifications in a way that looks like an algorithmic shortfall rather than a
bug — exactly the failure mode both hand-off documents warn about.

---

*Revision note: §§1–3 and §7 survived adversarial review unchanged in substance.
§4 previously gave `instrument_indices` the wrong rank, omitted both scale
factors, and §5 did not exist — an implementation built on that version would
have produced spectra with cosine 0.0028 against the truth. §§1 and 8 previously
claimed an independence and a validation strength neither had.*
