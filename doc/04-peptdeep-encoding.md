# PeptDeep input encoding, as ODIA must implement it

**Why this document exists.** OpenMS's PeptDeep bindings reject modified
peptides outright (`PeptDeepUtils::validatePeptide` throws *"Modified peptides
are not currently supported in this engine"*, and `generateUnmodifiedModXTensor`
returns a zero tensor). A library that cannot carry carbamidomethyl-C is not
usable, so ODIA has to build the `mod_x` tensor itself. Getting it wrong does
not fail loudly — it produces plausible, quietly wrong predictions. This records
the encoding and where each part was established.

**Sources.** All from AlphaPeptDeep upstream, since the models are its
`pretrained_models_v3/generic` weights converted by OpenMS's
`tools/scripts/export_peptdeep_models_to_onnx.py`:

| item | source |
|---|---|
| element list and order | `peptdeep/constants/model_const.yaml`, key `mod_elements` |
| composition → vector | `peptdeep/settings.py`, `_parse_mod_formula` |
| placement by site | `peptdeep/model/featurize.py`, `parse_mod_feature` |
| tensor shapes | `export_peptdeep_models_to_onnx.py` and `PeptDeepMS2Inference.cpp` |

---

## 1. The element list is 109 entries, and the count corroborates itself

`mod_elements` in `model_const.yaml` has exactly **109** entries, which matches
`PEPTDEEP_MOD_ELEMENTS = 109` in OpenMS's `PeptDeepUtils.h` independently. The
order is load-bearing — the index *is* the feature position:

```
C H N O P S            <- first six, the yaml comments "do not change them"
B F I K U V W X Y Ac Ag Al Am Ar As At Au Ba Be Bi Bk Br ...
... Zr 2H 13C 15N 18O ?
```

The final entry `?` is the catch-all bucket.

## 2. A modification becomes its elemental composition, unnormalised

From `_parse_mod_formula`:

```python
feature = np.zeros(109)
elems = formula.strip(")").split(")")     # "H(3)C(2)N(1)O(1)" -> ["H(3","C(2","N(1","O(1"]
for elem in elems:
    chem, num = elem.split("(")
    if chem in mod_elem_to_idx: feature[mod_elem_to_idx[chem]] = int(num)
    else:                       feature[-1] += int(num)     # '?' bucket
```

Two things to note, both easy to get wrong:

- **There is no scaling or normalisation.** The feature is raw element counts.
  Dividing by anything would be wrong.
- **Unknown elements accumulate into the last slot** (`+=`), while known
  elements are *assigned* (`=`). A formula repeating a known element would keep
  only the last count — faithful reproduction means matching that, not
  "fixing" it.

## 3. Placement is by site, on an (nAA + 2) axis

From `parse_mod_feature`: the tensor is `(nAA + 2, 109)` per peptide and
modifications are **added** at their site:

| site | index | meaning |
|---|---|---|
| `0` | `0` | peptide N-terminal modification |
| `1 … nAA` | same | residue modification (1-based) |
| `-1` | `nAA + 1` | peptide C-terminal modification |

`+=`, not `=`, so two modifications on one site accumulate. The `nAA + 2` axis
is the same padding the amino-acid tensor uses — index 0 and index `nAA+1` are
the terminal tokens, which is why `PeptDeepInputConfig::add_terminal_tokens`
writes a 0 at each end.

## 4. Tensor shapes

| input | shape | type | notes |
|---|---|---|---|
| `aa_indices` / `input_sequences` | `[batch, seq_len]` | int64 | `A→1 … Z→26`, 0 = pad/terminal |
| `mod_x` | `[batch, seq_len, 109]` | float32 | as above |
| `charges` | `[batch, 1]` | float32 | scaled by PeptDeep's charge factor |
| `nce` | `[batch, 1]` | float32 | MS2 model only |
| `instrument_indices` | `[batch, 1]` | int64 | MS2 model only; `max_instrument_num: 8` |

`seq_len` is `nAA + 2` with terminal tokens. `aa_embedding_size: 27` in
`model_const.yaml` corroborates the 1–26 residue encoding plus the zero pad.

## 5. Getting compositions without adding a dependency

AlphaPeptDeep reads compositions from alphabase's modification table. ODIA does
not need it: OpenMS's `ModificationsDB` already carries UniMod, and
`ResidueModification::getDiffFormula()` gives an `EmpiricalFormula` that
decomposes to (element symbol, count) — which is precisely what
`_parse_mod_formula` produces from its string. Using OpenMS keeps the
dependency set unchanged and uses the same UniMod data the rest of the tool
does.

**One mapping hazard:** the element list contains the isotopes `2H`, `13C`,
`15N`, `18O`, whereas OpenMS writes isotopes as `(13)C`. That translation has to
be explicit, or isotope-labelled modifications land silently in the `?` bucket.

## 6. What this can and cannot be validated against

- **Unmodified peptides: fully validatable.** ODIA's encoder and OpenMS's
  `PeptDeepInputBuilder` must produce identical tensors and therefore identical
  predictions from the same model. Any disagreement is an ODIA bug. This covers
  the residue encoding, the terminal-token padding, the charge scaling and the
  output decoding — everything except the modification path itself.
- **Modified peptides: not validatable here.** Nothing on this machine can say
  what AlphaPeptDeep would predict for a modified peptide. Agreement with the
  spec above is an argument, not evidence.

**Therefore: no library generated with modifications should be trusted until the
modified path is checked against AlphaPeptDeep reference values.** That is a
blocking backlog item, not a nicety — a silently wrong intensity pattern would
degrade identifications in a way that looks like an algorithmic shortfall rather
than a bug, which is exactly the failure mode both hand-off documents warn about.
