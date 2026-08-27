"""Canonical identifiers for joining ODIA, DIA-NN and OpenSWATH results.

Every silent-failure mode this project has hit at a join boundary is encoded
here once, so that no script has to remember it:

  * The modification alias. DIA-NN writes `C(UniMod:4)`; ODIA's library writes
    `C(Carbamidomethyl)`. A verbatim join drops 3,902 of DIA-NN's 39,149
    confident S08 precursors -- 10.0% -- and reports nothing. The same alias was
    copy-pasted into 25 scripts under 4 different variable names. Measured
    2026-08-27 on the XIC: 89.74% overlap verbatim, 100.000% aliased.

  * Decoy disambiguation. ODIA's `-out` carries no `_decoy` suffix on
    Precursor.Id while `-out_lib` does, so a name-only join sums a decoy into
    its target's trace. Join on (Precursor.Id, Decoy), always.

  * Transition identity. Joining on product m/z rounded to 3 decimals collides
    across ion series and is not stable across calibration. Use the canonical
    tuple (type, ordinal, charge, loss).

  * Acquisition-cell identity. Rounded retention time is not a cell key in
    diaPASEF: one RT can carry several mobility windows. Use
    (frame/cycle, window) where available and say so when it is not.

Nothing here rounds, fills or drops silently. Functions that cannot answer
raise instead of guessing.
"""
from __future__ import annotations

# ODIA's Fragment.Type codes. Verified against shared/libv2/coh_matched.py:4 and
# the library parquet, which contains only these two values.
FRAGMENT_TYPE_CODE = {2: 'b', 5: 'y'}
FRAGMENT_CODE_BY_LETTER = {v: k for k, v in FRAGMENT_TYPE_CODE.items()}

# The only alias that has actually bitten. Kept as a table so a second one can
# be added without another 25-script sweep.
MOD_ALIASES = (
    ('(Carbamidomethyl)', '(UniMod:4)'),
)


def to_diann_dialect(s: str) -> str:
    """ODIA spelling -> DIA-NN spelling. Safe on already-DIA-NN input."""
    for odia, diann in MOD_ALIASES:
        s = s.replace(odia, diann)
    return s


def to_odia_dialect(s: str) -> str:
    """DIA-NN spelling -> ODIA spelling. Safe on already-ODIA input."""
    for odia, diann in MOD_ALIASES:
        s = s.replace(diann, odia)
    return s


def canonical_precursor(precursor_id: str, decoy=0) -> tuple:
    """The join key for a precursor: dialect-normalised id plus decoy status.

    `decoy` is taken from the Decoy column, never inferred from a name suffix --
    a `_decoy` suffix is present in -out_lib and absent in -out for the same row.
    A trailing `_decoy` is stripped so the two spellings agree, but it also sets
    the flag when no column was supplied.
    """
    pid = to_diann_dialect(precursor_id)
    if pid.endswith('_decoy'):
        pid = pid[:-len('_decoy')]
        if decoy is None:
            decoy = 1
    if decoy is None:
        raise ValueError(
            f'decoy status unknown for {precursor_id!r}: pass the Decoy column. '
            'A name-only join merges a decoy into its target.')
    return (pid, int(bool(decoy)))


def canonical_transition(frag_type, ordinal: int, charge: int, loss: str = 'noloss') -> str:
    """Stable transition id, e.g. 'y7^1' or 'b4^2-H2O'.

    `frag_type` accepts either ODIA's numeric code or a letter. Product m/z is
    deliberately NOT part of the key: it moves under calibration and collides
    across series at 3 decimals.
    """
    if isinstance(frag_type, str):
        letter = frag_type
    else:
        letter = FRAGMENT_TYPE_CODE.get(frag_type)
        if letter is None:
            raise ValueError(f'unknown ODIA Fragment.Type code {frag_type!r}; '
                             f'known: {sorted(FRAGMENT_TYPE_CODE)}')
    tag = f'{letter}{int(ordinal)}^{int(charge)}'
    if loss and loss != 'noloss':
        tag += f'-{loss}'
    return tag


def canonical_cell(frame=None, window=None, rt=None) -> tuple:
    """Acquisition-cell key. Prefers (frame, window); refuses to fake one.

    Retention time alone is not a cell identity in diaPASEF. Callers that only
    have RT get a key tagged as approximate, so a downstream assertion can tell
    the difference instead of silently comparing incomparable cells.
    """
    if frame is not None and window is not None:
        return ('cell', int(frame), int(window))
    if rt is None:
        raise ValueError('need (frame, window) or at least rt')
    return ('rt~', float(rt))


def assert_schema(table, expected_columns, name='table', rows=None, sha256=None):
    """Fail loudly on the artefact swaps that have silently changed results.

    Two files in shared/libv2 have IDENTICAL byte counts (92,228,109) and
    different content -- sub_apex.parquet and sub_apex_seeded.parquet -- and two
    scripts behind the same document each read a different one. A size check
    cannot catch that; a hash can.
    """
    have = list(getattr(table, 'column_names', None) or table.schema.names)
    missing = [c for c in expected_columns if c not in have]
    if missing:
        raise AssertionError(f'{name}: missing columns {missing}; has {have}')
    if rows is not None and table.num_rows != rows:
        raise AssertionError(f'{name}: expected {rows:,} rows, got {table.num_rows:,}')
    if sha256 is not None:
        raise NotImplementedError(
            f'{name}: hash the file on disk with file_sha256(), not the table')
    return True


def file_sha256(path, _bufsize=1 << 20) -> str:
    import hashlib
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(_bufsize), b''):
            h.update(chunk)
    return h.hexdigest()
