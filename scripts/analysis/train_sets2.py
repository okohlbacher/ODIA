#!/usr/bin/env python3
"""The information-matched rerun (doc/59), with three defects fixed.

1. LIBRARY PRIOR. All three reviewers found that the trace models never saw
   which fragments SHOULD be strong, while rung 1's var_library_corr /
   _dotprod / _rmsd compare observed against expected. `--feat` ablates it:
       anon  traces only, as doc/58 ran
       lib   + per-fragment library relative intensity
       full  + product m/z, ion series, ordinal, fragment charge
   The ladder anon -> lib -> full isolates exactly the signal the reviewers say
   was missing, and it is the var_library_corr signal that collapses hardest at
   low abundance.

2. TEMPORAL BOTTLENECK. doc/58's encoder ended in AdaptiveAvgPool1d(1), which
   averages each fragment's cycles away BEFORE any cross-fragment attention --
   kimi's objection, and the literature agrees: Alpha-XIC (the only trace-level
   model in doc/48 whose claim survived review, +9.4-16.2% appended to DIA-NN)
   uses Bi-GRU x2 with self-attention and keeps the time axis. Replaced with a
   bidirectional GRU over the convolved sequence.

3. EPOCH SELECTION ON TEST. doc/58 evaluated test every epoch and reported the
   final. Now a three-way PROTEIN split: folds 0-2 test, 3-4 validation, 5-9
   train. The epoch is chosen on validation; test is read once, at that epoch.
   Multiple seeds give a spread instead of one number, because 0.0004 was never
   interpretable without one.
"""
import sys, hashlib, numpy as np, pyarrow.parquet as pq, torch, torch.nn as nn
from sklearn.metrics import roc_auc_score

D      = sys.argv[1] if len(sys.argv) > 1 else '/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer/shared/corpus'
FEAT   = sys.argv[2] if len(sys.argv) > 2 else 'lib'
EPOCHS = int(sys.argv[3]) if len(sys.argv) > 3 else 30
CENTRE = sys.argv[4] if len(sys.argv) > 4 else 'apex'
SEEDS  = [int(x) for x in (sys.argv[5].split(',') if len(sys.argv) > 5 else ['0', '1', '2'])]
# 'shift'  positive = own window, negative = THE SAME precursor +300 s.
#          Metadata is arm-invariant, so leakage is impossible by construction.
#          Asks: is the peptide HERE rather than there. Localisation-flavoured.
# 'ent'    positive = DIA-NN-confident target, negative = ENTRAPMENT target,
#          both in their own window. Asks the PRESENCE question the project
#          actually cares about -- at the cost that the two are different
#          precursors of different organisms, so metadata CAN separate them and
#          the leakage floor has to be measured rather than assumed.
CONTRAST = sys.argv[6] if len(sys.argv) > 6 else 'shift'
dev = 'cuda' if torch.cuda.is_available() else 'cpu'
print(f'device {dev}  feat {FEAT}  centre {CENTRE}  epochs {EPOCHS}  seeds {SEEDS}')

lab = pq.read_table(f'{D}/tensor_s08_labels.parquet')
ids = np.array(lab.column('Precursor.Id').to_pylist())
dec = lab.column('Decoy').to_numpy()
label = np.array(lab.column('Label').to_pylist())
m2 = pq.read_table(f'{D}/tensor_s08shift_meta.parquet')
assert (np.array(m2.column('Precursor.Id').to_pylist()) == ids).all() and \
       (m2.column('Decoy').to_numpy() == dec).all(), 'arms not row-aligned'

X0 = np.ascontiguousarray(np.load(f'{D}/tensor_s08_traces.npy'))
X1 = np.ascontiguousarray(np.load(f'{D}/tensor_s08shift_traces.npy'))
M  = np.load(f'{D}/tensor_s08_mask.npy')
RI = np.load(f'{D}/desc_relint.npy')
PM = np.load(f'{D}/desc_prodmz.npy') / 1000.0
FC = np.load(f'{D}/desc_frcharge.npy')
OR = np.load(f'{D}/desc_ordinal.npy') / 20.0
SE = np.load(f'{D}/desc_series.npy')
A0 = np.load(f'{D}/apex_s08.npy'); A1 = np.load(f'{D}/apex_s08shift.npy')

lib = pq.read_table(f'{D}/corpus_lib.parquet', columns=['Precursor.Id','Decoy','Protein.Group'])
pgm = {(p, int(d)): str(g) for p, d, g in zip(lib.column('Precursor.Id').to_pylist(),
                                              lib.column('Decoy').to_numpy(),
                                              lib.column('Protein.Group').to_pylist())}
fold = np.array([int(hashlib.md5(pgm.get((p, int(d)), p).encode()).hexdigest(), 16) % 10
                 for p, d in zip(ids, dec)])
pos = label == 'pos'
neg = label == 'ent'
if CENTRE == 'apex':
    ok = (A0 >= 0) if CONTRAST == 'ent' else ((A0 >= 0) & (A1 >= 0))
    pos = pos & ok
    neg = neg & ok
if CONTRAST == 'shift':
    sel = pos
    print(f'contrast SHIFT: {pos.sum():,} precursors, each its own positive and negative')
else:
    sel = pos | neg
    print(f'contrast ENT: {pos.sum():,} positives against {neg.sum():,} entrapment '
          f'negatives -- DIFFERENT precursors, so the metadata floor must be measured')
tr_i = np.flatnonzero(sel & (fold >= 5))
va_i = np.flatnonzero(sel & (fold >= 3) & (fold < 5))
te_i = np.flatnonzero(sel & (fold < 3))
ispos = pos.astype(np.float32)
print(f'train {len(tr_i):,}  val {len(va_i):,}  test {len(te_i):,}')

C, W = 128, 96
NDESC = {'anon': 1, 'lib': 2, 'full': 6}[FEAT]

def batch(idx, arm, jitter):
    t = torch.from_numpy(np.asarray((X0 if arm == 0 else X1)[idx], dtype=np.float32))
    if CENTRE == 'apex':
        cen = (A0 if arm == 0 else A1)[idx].astype(int)
        lo = np.clip(cen - W // 2, 0, C - W)
        if jitter: lo = np.clip(lo + np.random.randint(-8, 9, size=len(lo)), 0, C - W)
        t = torch.stack([t[k, :, a:a + W] for k, a in enumerate(lo)])
    else:
        o = np.random.randint(0, C - W + 1) if jitter else (C - W) // 2
        t = t[:, :, o:o + W]
    m = torch.from_numpy(M[idx].astype(np.float32))
    bg = t.median(dim=2, keepdim=True).values
    ch = torch.stack([torch.log1p(t), torch.clamp(t - bg, min=0).log1p(), (t > 0).float()], 2)
    d = [torch.log1p(t.sum(2))]                               # observed scale
    if FEAT in ('lib', 'full'):
        d.append(torch.from_numpy(RI[idx]))                   # the EXPECTED pattern
    if FEAT == 'full':
        d += [torch.from_numpy(PM[idx]), torch.from_numpy(FC[idx]),
              torch.from_numpy(OR[idx]), torch.from_numpy(SE[idx])]
    return ch, m, torch.stack(d, -1)

class FragEnc(nn.Module):
    """Conv front end then a bidirectional GRU over time -- the time axis is
    kept, not averaged away. Alpha-XIC's shape, and the fix for doc/58's
    AdaptiveAvgPool1d bottleneck."""
    def __init__(s, d=64):
        super().__init__()
        s.conv = nn.Sequential(nn.Conv1d(3, 32, 7, padding=3), nn.GELU(), nn.MaxPool1d(2),
                               nn.Conv1d(32, 48, 5, padding=2), nn.GELU(), nn.MaxPool1d(2))
        s.gru = nn.GRU(48, d // 2, num_layers=1, batch_first=True, bidirectional=True)
        s.d = d
    def forward(s, ch, desc):
        B, F, K, L = ch.shape
        z = s.conv(ch.reshape(B * F, K, L)).transpose(1, 2)
        o, _ = s.gru(z)
        z = o[:, -1, :].reshape(B, F, s.d)
        return torch.cat([z, desc], -1)

class SetTransformer(nn.Module):
    def __init__(s, d=64, heads=4, blocks=2, ndesc=1):
        super().__init__()
        s.enc = FragEnc(d); s.proj = nn.Linear(d + ndesc, d)
        s.blocks = nn.ModuleList([nn.TransformerEncoderLayer(d, heads, 4 * d, 0.1,
                                  batch_first=True, norm_first=True) for _ in range(blocks)])
        s.q = nn.Parameter(torch.randn(1, 1, d))
        s.att = nn.MultiheadAttention(d, heads, batch_first=True)
        s.head = nn.Sequential(nn.Linear(d + 1, 128), nn.GELU(), nn.Linear(128, 1))
    def forward(s, ch, m, desc):
        z = s.proj(s.enc(ch, desc)); pad = m == 0
        for b in s.blocks: z = b(z, src_key_padding_mask=pad)
        p, _ = s.att(s.q.expand(z.size(0), -1, -1), z, z, key_padding_mask=pad)
        return s.head(torch.cat([p.squeeze(1), m.sum(1, keepdim=True)], -1)).squeeze(-1)

def evaluate(model, idx):
    model.eval(); sc, ys = [], []
    with torch.no_grad():
        for a in range(0, len(idx), 512):
            s_ = idx[a:a + 512]
            if CONTRAST == 'ent':
                ch, m, d = batch(s_, 0, False)
                o = model(ch.to(dev), m.to(dev), d.to(dev))
                sc.append(o.float().cpu().numpy()); ys.append(ispos[s_])
            else:
                for arm in (0, 1):
                    ch, m, d = batch(s_, arm, False)
                    o = model(ch.to(dev), m.to(dev), d.to(dev))
                    sc.append(o.float().cpu().numpy()); ys.append(np.full(len(s_), 1.0 - arm))
    return roc_auc_score(np.concatenate(ys), np.concatenate(sc))

results = []
for seed in SEEDS:
    torch.manual_seed(seed); np.random.seed(seed)
    model = SetTransformer(ndesc=NDESC).to(dev)
    opt = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, EPOCHS)
    lossf = nn.BCEWithLogitsLoss()
    best_va, best_te, best_ep = 0.0, 0.0, 0
    for e in range(EPOCHS):
        model.train()
        order = np.random.permutation(len(tr_i))
        for a in range(0, len(order), 256):
            s_ = np.sort(tr_i[order[a:a + 256]])
            outs, tgt = [], []
            if CONTRAST == 'ent':
                ch, m, d = batch(s_, 0, True)
                outs.append(model(ch.to(dev), m.to(dev), d.to(dev)))
                tgt.append(torch.from_numpy(ispos[s_]).to(dev))
            else:
                for arm in (0, 1):
                    ch, m, d = batch(s_, arm, True)
                    o = model(ch.to(dev), m.to(dev), d.to(dev))
                    outs.append(o); tgt.append(torch.full_like(o, 1.0 - arm))
            opt.zero_grad(); lossf(torch.cat(outs), torch.cat(tgt)).backward(); opt.step()
        sched.step()
        va = evaluate(model, va_i)
        # Test is read at the validation-selected epoch only. It is never used
        # to choose anything.
        if va > best_va:
            best_va, best_ep = va, e + 1
            best_te = evaluate(model, te_i)
    print(f'  seed {seed}: val {best_va:.4f} at epoch {best_ep}  ->  TEST {best_te:.4f}')
    results.append(best_te)

r = np.array(results)
print(f'\n{FEAT.upper():>5} [{CENTRE}/{CONTRAST}]  test AUC {r.mean():.4f} +- {r.std():.4f} '
      f'over {len(r)} seeds   (rung 1 = 0.8649)')
