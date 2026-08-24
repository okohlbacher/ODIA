#!/usr/bin/env python3
"""Rungs 3 and 4: Deep Sets, then a transformer, on the same contrast and split.

Both must clear RUNG 1 = 0.8649 (the shipped 19 sub-scores), which is the bar
because rung 2's order statistics came in BELOW it at 0.8422. Codex: "the
transformer is not yet justified" until Deep Sets is measured.

Fragments are a SET -- their order in the library is arbitrary -- so every
architecture here is permutation-equivariant. A convolution over fragment index
would learn an ordering that does not exist.

Jitter augmentation (doc/56 A12): the window centre shifts randomly each epoch.
Training windows are centred on the PREDICTED retention time while the mandated
evaluation is on ODIA's picked apex, so a model keying on absolute position
inside the window would not transfer. A model that cannot use position cannot be
fooled by it either, which attacks A1's centring leak from the other side.

Representation (doc/56 A7): NOT per-group normalised -- normalising away
intensity removes exactly the low-abundance evidence under investigation.
Channels are log1p raw, background-subtracted and a nonzero mask, with
per-fragment scale supplied as a scalar.
"""
import sys, hashlib, numpy as np, pyarrow.parquet as pq, torch, torch.nn as nn
from sklearn.metrics import roc_auc_score

D = sys.argv[1] if len(sys.argv) > 1 else '/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer/shared/corpus'
ARCH = sys.argv[2] if len(sys.argv) > 2 else 'deepsets'
EPOCHS = int(sys.argv[3]) if len(sys.argv) > 3 else 12
# 'apex' centres each arm's crop on ODIA's OWN picked apex for that arm, which
# removes rung 1's advantage (it scores an already-located candidate) and is the
# deployment-distribution evaluation doc/56 A1 makes mandatory. No asymmetry:
# both arms get the same picker applied to their own window.
CENTRE = sys.argv[4] if len(sys.argv) > 4 else 'window'
dev = 'cuda' if torch.cuda.is_available() else 'cpu'
torch.manual_seed(20260824); np.random.seed(20260824)
print(f'device {dev}  arch {ARCH}')

lab = pq.read_table(f'{D}/tensor_s08_labels.parquet')
ids = np.array(lab.column('Precursor.Id').to_pylist())
dec = lab.column('Decoy').to_numpy()
label = np.array(lab.column('Label').to_pylist())
# Into RAM, NOT mmap. These live on /ceph, a network filesystem, and the
# training loop fancy-indexes a random batch every step -- mmap turns that into
# thousands of small remote reads and the GPU sits at 0% while it waits. 0.87 GB
# per arm against 2.3 TB of node memory.
X0 = np.ascontiguousarray(np.load(f'{D}/tensor_s08_traces.npy'))
X1 = np.ascontiguousarray(np.load(f'{D}/tensor_s08shift_traces.npy'))
M  = np.load(f'{D}/tensor_s08_mask.npy')

lib = pq.read_table(f'{D}/corpus_lib.parquet', columns=['Precursor.Id','Decoy','Protein.Group'])
pgm = {(p, int(d)): str(g) for p, d, g in zip(lib.column('Precursor.Id').to_pylist(),
                                              lib.column('Decoy').to_numpy(),
                                              lib.column('Protein.Group').to_pylist())}
fold = np.array([int(hashlib.md5(pgm.get((p, int(d)), p).encode()).hexdigest(), 16) % 10
                 for p, d in zip(ids, dec)])
pos = label == 'pos'
if CENTRE == 'apex':
    A0 = np.load(f'{D}/apex_s08.npy'); A1 = np.load(f'{D}/apex_s08shift.npy')
    # Both arms must have a located candidate, or the pair is not on the
    # deployment distribution and the two arms would be centred differently --
    # which is the very asymmetry this evaluation exists to rule out.
    have = (A0 >= 0) & (A1 >= 0)
    print(f'apex-centred: {have.sum():,} of {len(have):,} have a candidate in BOTH arms; '
          f'positives with both {int((pos & have).sum()):,} of {int(pos.sum()):,}')
    pos = pos & have
tr_idx = np.flatnonzero(pos & (fold >= 3))
te_idx = np.flatnonzero(pos & (fold < 3))
print(f'positives {pos.sum():,}  train {len(tr_idx):,}  test {len(te_idx):,} '
      f'(each contributes one positive and one shifted negative)')

C, W = 128, 96
def crop(t, centres, jitter):
    """Per-row crop, so each precursor is cut around its own apex."""
    lo = np.clip(centres - W // 2, 0, C - W)
    if jitter:
        lo = np.clip(lo + np.random.randint(-8, 9, size=len(lo)), 0, C - W)
    out = torch.empty((t.shape[0], t.shape[1], W), dtype=t.dtype)
    for k, a in enumerate(lo):
        out[k] = t[k, :, a:a + W]
    return out

def batch(idx, arm, jitter):
    T = (X0 if arm == 0 else X1)[idx]
    t = torch.from_numpy(np.asarray(T, dtype=np.float32))
    if CENTRE == 'apex':
        t = crop(t, (A0 if arm == 0 else A1)[idx].astype(int), jitter)
    else:
        o = np.random.randint(0, C - W + 1) if jitter else (C - W) // 2
        t = t[:, :, o:o + W]
    m = torch.from_numpy(M[idx].astype(np.float32))
    bg = t.median(dim=2, keepdim=True).values
    ch = torch.stack([torch.log1p(t), torch.clamp(t - bg, min=0).log1p(), (t > 0).float()], 2)
    return ch, m, torch.log1p(t.sum(2))

class FragEnc(nn.Module):
    """Shared 1-D encoder: every fragment goes through the SAME weights, which
    is what makes the set treatment permutation-equivariant."""
    def __init__(s, d=64):
        super().__init__()
        s.net = nn.Sequential(
            nn.Conv1d(3, 32, 7, padding=3), nn.GELU(), nn.MaxPool1d(2),
            nn.Conv1d(32, 64, 5, padding=2), nn.GELU(), nn.MaxPool1d(2),
            nn.Conv1d(64, d, 3, padding=1), nn.GELU(), nn.AdaptiveAvgPool1d(1))
        s.d = d
    def forward(s, ch, scale):
        B, F, K, L = ch.shape
        z = s.net(ch.reshape(B * F, K, L)).reshape(B, F, s.d)
        return torch.cat([z, scale.unsqueeze(-1)], -1)

class DeepSets(nn.Module):
    def __init__(s, d=64):
        super().__init__()
        s.enc = FragEnc(d)
        s.head = nn.Sequential(nn.Linear(3 * (d + 1) + 1, 128), nn.GELU(),
                               nn.Dropout(0.1), nn.Linear(128, 1))
    def forward(s, ch, m, scale):
        z = s.enc(ch, scale) * m.unsqueeze(-1)
        n = m.sum(1, keepdim=True).clamp(min=1)
        mean = z.sum(1) / n
        mx = z.masked_fill(m.unsqueeze(-1) == 0, -1e9).max(1).values
        return s.head(torch.cat([mean, mx, z.sum(1), n], -1)).squeeze(-1)

class SetTransformer(nn.Module):
    def __init__(s, d=64, heads=4, blocks=2):
        super().__init__()
        s.enc = FragEnc(d); s.proj = nn.Linear(d + 1, d)
        s.blocks = nn.ModuleList([nn.TransformerEncoderLayer(d, heads, 4 * d, 0.1,
                                  batch_first=True, norm_first=True) for _ in range(blocks)])
        s.q = nn.Parameter(torch.randn(1, 1, d))
        s.att = nn.MultiheadAttention(d, heads, batch_first=True)
        s.head = nn.Sequential(nn.Linear(d + 1, 128), nn.GELU(), nn.Linear(128, 1))
    def forward(s, ch, m, scale):
        z = s.proj(s.enc(ch, scale))
        pad = m == 0
        for b in s.blocks: z = b(z, src_key_padding_mask=pad)
        p, _ = s.att(s.q.expand(z.size(0), -1, -1), z, z, key_padding_mask=pad)
        return s.head(torch.cat([p.squeeze(1), m.sum(1, keepdim=True)], -1)).squeeze(-1)

model = (DeepSets() if ARCH == 'deepsets' else SetTransformer()).to(dev)
print(f'parameters {sum(p.numel() for p in model.parameters()):,}')
opt = torch.optim.AdamW(model.parameters(), lr=2e-3, weight_decay=1e-4)
lossf = nn.BCEWithLogitsLoss(); BS = 256

def run_epoch(idx, train, jitter):
    model.train(train)
    order = np.random.permutation(len(idx)) if train else np.arange(len(idx))
    sc, ys = [], []
    for a in range(0, len(order), BS):
        sel = np.sort(idx[order[a:a + BS]])
        outs, tgt = [], []
        for arm in (0, 1):
            ch, m, s_ = batch(sel, arm, jitter)
            o = model(ch.to(dev), m.to(dev), s_.to(dev))
            outs.append(o); tgt.append(torch.full_like(o, 1.0 - arm))
        o = torch.cat(outs); y = torch.cat(tgt)
        if train:
            opt.zero_grad(); lossf(o, y).backward(); opt.step()
        sc.append(o.detach().float().cpu().numpy()); ys.append(y.detach().cpu().numpy())
    return roc_auc_score(np.concatenate(ys), np.concatenate(sc)), np.concatenate(sc)

best = 0.0
for e in range(EPOCHS):
    a_tr, _ = run_epoch(tr_idx, True, True)
    with torch.no_grad():
        a_te, sc = run_epoch(te_idx, False, False)
    best = max(best, a_te)
    print(f'  epoch {e+1:>2}  train AUC {a_tr:.4f}   test AUC {a_te:.4f}')
print(f'\nRUNG {"3 Deep Sets" if ARCH=="deepsets" else "4 transformer"}  '
      f'[{CENTRE}-centred] final {a_te:.4f}  best {best:.4f}'
      f'   (bar: rung 1 = 0.8649, rung 2 = 0.8422)')
np.save(f'/tmp/scores_{ARCH}.npy', sc)
