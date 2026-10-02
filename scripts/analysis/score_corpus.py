#!/usr/bin/env python3
"""Rank the whole corpus by the trace model, and by ODIA's shipped DScore, and
compare identifications at MATCHED ENTRAPMENT FDP.

AUC on a paired contrast is not identifications. This project's currency is
precursors at matched entrapment FDP -- (e/r)/t over ordinary targets, doc/43's
definition -- and nothing measured so far is in it.

The model is trained on positive-vs-its-own-shifted-window and never sees an
entrapment precursor, so doc/56 A4's warning (a model trained AGAINST entrapment
suppresses it in production and makes entrapment FDP read optimistic) should not
apply. That is a claim to EARN, not assert, so the entrapment rate among
top-ranked precursors is reported for both rankings.

Everything is scored on the OWN-window arm only; the shifted arm exists to train
the contrast and has no place in a ranking of real precursors.
"""
import sys, numpy as np, pyarrow.parquet as pq, torch, torch.nn as nn, hashlib

D    = sys.argv[1]
CKPT = sys.argv[2]
dev  = 'cuda' if torch.cuda.is_available() else 'cpu'
ck   = torch.load(CKPT, map_location=dev, weights_only=False)
FEAT, NDESC, NSCAL = ck['feat'], ck['ndesc'], ck['nscal']
print(f"checkpoint: feat={FEAT} contrast={ck['contrast']} seed={ck['seed']} "
      f"val={ck['val']:.4f} epoch={ck['epoch']}")

lab = pq.read_table(f'{D}/tensor_ih1_labels.parquet')
ids = np.array(lab.column('Precursor.Id').to_pylist())
dec = lab.column('Decoy').to_numpy(); label = np.array(lab.column('Label').to_pylist())
A0  = np.load(f'{D}/apex_ih1.npy')
X0  = np.ascontiguousarray(np.load(f'{D}/tensor_ih1_traces.npy'))
M   = np.load(f'{D}/tensor_ih1_mask.npy')
RI  = np.load(f'{D}/desc_relint.npy'); PM = np.load(f'{D}/desc_prodmz.npy') / 1000.0
FC  = np.load(f'{D}/desc_frcharge.npy'); OR = np.load(f'{D}/desc_ordinal.npy') / 20.0
SE  = np.load(f'{D}/desc_series.npy');  PD = np.load(f'{D}/desc_precursor.npy')
SC0 = np.load(f'{D}/scal_ih1.npy'); OK0 = np.load(f'{D}/scal_ih1_ok.npy')
DS  = np.load(f'{D}/dscore_ih1.npy'); DSOK = np.load(f'{D}/dscore_ih1_ok.npy')

C, W = 128, 96
GAUSS = np.exp(-0.5 * ((np.arange(W) - (W - 1) / 2.0) / 1.07) ** 2).astype(np.float32)
lib = pq.read_table(f'{D}/corpus_lib.parquet', columns=['Precursor.Id','Decoy','Protein.Group'])
pgm = {(p, int(d)): str(g) for p, d, g in zip(lib.column('Precursor.Id').to_pylist(),
        lib.column('Decoy').to_numpy(), lib.column('Protein.Group').to_pylist())}
fold = np.array([int(hashlib.md5(pgm.get((p, int(d)), p).encode()).hexdigest(), 16) % 10
                 for p, d in zip(ids, dec)])

# Scalars were standardised on TRAIN rows in training; reproduce exactly.
trmask = (label == 'pos') & (A0 >= 0) & (fold >= 5)
_st = SC0[trmask]; SC_MU = _st.mean(0).astype(np.float32)
SC_SD = np.maximum(_st.std(0), 1e-6).astype(np.float32)

exec(open(__file__.replace('score_corpus.py', 'train_sets2_models.py')).read()) \
    if False else None

class FragEnc(nn.Module):
    def __init__(s, d=64):
        super().__init__()
        s.conv = nn.Sequential(nn.Conv1d(4, 32, 7, padding=3), nn.GELU(), nn.MaxPool1d(2),
                               nn.Conv1d(32, 48, 5, padding=2), nn.GELU(), nn.MaxPool1d(2))
        s.gru = nn.GRU(48, d // 2, num_layers=1, batch_first=True, bidirectional=True)
        s.d = d
    def forward(s, ch, desc):
        B, F, K, L = ch.shape
        z = s.conv(ch.reshape(B * F, K, L)).transpose(1, 2)
        o, _ = s.gru(z)
        return torch.cat([o[:, -1, :].reshape(B, F, s.d), desc], -1)

class SetTransformer(nn.Module):
    def __init__(s, d=64, heads=4, blocks=2, ndesc=1, nscal=0):
        super().__init__()
        s.enc = FragEnc(d); s.proj = nn.Linear(d + ndesc, d)
        s.blocks = nn.ModuleList([nn.TransformerEncoderLayer(d, heads, 4 * d, 0.1,
                                  batch_first=True, norm_first=True) for _ in range(blocks)])
        s.q = nn.Parameter(torch.randn(1, 1, d))
        s.att = nn.MultiheadAttention(d, heads, batch_first=True)
        s.head = nn.Sequential(nn.Linear(d + 1 + nscal, 128), nn.GELU(),
                               nn.Dropout(0.1), nn.Linear(128, 1))
    def forward(s, ch, m, desc, scal=None):
        z = s.proj(s.enc(ch, desc)); pad = m == 0
        for b in s.blocks: z = b(z, src_key_padding_mask=pad)
        p, _ = s.att(s.q.expand(z.size(0), -1, -1), z, z, key_padding_mask=pad)
        h = [p.squeeze(1), m.sum(1, keepdim=True)]
        if scal is not None: h.append(scal)
        return s.head(torch.cat(h, -1)).squeeze(-1)

model = SetTransformer(ndesc=NDESC, nscal=NSCAL).to(dev)
model.load_state_dict(ck['state']); model.eval()

sel = np.flatnonzero((label != 'unlabelled') & (A0 >= 0) & DSOK)
print(f'scoring {len(sel):,} precursors with a located candidate and a DScore')
scores = np.zeros(len(sel), dtype=np.float32)
with torch.no_grad():
    for a in range(0, len(sel), 512):
        idx = sel[a:a + 512]
        t = torch.from_numpy(X0[idx])
        cen = A0[idx].astype(int); lo = np.clip(cen - W // 2, 0, C - W)
        t = torch.stack([t[k, :, x:x + W] for k, x in enumerate(lo)])
        m = torch.from_numpy(M[idx].astype(np.float32))
        bg = t.median(dim=2, keepdim=True).values
        gt = torch.from_numpy(GAUSS).view(1, 1, -1).expand(t.shape[0], t.shape[1], -1)
        ch = torch.stack([torch.log1p(t), torch.clamp(t - bg, min=0).log1p(),
                          (t > 0).float(), gt], 2)
        tb = torch.clamp(t - bg, min=0)
        tc_ = tb - tb.mean(2, keepdim=True); gc_ = gt - gt.mean(2, keepdim=True)
        den = torch.sqrt((tc_ ** 2).sum(2) * (gc_ ** 2).sum(2)).clamp(min=1e-9)
        d = [torch.log1p(t.sum(2)), (tc_ * gc_).sum(2) / den,
             torch.from_numpy(RI[idx]), torch.from_numpy(PM[idx]),
             torch.from_numpy(FC[idx]), torch.from_numpy(OR[idx]),
             torch.from_numpy(SE[idx])]
        pd_ = torch.from_numpy(PD[idx])
        d += [pd_[:, k:k + 1].expand(-1, t.shape[1]) for k in range(pd_.shape[1])]
        sc = None
        if NSCAL:
            sc = torch.from_numpy(np.concatenate(
                [(SC0[idx] - SC_MU) / SC_SD, OK0[idx].astype(np.float32)[:, None]],
                axis=1).astype(np.float32)).to(dev)
        scores[a:a + len(idx)] = model(ch.to(dev), m.to(dev),
                                       torch.stack(d, -1).to(dev), sc).float().cpu().numpy()

is_ent = (label == 'ent')[sel]; is_tgt = (dec == 0)[sel]
n_ent = int((label[sel] == 'ent').sum()); n_ord = int((is_tgt & ~is_ent).sum())
r = n_ent / max(n_ord, 1)
print(f'entrapment ratio on the scored set: {n_ent:,}/{n_ord:,} = {r:.5f}')

def curve(score, name):
    o = np.argsort(-score)
    e = np.cumsum(is_ent[o]); t = np.cumsum(is_tgt[o] & ~is_ent[o])
    fdp = (e / r) / np.maximum(t, 1)
    print(f'\n  {name}')
    for g in (0.01, 0.02, 0.05, 0.10):
        k = np.flatnonzero((fdp <= g) & (np.arange(len(fdp)) > 300))
        if len(k):
            i = k[-1]
            print(f'    FDP <= {g:5.0%}   targets {int(t[i]):>7,}   entrapment {int(e[i]):>5,}'
                  f'   depth {i+1:>7,}')
        else:
            print(f'    FDP <= {g:5.0%}   (never reached)')
    return t, e

curve(DS[sel], 'ODIA shipped DScore')
curve(scores, 'trace model')
