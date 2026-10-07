# -*- coding: utf-8 -*-
"""smoke test: one packing of a family, timing, anomalies, termination measures of both flows."""
import sys, time, random
from fractions import Fraction as Fr
import tbx_core as C
import tbx_checks as K
from tbx_run import make, resolve_theta

fam = sys.argv[1] if len(sys.argv) > 1 else 'random'
k = int(sys.argv[2]) if len(sys.argv) > 2 else 12
seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
tmode = sys.argv[4] if len(sys.argv) > 4 else '0.08'
rng = random.Random(seed)
t0 = time.time()
P = make(fam, k, rng)
print('gen', fam, 'k', k, 'N', P.N, 'W', P.W, 'meta', P.meta, 'verify', P.verify(), '%.2fs' % (time.time() - t0))
h = P.k / 2 - 1
sinT = resolve_theta(tmode, P)
for side, PP in (('floor', P), ('ceil', P.reflect())):
    LD = K.LineData(PP)
    t1 = time.time()
    F = C.Flow(PP, sinT, Fr(1, 20), h).run()
    tot = K.losses(F)
    st, viol = K.path_checks(F, LD, side)
    print(side, 'recs', len(F.recs), 'R1multi', F.nR1multi, 'anom', F.anom[:3], '%.2fs' % (time.time() - t1),
          {t: round(float(v), 4) for t, v in tot.items()}, 'npass', st['npass'], 'maxpass', st['maxpass'],
          'Dmax %.5f' % float(st['L3_Dmax']), 'viol', viol[:3])
print('rss', C.rss_mb())
