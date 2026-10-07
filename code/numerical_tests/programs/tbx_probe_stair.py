# -*- coding: utf-8 -*-
"""probe: does the stair target exist, is it type-Z, and which paths reach it?"""
import sys, random
from fractions import Fraction as Fr
import tbx_core as C
import tbx_gen as G
import tbx_checks as K

seed = int(sys.argv[1]) if len(sys.argv) > 1 else 1
a = float(sys.argv[2]) if len(sys.argv) > 2 else 0.3
m = int(sys.argv[3]) if len(sys.argv) > 3 else 7
delta = Fr(sys.argv[4]) if len(sys.argv) > 4 else Fr(1, 10)
rng = random.Random(seed)
P = G.two_sided('stair', rng, 16, dict(a=a, m=m, ab=0.02))
print('N', P.N, P.verify())
c, s = G.rot(a)
steep = [S for S in P.S if S.sa == abs(s) and S.cy < 8]
small = [S for S in P.S if 0 < S.sa < Fr(1, 20) and S.cy < 8]
print('steep squares (bottom half):', len(steep), ' slightly tilted (bottom half):', [(S.i, float(S.v), float(S.sa)) for S in small])
for S in small:
    print('  target', S.i, 'v', float(S.v), 'Phi', [(float(x), float(y)) for x, y in S.phi_set()])
sinT = abs(s) + Fr(1, 10**6)
F = C.Flow(P, sinT, delta, P.k / 2 - 1).run()
print('losses', {t: round(float(v), 4) for t, v in K.losses(F).items()})
# entries on the small squares
for rec in F.recs:
    segs = C.path_segments(P.S, rec)
    js = [sg[1] for sg in segs if sg[0] == 'P']
    ent = js + ([rec[6][0]] if rec[2] in ('M', 'T') else [])
    hit = [j for j in ent if j in [S.i for S in small]]
    if len([j for j in js if j in [S.i for S in steep]]) >= m - 1 or hit:
        print(' rec', rec[2], '(%.4f, %.4f)' % (float(rec[0]), float(rec[1])), 'passes', len(js), 'steep', sum(1 for j in js if j in [S.i for S in steep]),
              'term height %.4f' % float(rec[4][1] + rec[5][1] * (rec[0] + rec[1]) / 2), 'small hit', hit)
for tau0 in (Fr(3, 10), Fr(2, 10)):
    print('typeZ tau0', float(tau0), K.typeZ(P, C.sin_rat(0.08), tau0)[:10])
