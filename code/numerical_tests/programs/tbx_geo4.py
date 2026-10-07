# -*- coding: utf-8 -*-
"""Geometric core of Lemma 4 (proof of the draft), by case, on all slightly tilted squares of generated packings.

For a square S (0 < a < beta) and a height y in Phi(S) with dist(y, Z) >= tau0, the proof claims that every point
z of Bot(S) has dist(z_y, Z) >= tau0 - w, with
   bottom-ramp case  y in (v, v + sin a cos a):            |z_y - y|     <  sin a
   top-ramp case     y in (top - sin a cos a, top):        |z_y - (y-1)| <  (1 - cos a) + sin a
so w = (1 - cos a) + sin a works in both cases (the draft writes 1.0001 beta, valid for beta <= 2e-4).
Checked exactly at the level of intervals: for each case and each tau0, the worst y is the point of
Phi-part cap H closest to the bottom-side height range; we compute
   margin(case) = min_{z in Bot} dist(z_y, Z) - (tau0 - w_case(a))   over the part of Phi(S) cap H of that case.
Control: the same with w replaced by 1.0001 a (the draft's constant): counts squares where it would fail.
usage: python tbx_geo4.py OUT.json NPACK_PER_FAMILY SEED
"""
import sys, json, random, time, math
from fractions import Fraction as Fr
import mpmath as mp
import tbx_core as C
import tbx_checks as K
from tbx_run import make

out = sys.argv[1]; npk = int(sys.argv[2]); seed = int(sys.argv[3])
rng = random.Random(seed)
fams = ['random', 'rows', 'column', 'deepcol', 'phase', 'lshape', 'rotgrid', 'rotband', 'stair', 'converge', 'wall']
taus = [Fr(1, 10), Fr(2, 10), Fr(3, 10), Fr(4, 10), Fr(49, 100)]
res = dict(packings=0, squares=0, cases={'bottom': 0, 'top': 0}, min_margin={'bottom': None, 'top': None},
           viol=0, lin_fail={'bottom': 0, 'top': 0}, max_a=0.0)
t0 = time.time()


def H_points(a, b, tau0):
    """closed set (a, b) cap H as list of closed intervals [c, d] (a, b open ends kept as limits)."""
    out_ = []
    for n in range(math.floor(a) - 1, math.floor(b) + 2):
        c = max(a, n + tau0); d = min(b, n + 1 - tau0)
        if d >= c and not (c == b or d == a):
            out_.append((c, d))
    return out_


for fam in fams:
    for k in (14, 16):
        for i in range(npk):
            if time.time() - t0 > 480:
                break
            P = make(fam, k, rng)
            res['packings'] += 1
            for S in P.S:
                if S.sa == 0:
                    continue
                res['squares'] += 1
                a = S.a
                res['max_a'] = max(res['max_a'], float(a))
                lo_b, hi_b = S.v, S.v + S.sa                 # heights of Bot(S)
                mdb = K.mindist_range(lo_b, hi_b)
                (p0, p1), (q0, q1) = S.phi_set()
                for tau0 in taus:
                    for case, (u0, u1) in (('bottom', (p0, p1)), ('top', (q0, q1))):
                        if not H_points(u0, u1, tau0):
                            continue
                        res['cases'][case] += 1
                        w = (1 - mp.cos(a)) + mp.sin(a)
                        marg = C.mpq(mdb) - (C.mpq(tau0) - w)
                        mm = res['min_margin'][case]
                        if mm is None or marg < mm:
                            res['min_margin'][case] = marg
                        if marg < 0:
                            res['viol'] += 1
                        if C.mpq(mdb) < C.mpq(tau0) - mp.mpf('1.0001') * a:
                            res['lin_fail'][case] += 1
res['min_margin'] = {k: (float(v) if v is not None else None) for k, v in res['min_margin'].items()}
res['seconds'] = time.time() - t0
json.dump(res, open(out, 'w'), indent=1)
print(res)
