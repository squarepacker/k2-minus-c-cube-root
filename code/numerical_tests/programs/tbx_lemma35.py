# -*- coding: utf-8 -*-
"""[R, Lemma 3.5] (used in Lemma 5 of the draft) on every elementary line interval of generated packings, exactly:
    sum_{0 < c_S(y) < 1} c_S(y) >= 1 - omega(y) - E(y).
On an elementary interval the three functions are affine; the check is done at both ends (closure) and at the
midpoint.  Negative control: the same with 2 in place of 1 (fails exactly on lines with k - 1 long chords).
usage: python tbx_lemma35.py OUT.json NPACK_PER_FAMILY SEED
"""
import sys, json, random, time
from fractions import Fraction as Fr
import tbx_checks as K
from tbx_run import make

out = sys.argv[1]; npk = int(sys.argv[2]); seed = int(sys.argv[3])
fams = ['random', 'rows', 'column', 'deepcol', 'phase', 'lshape', 'rotgrid', 'rotband', 'stair', 'converge', 'wall']
rng = random.Random(seed)
res = dict(packings=0, intervals=0, viol=0, min_slack=None, max_nlong=0, neg2_fail_intervals=0, neg2_fail_packings=0,
           kmax=0, worst=None)
t0 = time.time()
for fam in fams:
    for k in (14, 16):
        for i in range(npk):
            if time.time() - t0 > 480:
                break
            P = make(fam, k, rng)
            LD = K.LineData(P)
            res['packings'] += 1
            negp = False
            for e in LD.el:
                res['intervals'] += 1
                om0, om1 = e['om']; E0, E1 = e['E']; s0, s1 = e['sh']
                for y in (e['lo'], e['hi'], (e['lo'] + e['hi']) / 2):
                    sl = (s0 + s1 * y) - (1 - (om0 + om1 * y) - (E0 + E1 * y))
                    if res['min_slack'] is None or sl < res['min_slack']:
                        res['min_slack'] = sl
                        res['worst'] = (fam, k, i, float(y), e['nlong'])
                    if sl < 0:
                        res['viol'] += 1
                if (s0 + s1 * (e['lo'] + e['hi']) / 2) < 2 - (om0 + om1 * (e['lo'] + e['hi']) / 2) - (E0 + E1 * (e['lo'] + e['hi']) / 2):
                    res['neg2_fail_intervals'] += 1
                    negp = True
                res['max_nlong'] = max(res['max_nlong'], e['nlong'])
            res['neg2_fail_packings'] += 1 if negp else 0
res['min_slack'] = float(res['min_slack'])
res['seconds'] = time.time() - t0
json.dump(res, open(out, 'w'), indent=1)
print(res)
