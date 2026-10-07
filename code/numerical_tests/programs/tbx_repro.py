# -*- coding: utf-8 -*-
"""Exact reproducer: regenerate packing number IDX of a campaign job (same seed => same random stream), write its
squares as exact rationals, and list the entry candidates on type-Z squares (Lemma 4) with exact data.

usage: python tbx_repro.py OUT.json FAMILY K SEED IDX DELTA THETAMODE TAU0 TAU1 BETA
  e.g. python tbx_repro.py out/repro_stair.json stair 16 2027 0 1/20 a+ 3/10 15/100 0.08
"""
import sys, json, random
from fractions import Fraction as Fr
import tbx_core as C
import tbx_checks as K
from tbx_run import make, resolve_theta, chain_theta, chain_theta_fp, TAU1_CHAIN, THETA_BAR

outp, fam, k, seed, idx = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
delta = Fr(sys.argv[6]); tmode = sys.argv[7]
tau0 = Fr(sys.argv[8]); tau1 = Fr(sys.argv[9]); beta = float(sys.argv[10])
rng = random.Random(seed)
for i in range(idx + 1):
    P = make(fam, k, rng)
ok, msg = P.verify()
LD = K.LineData(P)
if tmode == 'chain':
    sinT, _ = chain_theta(LD, TAU1_CHAIN, THETA_BAR)
elif tmode == 'chainfp':
    sinT, _ = chain_theta_fp(LD, TAU1_CHAIN)
else:
    sinT = resolve_theta(tmode, P)
h = P.k / 2 - 1
sinB = C.sin_rat(beta)
Z = K.typeZ(P, sinB, tau0)
out = dict(family=fam, k=k, seed=seed, idx=idx, verified=ok, verify_msg=msg, N=P.N, W=str(P.W), meta=str(P.meta),
           sinT=str(sinT), delta=str(delta), tau0=str(tau0), tau1=str(tau1), sinB=str(sinB),
           w_beta=str(K.wfun(sinB)), squares=[[str(v) for v in t] for t in P.raw], typeZ=Z, entries=[])
for side, PP in (('floor', P), ('ceil', P.reflect())):
    F = C.Flow(PP, sinT, delta, h).run()
    for rec in F.recs:
        segs = C.path_segments(PP.S, rec)
        passes = [sg for sg in segs if sg[0] == 'P']
        ents = [(m, sg[1], sg[2], sg[3]) for m, sg in enumerate(passes)]
        if rec[2] in ('M', 'T'):
            ents.append((len(passes), rec[6][0], rec[6][1], rec[6][2]))
        for (m, j, Q0, Q1) in ents:
            if j not in Z:
                continue
            D = sum((1 - PP.S[passes[i][1]].c for i in range(m)), Fr(0))
            xm = (rec[0] + rec[1]) / 2
            qy = Q0[1] + Q1[1] * xm
            out['entries'].append(dict(side=side, square=j, x_interval=[str(rec[0]), str(rec[1])], x_mid=str(xm),
                                       passes_before=m, passed=[sg[1] for sg in passes[:m]], D=str(D), D_float=float(D),
                                       D_le_tau1=bool(D <= tau1), q_y_mid=str(qy), dist_qy_Z=float(K.fdist(qy)),
                                       termination=rec[2]))
json.dump(out, open(outp, 'w'), indent=1)
print('verified', ok, msg, '| type-Z squares', len(Z), '| entries on type-Z squares', len(out['entries']))
for e in out['entries'][:10]:
    print(' ', e['side'], 'square', e['square'], 'x in', [round(float(Fr(v)), 6) for v in e['x_interval']], 'm', e['passes_before'],
          'D = %.6f' % e['D_float'], 'D <= tau1:', e['D_le_tau1'], 'dist(q_y, Z) = %.4f' % e['dist_qy_Z'], e['termination'])
