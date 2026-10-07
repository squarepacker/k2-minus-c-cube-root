# -*- coding: utf-8 -*-
"""tbx_run.py -- driver for the first-draft numerical test.

usage:  python tbx_run.py OUT.jsonl FAMILY K SEED NPACK DELTA THETAS [NY] [TLIMIT]
  FAMILY  random | column | phase | lshape | rotgrid | stair
  THETAS  comma list of 'chain' (theta from the measured chain, tau1 = 0.15, theta_bar = 0.25)
          and/or fixed angles in radians (e.g. chain,0.06,0.4)
Every packing is verified exactly before use.  Each process stops starting new work after TLIMIT seconds
(default 420) and aborts a flow that exceeds its own limits (time, 380 MB working set).
"""
import sys, json, time, random, math, traceback
from fractions import Fraction as Fr
import mpmath as mp
import tbx_core as C
import tbx_gen as G
import tbx_checks as K

TAU_SETS = [  # (tau0, tau1, beta, nu0)
    (Fr(3, 10), Fr(15, 100), 0.03, Fr(3, 10)),
    (Fr(3, 10), Fr(15, 100), 0.08, Fr(3, 10)),
    (Fr(3, 10), Fr(15, 100), 0.14, Fr(3, 10)),
    (Fr(2, 10), Fr(8, 100), 0.10, Fr(3, 10)),
    (Fr(49, 100), Fr(8, 100), 0.34, Fr(3, 10)),
    (Fr(49, 100), Fr(8, 100), 0.34, Fr(49, 100)),
    (Fr(3, 10), Fr(15, 100), 0.14, Fr(49, 100)),
    (Fr(2, 10), Fr(8, 100), 0.10, Fr(49, 100)),
]
THETA_BAR = 0.25
TAU1_CHAIN = Fr(15, 100)


def jf(x):
    """JSON-friendly conversion."""
    if x is None:
        return None
    if isinstance(x, Fr):
        return float(x)
    if isinstance(x, mp.mpf):
        return float(x)
    if isinstance(x, (bool, int, float, str)):
        return x
    if isinstance(x, dict):
        return {str(k): jf(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [jf(v) for v in x]
    return str(x)


def make(family, k, rng):
    if family == 'random':
        kw = dict(amax=rng.choice([0.03, 0.045, 0.06, 0.1]), pdel=rng.choice([0.0, 0.05, 0.15]),
                  pbig=rng.choice([0.05, 0.15]), ptilt=rng.choice([0.3, 0.6]))
        kw['a'] = kw['amax']
    elif family == 'rows':
        kw = dict(amax=rng.choice([0.02, 0.04]), pdel=rng.choice([0.02, 0.08]))
        kw['a'] = kw['amax']
    elif family == 'column':
        kw = dict(a=rng.choice([0.02, 0.035, 0.045]), ncols=rng.choice([1, 2, 3]),
                  fill_angles=rng.choice([(0,), (0, 0.01, -0.01), (0, 0.004)]))
    elif family == 'phase':
        kw = dict(b=rng.choice([0.08, 0.12, 0.2]), abeta=rng.choice([0.01, 0.025]))
        kw['a'] = kw['b']
    elif family == 'lshape':
        kw = dict(w=rng.choice([1.15, 1.25, 1.4]), a=rng.choice([0.03, 0.05, 0.08]),
                  abeta=rng.choice([0.01, 0.02]))
    elif family == 'rotgrid':
        kw = dict(a=rng.choice([0.05, 0.12, 0.2, 0.3, 0.3, 0.32]))
    elif family == 'stair':
        kw = dict(a=rng.choice([0.3, 0.34]), m=rng.choice([5, 6, 7]), ab=rng.choice([0.01, 0.02]))
    elif family == 'converge':
        kw = dict(a=rng.choice([0.04, 0.08, 0.12]), pskip=rng.choice([0.0, 0.2]))
    elif family == 'wall':
        kw = dict(a=rng.choice([0.04, 0.08, 0.12]))
    elif family == 'rotband':
        kw = dict(a=rng.choice([0.28, 0.29, 0.3, 0.31, 0.32]), oy=Fr(rng.randrange(1, 5000), 10000))
    elif family == 'deepcol':
        kw = dict(a=rng.choice([0.08, 0.1, 0.12, 0.13, 0.14]), ncols=rng.choice([1, 2, 4]),
                  fill_angles=rng.choice([(0,), (0, 0.01, -0.01)]))
    else:
        raise ValueError(family)
    a_meta = kw.get('a')
    kwg = {kk: v for kk, v in kw.items() if not (kk == 'a' and family in ('random', 'rows', 'phase'))}
    gfam = {'deepcol': 'column', 'rotband': 'rotgrid'}.get(family, family)
    P = G.two_sided(gfam, rng, k, kwg)
    P.meta['a'] = a_meta
    return P


def resolve_theta(tmode, P):
    """'a+': just above the family angle a (|s| of rot(a) + 1e-6); a number: that angle (radians)."""
    if tmode == 'a+':
        c, s = G.rot(P.meta['a'])
        return abs(s) + Fr(1, 10**6)
    return C.sin_rat(float(tmode))


def chain_theta_fp(LD, tau1):
    """largest theta (bisection) with theta <= f(theta) := 2 tau1 cos(theta) / (A_meas(theta) W/2 + theta m_ge),
    m_ge = |{omega >= 1/2}|.  Then for every path: D <= (theta/2) sum a <= (theta/2)(A_meas W/2 + theta m_ge)/cos theta
    <= tau1  (Lemmas 1 and 3 with the measured constant; f is decreasing in theta)."""
    mge = C.mpq(LD.m_ge_half)
    lo, hi = 1e-6, 0.7
    for _ in range(50):
        mid = (lo + hi) / 2
        s = C.sin_rat(mid)
        th = C.asin_q(s)
        val = 2 * C.mpq(tau1) * mp.cos(th) / (LD.A_meas_W(s) / 2 + th * mge)
        if th <= val:
            lo = mid
        else:
            hi = mid
    s = Fr(int(mp.floor(mp.sin(mp.mpf(lo)) * 10**6)), 10**6)
    th = C.asin_q(s)
    fval = 2 * C.mpq(tau1) * mp.cos(th) / (LD.A_meas_W(s) / 2 + th * mge)
    return s, dict(theta_fp=th, f_at_theta=fval, ok=bool(th <= fval))


def chain_theta(LD, tau1, theta_bar):
    sb = C.sin_rat(theta_bar)
    thb = C.asin_q(sb)
    AW = LD.A_meas_W(sb)
    B1W = (AW / 2 + thb * C.mpq(LD.m_ge_half)) / mp.cos(thb)
    th = min(thb, 2 * C.mpq(tau1) / B1W)
    s = Fr(int(mp.floor(mp.sin(th) * 10**6)), 10**6)
    return s, dict(AW_bar=AW, B1W=B1W, theta_chain=th)


def evaluate(P, delta, thetas, ny, t_end):
    k = P.k; h = k / 2 - 1
    Pr = P.reflect()
    LD = K.LineData(P)
    LDr = K.LineData(Pr)
    out = []
    for tmode in thetas:
        if time.time() > t_end:
            break
        rec = dict(name=P.name, k=int(k), N=P.N, W=int(P.W), delta=float(delta), tmode=tmode)
        if tmode == 'chain':
            sinT, info = chain_theta(LD, TAU1_CHAIN, THETA_BAR)
            rec['chain'] = jf(info)
        elif tmode == 'chainfp':
            sinT, info = chain_theta_fp(LD, TAU1_CHAIN)
            rec['chain'] = jf(info)
        else:
            sinT = resolve_theta(tmode, P)
        rec['sinT'] = float(sinT)
        t0 = time.time()
        remaining = max(30.0, t_end - time.time())
        FF = C.Flow(P, sinT, delta, h, tlimit=remaining).run()
        FC = C.Flow(Pr, sinT, delta, h, tlimit=remaining).run()
        rec['t_flows'] = time.time() - t0
        rec['nrec'] = [len(FF.recs), len(FC.recs)]
        rec['anom'] = [FF.anom[:5], FC.anom[:5]]
        viol = []
        stF, v = K.path_checks(FF, LD, 'floor'); viol += v
        stC, v = K.path_checks(FC, LDr, 'ceil'); viol += v
        rec['paths'] = [jf(stF), jf(stC)]
        maxD = max(stF['L3_Dmax'], stC['L3_Dmax'])
        rec['maxD'] = float(maxD)
        ys = K.sample_heights(h, ny)
        sF = K.line_samples(FF, ys)
        sC = K.line_samples(FC, ys)
        # Lemma 6 sets: type-Z for each tau set; and all unpassed squares below h
        Zsets_F = {'all': list(range(P.N))}
        Zsets_C = {'all': list(range(P.N))}
        tz = {}
        for (tau0, tau1, beta, nu0) in TAU_SETS:
            sinB = C.sin_rat(beta)
            key = '%s_%s_%s' % (float(tau0), float(tau1), beta)
            if key in tz:
                continue
            Z = K.typeZ(P, sinB, tau0)       # H is reflection invariant: same indices for the ceiling
            tz[key] = Z
            Zsets_F['tZ_' + key] = Z
            Zsets_C['tZ_' + key] = Z
        l6F, v = K.lemma6_checks(FF, sF, Zsets_F, h, 'floor'); viol += v
        l6C, v = K.lemma6_checks(FC, sC, Zsets_C, h, 'ceil'); viol += v
        rec['L6'] = [jf(l6F), jf(l6C)]
        LtotF = l6F['LT'] + l6F['Lam_b']
        LtotC = l6C['LT'] + l6C['Lam_b']
        rec['Ltot'] = [float(LtotF), float(LtotC)]
        rec['Aprime_meas'] = float(C.mpq(l6F['LT'] + l6C['LT']) * C.asin_q(sinT) / C.mpq(P.W))
        l4 = {}; ta = {}
        for (tau0, tau1, beta, nu0) in TAU_SETS:
            sinB = C.sin_rat(beta)
            key = '%s_%s_%s' % (float(tau0), float(tau1), beta)
            Z = tz[key]
            wb = K.wfun(sinB)
            cond = bool(C.mpq(tau0) - wb > C.mpq(tau1) and tau1 > delta)
            r4F, v1 = K.lemma4_checks(FF, Z, tau0, tau1, 'floor')
            r4C, v2 = K.lemma4_checks(FC, Z, tau0, tau1, 'ceil')
            if cond:
                viol += v1 + v2
            l4[key] = dict(cond=cond, maxD_le_tau1=bool(maxD <= tau1), floor=jf(r4F), ceil=jf(r4C))
            if cond and maxD <= tau1 and (r4F['n_entry_pieces'] + r4C['n_entry_pieces']) > 0:
                viol.append(('L4-conclusion', key))
            keyA = key + '_nu%s' % float(nu0)
            stA, v = K.theoremA(P, LD, tau0, tau1, sinB, nu0, h, LtotF, LtotC, FF, FC, keyA)
            stA['cond'] = cond; stA['maxD_le_tau1'] = bool(maxD <= tau1)
            if cond and maxD <= tau1:
                viol += v
                if stA['ZKb_passed'] or stA['ZKt_passed'] or stA['ZK2b_passed'] or stA['ZK2t_passed']:
                    viol.append(('TA-ZK-passed', keyA))
            else:
                stA['viol_if_unconditional'] = [str(x) for x in v]
            ta[keyA] = jf(stA)
        rec['L4'] = l4
        rec['TA'] = ta
        rec['viol'] = [str(x) for x in viol]
        rec['t_total'] = time.time() - t0
        rec['rss'] = C.rss_mb()
        out.append(rec)
    return out


def main():
    outp = sys.argv[1]; fam = sys.argv[2]; k = int(sys.argv[3]); seed = int(sys.argv[4])
    npack = int(sys.argv[5]); delta = Fr(sys.argv[6]); thetas = sys.argv[7].split(',')
    ny = int(sys.argv[8]) if len(sys.argv) > 8 else 24
    tl = float(sys.argv[9]) if len(sys.argv) > 9 else 420.0
    t_start = time.time(); t_end = t_start + tl
    rng = random.Random(seed)
    with open(outp, 'a', encoding='utf-8') as f:
        for ip in range(npack):
            if time.time() > t_end - 20:
                break
            try:
                P = make(fam, k, rng)
            except Exception as e:
                f.write(json.dumps(dict(error='gen', fam=fam, msg=str(e))) + '\n'); f.flush()
                continue
            ok, msg = P.verify()
            hdr = dict(fam=fam, k=k, seed=seed, idx=ip, N=P.N, W=int(P.W), verified=ok, vmsg=msg)
            if not ok:
                f.write(json.dumps(hdr) + '\n'); f.flush()
                continue
            try:
                res = evaluate(P, delta, thetas, ny, t_end)
                for r in res:
                    r.update(hdr)
                    f.write(json.dumps(r) + '\n')
                f.flush()
            except C.Budget as e:
                hdr['budget'] = str(e)
                f.write(json.dumps(hdr) + '\n'); f.flush()
            except Exception as e:
                hdr['error'] = traceback.format_exc()[-1500:]
                f.write(json.dumps(hdr) + '\n'); f.flush()
    print('done', fam, k, seed, 'elapsed %.1f s' % (time.time() - t_start), 'rss', C.rss_mb())


if __name__ == '__main__':
    main()
