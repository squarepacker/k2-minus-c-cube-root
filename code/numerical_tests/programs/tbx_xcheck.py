# -*- coding: utf-8 -*-
"""tbx_xcheck.py -- cross-checks of the exact piecewise tracer (tbx_core.Flow).

(1) Independent exact point tracer (tbx_point.free_path, Cyrus-Beck, no pieces) at random rational floor
    points of every record (and at record endpoints = degenerate points):
    - non-M records: the free path coincides with the flow path (same passes and entry points, same
      termination type, same termination point, exactly);
    - M records: the free path reaches the merge point exactly with the same gap; the R1 winner x' at that
      point is located in the winner records, its free path reaches the SAME point exactly, and
      (g', x') < (g, x) lexicographically;
    - at all points (including the degenerate record endpoints): Lemma 3 quantization at every entry
      candidate of the free path: q_y - m in [-D, delta].
(2) Comparison of the termination-type measures with the float tracer asmx_core.py of the earlier assembly
    test (independent code by another session), if it is importable.

usage: python tbx_xcheck.py OUT.json FAMILY K SEED NPACK DELTA THETA [NPTS]
"""
import sys, os, json, random, time
from fractions import Fraction as Fr
import math
import tbx_core as C
import tbx_gen as G
import tbx_point as PT
from tbx_run import make, jf

ASMX = os.environ.get('ASMX_DIR', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'asmx'))   # folder holding asmx_core.py of the k^{1/4} repository (code/numerical_tests/assembly/)


def rand_in(rng, a, b):
    t = Fr(rng.randrange(1, 10**6), 10**6 + 3)
    return a + (b - a) * t


def winners_at(flow, j):
    """records that R1-won at square j: (rec, Q0, Q1, gq0, gq1) for pass entries at j and T at j."""
    out = []
    for rec in flow.recs:
        if rec[2] == 'T' and rec[6][0] == j:
            out.append((rec, rec[6][1], rec[6][2], rec[6][3], rec[6][4]))
        for nd in C.node_chain(rec[3]):
            if nd[0] == 'P' and nd[2][0] == j:
                out.append((rec, nd[2][1], nd[2][2], nd[2][3], nd[2][4]))
    return out


def quant_check(P, fp, delta):
    """Lemma 3 on the free path: at each candidate q (after m passes): q_y - m in [-D, delta]."""
    bad = []
    D = Fr(0)
    m = 0
    passes = fp['passes']
    for (j, q, g) in fp['cands']:
        val = q[1] - m + D
        if val < 0 or val > delta:
            bad.append((j, float(val)))
        # this candidate is passed iff it appears in passes at position m
        if m < len(passes) and passes[m][0] == j:
            D += 1 - P.S[j].c
            m += 1
    return bad


def check_flow(P, flow, rng, npts, delta, sinA, hF):
    res = dict(points=0, agree=0, disagree=[], merges=0, merge_ok=0, merge_bad=[], quant_bad=[], degenerate=0)
    recs = flow.recs
    idxs = list(range(len(recs)))
    rng.shuffle(idxs)
    for ri in idxs[:npts]:
        rec = recs[ri]
        xa, xb, typ = rec[0], rec[1], rec[2]
        x = rand_in(rng, xa, xb)
        res['points'] += 1
        fp = PT.free_path(P, x, sinA, delta, hF)
        qb = quant_check(P, fp, delta)
        if qb:
            res['quant_bad'].append((float(x), qb[:2]))
        segs = C.path_segments(P.S, rec)
        fpasses = [sg[1] for sg in segs if sg[0] == 'P']
        if typ != 'M':
            tp = fp['term']
            same = [p[0] for p in fp['passes']] == fpasses
            T = (rec[4][0] + rec[5][0] * x, rec[4][1] + rec[5][1] * x)
            same = same and tp[0] == typ and tp[1] == T
            if typ in ('T', 'E'):
                jj = rec[6][0] if typ == 'T' else rec[6]
                same = same and tp[2] == jj
            if same:
                res['agree'] += 1
            else:
                res['disagree'].append(dict(x=str(x), typ=typ, tracer=tp[0], fpasses=fpasses,
                                            ppasses=[p[0] for p in fp['passes']]))
        else:
            res['merges'] += 1
            j, Q0, Q1, gq0, gq1 = rec[6]
            q = (Q0[0] + Q1[0] * x, Q0[1] + Q1[1] * x)
            g = gq0 + gq1 * x
            fp2 = PT.free_path(P, x, sinA, delta, hF, stop_at=j)
            ok = fp2['term'][0] == 'STOP' and fp2['term'][1] == q and fp2['term'][3] == g
            # locate the winner
            Y = P.S[j]
            xi = (q[0] - Y.cx) * Y.u[0] + (q[1] - Y.cy) * Y.u[1]
            found = None
            for (wr, W0, W1, wg0, wg1) in winners_at(flow, j):
                xi0 = (W0[0] - Y.cx) * Y.u[0] + (W0[1] - Y.cy) * Y.u[1]
                xi1 = W1[0] * Y.u[0] + W1[1] * Y.u[1]
                xw = (xi - xi0) / xi1
                if wr[0] < xw < wr[1]:
                    found = (xw, wg0 + wg1 * xw)
                    break
            if found is None:
                ok = False
                res['merge_bad'].append(dict(x=str(x), why='winner not found'))
            else:
                xw, gw = found
                fpw = PT.free_path(P, xw, sinA, delta, hF, stop_at=j)
                okw = fpw['term'][0] == 'STOP' and fpw['term'][1] == q and fpw['term'][3] == gw
                lex = (gw, xw) < (g, x)
                if not (okw and lex):
                    ok = False
                    res['merge_bad'].append(dict(x=str(x), xw=str(xw), okw=okw, lex=lex))
            if ok:
                res['merge_ok'] += 1
    # degenerate points: record endpoints
    ends = sorted(set([r[0] for r in recs] + [r[1] for r in recs]))
    for x in ends[1:-1][:npts]:
        fp = PT.free_path(P, x, sinA, delta, hF)
        res['degenerate'] += 1
        qb = quant_check(P, fp, delta)
        if qb:
            res['quant_bad'].append(('endpoint', float(x), qb[:2]))
    return res


def asmx_compare(P, delta, sinA, hF):
    try:
        sys.path.insert(0, ASMX)
        import asmx_core as AX
    except Exception as e:
        return dict(error='asmx_core not importable: %s' % e)
    sq = [(float(s.cx), float(s.cy), math.atan2(float(s.s), float(s.c))) for s in P.S]
    AP = AX.Packing(float(P.k), sq)
    alpha = float(C.asin_q(sinA))
    F = AX.Flow(AP, float(delta), alpha, float(hF))
    F.run(tlimit=120.0)
    tot = {}
    for r in F.recs:
        tot[r[2]] = tot.get(r[2], 0.0) + (r[1] - r[0])
    return dict(measures=tot, anom=F.anom[:5], dropped=F.dropped)


def main():
    outp = sys.argv[1]; fam = sys.argv[2]; k = int(sys.argv[3]); seed = int(sys.argv[4])
    npack = int(sys.argv[5]); delta = Fr(sys.argv[6]); theta = sys.argv[7]
    npts = int(sys.argv[8]) if len(sys.argv) > 8 else 60
    rng = random.Random(seed)
    out = []
    t0 = time.time()
    for ip in range(npack):
        if time.time() - t0 > 400:
            break
        P = make(fam, k, rng)
        hF = P.k / 2 - 1
        if theta == 'a+':
            c, s = G.rot(P.meta['a'])
            sinA = abs(s) + Fr(1, 10**6)
        else:
            sinA = C.sin_rat(float(theta))
        r = dict(fam=fam, k=k, seed=seed, idx=ip, N=P.N, sinA=float(sinA), delta=float(delta))
        for side, PP in (('floor', P), ('ceil', P.reflect())):
            F = C.Flow(PP, sinA, delta, hF).run()
            mine = {}
            for rec in F.recs:
                mine[rec[2]] = mine.get(rec[2], Fr(0)) + (rec[1] - rec[0])
            rr = check_flow(PP, F, rng, npts, delta, sinA, hF)
            ax = asmx_compare(PP, delta, sinA, hF)
            diff = None
            if 'measures' in ax:
                diff = max(abs(float(mine.get(t, 0)) - ax['measures'].get(t, 0.0)) for t in C.TYPES)
            r[side] = dict(nrec=len(F.recs), anom=F.anom[:5], R1multi=F.nR1multi,
                           measures={t: float(v) for t, v in mine.items()}, point=rr, asmx=ax, asmx_maxdiff=diff)
        out.append(jf(r))
    with open(outp, 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=1)
    print('done', len(out), 'packings', '%.1f s' % (time.time() - t0), 'rss', C.rss_mb())


if __name__ == '__main__':
    main()
