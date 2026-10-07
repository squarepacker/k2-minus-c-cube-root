# -*- coding: utf-8 -*-
"""tbx_checks.py -- the links of the first draft of the argument evaluated on one packing with measured constants.

Item 1  Lemma 1 (tilt budget) per path piece: geometric core exactly; per-path RHS and the measured
        global bound A_meas(theta) W/2 + theta |{omega >= 1/2}|  (and the draft form (A_meas/2 + 2 theta) W).
        Lemma 2 (lateral drift) exactly.
Item 2  Lemma 3 (quantization) at every entry candidate exactly; D <= (theta/2) sum a  (mpmath, 50 digits).
Item 3  Lemma 4: entries on 'type-Z' squares (0 < a < beta, Phi meets H); geometric core of the proof.
Item 4  Lemma 6 (shadow budget): global, integrated and pointwise (exact at sampled heights).
Item 5  Theorem A with measured constants (sets I, II, K exact; Lemma 5 chain exact).
Negative controls are evaluated on the same data.
"""
from fractions import Fraction as Fr
import math, bisect
import mpmath as mp
from tbx_core import (F0, F1, HALF, asin_q, mpq, Lines, lin_set, H_ivs, iv_inter, iv_measure,
                      path_segments, measure_lt, line_full, density_excess, TYPES)


def fdist(y):
    """dist(y, Z) for a Fraction y (exact)."""
    f = y - math.floor(y)
    return min(f, 1 - f)


def mindist_range(a, b):
    """min of dist(y, Z) over y in [a, b] (exact)."""
    if math.floor(a) != math.floor(b) or a == math.floor(a) or b == math.floor(b):
        return F0
    return min(fdist(a), fdist(b))


def maxdist_range(a, b):
    """max of dist(y, Z) over y in [a, b] (exact)."""
    best = max(fdist(a), fdist(b))
    n = math.floor(a)
    for m in (n, n + 1, n + 2):
        hm = Fr(m) + HALF
        if a <= hm <= b:
            best = HALF
    return best


def wfun(sin_b):
    b = asin_q(sin_b)
    return 1 - mp.cos(b) + mp.sin(b)


class LineData:
    """per-packing line data: omega, max inclination, E, short chords on elementary intervals."""

    def __init__(self, P):
        self.P = P
        self.L = Lines(P)
        self.el = self.L.el
        self.k = P.k
        # {omega < 1/2} and {omega >= 1/2} measures per elementary interval
        self.m_lt_half = []
        self.m_ge_half = F0
        for e in self.el:
            iv = lin_set(e['lo'], e['hi'], e['om'][0], e['om'][1], '<', HALF)
            m = (iv[1] - iv[0]) if iv else F0
            self.m_lt_half.append(m)
            self.m_ge_half += (e['hi'] - e['lo']) - m
        self._asin_cache = {}

    def asin(self, q):
        if q not in self._asin_cache:
            self._asin_cache[q] = asin_q(q)
        return self._asin_cache[q]

    def g_values(self, sinT):
        """min(maxincl(y), theta) on each elementary interval (mp)."""
        th = self.asin(sinT)
        out = []
        for e in self.el:
            ms = e['maxsin']
            out.append(th if ms >= sinT else self.asin(ms))
        return out

    def A_meas_W(self, sinT):
        """int_{omega < 1/2} 2 min(maxincl, theta) dy  (the measured A(theta) * W)."""
        g = self.g_values(sinT)
        tot = mp.mpf(0)
        for gi, m in zip(g, self.m_lt_half):
            if m > 0:
                tot += 2 * gi * mpq(m)
        return tot

    def prefix(self, sinT):
        g = self.g_values(sinT)
        pre = [mp.mpf(0)]
        for gi, e in zip(g, self.el):
            pre.append(pre[-1] + gi * mpq(e['hi'] - e['lo']))
        return g, pre

    def cum(self, gp, y):
        g, pre = gp
        if y <= 0:
            return mp.mpf(0)
        if y >= self.k:
            return pre[-1]
        i = bisect.bisect_right(self.L.lo_list, y) - 1
        e = self.el[i]
        return pre[i] + g[i] * mpq(y - e['lo'])


# ---------------------------------------------------------------------------------------------
#  items 1 and 2 (and Lemma 2) on every record of a flow
# ---------------------------------------------------------------------------------------------
def path_checks(flow, LD, label=''):
    S = flow.S; k = flow.k; delta = flow.delta
    sinT = flow.sinA
    th = asin_q(sinT)
    costh = mp.cos(th)
    gp = LD.prefix(sinT)
    AW = LD.A_meas_W(sinT)
    W = mpq(flow.P.W)
    BL1 = AW / 2 + th * mpq(LD.m_ge_half)                 # sharp measured bound of Lemma 1
    BL1_draft = (AW / W / 2 + 2 * th) * W                  # the draft's form with A := A_meas(theta)
    viol = []
    st = dict(nrec=len(flow.recs), npass=0, nentry=0, maxpass=0,
              L1_lhs_max=mp.mpf(0), L1_ratio_path_max=mp.mpf(0), L1_ratio_glob_max=mp.mpf(0),
              L1_ratio_draft_max=mp.mpf(0), L1_neg_halfA=0, L1_neg_ext_overlap=0, L1_neg_capped=0,
              L3_Dmax=F0, L3_ratio_D_bound_max=mp.mpf(0), L3_gapfrac_max=F0, L3_distq_max=F0,
              L3_neg_noD=0, L3_neg_nogap=0, L2_ratio_max=mp.mpf(0), L2W_ratio_max=mp.mpf(0),
              AW=AW, BL1=BL1, BL1_draft=BL1_draft, theta=th)
    for rec in flow.recs:
        xa, xb, typ = rec[0], rec[1], rec[2]
        xm = (xa + xb) / 2
        segs = path_segments(S, rec)
        gaps = [sg for sg in segs if sg[0] == 'G']
        passes = [sg for sg in segs if sg[0] == 'P']
        st['npass'] += len(passes)
        st['maxpass'] = max(st['maxpass'], len(passes))
        # ---------------- Lemma 1: geometric core ----------------
        lhs = mp.mpf(0)
        sum_sin = F0
        for idx, (_, j, Q0, Q1) in enumerate(passes):
            Z = S[j]
            if not (Z.sa < sinT):
                viol.append(('L1-pass-tilt>=theta', label, j))
            lhs += Z.a * mpq(Z.c)
            sum_sin += Z.sa
            for xx in (xa, xb):
                qy = Q0[1] + Q1[1] * xx
                xi = (Q0[0] + Q1[0] * xx - Z.cx) * Z.u[0] + (qy - Z.cy) * Z.u[1]
                ze = (Q0[0] + Q1[0] * xx - Z.cx) * Z.n[0] + (qy - Z.cy) * Z.n[1]
                if abs(xi) > HALF or ze != -HALF:
                    viol.append(('L1-entry-not-on-Bot', label, j))
                if qy < Z.v or qy + Z.c > Z.top:
                    viol.append(('L1-I_j-not-in-extent', label, j))
                if idx + 1 < len(passes):
                    Qn0, Qn1 = passes[idx + 1][2], passes[idx + 1][3]
                    qn = Qn0[1] + Qn1[1] * xx
                    if qn < qy + Z.c:
                        viol.append(('L1-I_j-overlap', label, j))
                    # negative control: intervals of length = vertical extent (cos a + sin a) overlap?
                    if qn < qy + Z.ext:
                        st['L1_neg_ext_overlap'] += 1
            # strictness of disjointness on the open piece: next entry strictly above at the midpoint
            if idx + 1 < len(passes):
                Qn0, Qn1 = passes[idx + 1][2], passes[idx + 1][3]
                if not (Qn0[1] + Qn1[1] * xm > Q0[1] + Q1[1] * xm + Z.c):
                    viol.append(('L1-I_j-touch', label, j))
        if passes:
            # per-path RHS at the midpoint: int over union I_j of min(maxincl, theta)
            rhs = mp.mpf(0)
            rhs_cap = mp.mpf(0)
            for (_, j, Q0, Q1) in passes:
                Z = S[j]
                a_ = Q0[1] + Q1[1] * xm
                rhs += LD.cum(gp, a_ + Z.c) - LD.cum(gp, a_)
            if lhs > rhs * (1 + mp.mpf(10) ** -40):
                viol.append(('L1-path', label, float(lhs), float(rhs)))
            rp = (lhs / rhs) if rhs > 0 else (mp.mpf(0) if lhs == 0 else mp.inf)
            st['L1_ratio_path_max'] = max(st['L1_ratio_path_max'], rp)
            # negative control: drop the term theta*|{omega >= 1/2}| (i.e. claim lhs <= A_meas W / 2)
            if lhs > AW / 2:
                st['L1_neg_noY2'] = st.get('L1_neg_noY2', 0) + 1
            st['L1_ratio_glob_max'] = max(st['L1_ratio_glob_max'], lhs / BL1)
            st['L1_ratio_draft_max'] = max(st['L1_ratio_draft_max'], lhs / BL1_draft)
            st['L1_lhs_max'] = max(st['L1_lhs_max'], lhs)
            if lhs > BL1:
                viol.append(('L1-global', label, float(lhs), float(BL1)))
            # negative controls: A halved; maxincl capped at theta/2
            if lhs > AW / 4 + th * mpq(LD.m_ge_half):
                st['L1_neg_halfA'] += 1
            if lhs > rhs / 2:
                st['L1_neg_capped'] += 1
        # ---------------- Lemma 2: lateral drift ----------------
        sum_a_c = lhs
        bound_exact = sum_sin + delta            # sum |sin phi_j| + delta  (exact)
        bound_draft = sum_a_c / costh + mpq(delta)
        drift = F0
        px0 = F0; px1 = F1   # start point (x, 0): p_x(x) = x
        for sg in segs:
            if sg[0] == 'G':
                _, A0, A1, d, t0, t1, src = sg
                e0 = A0[0] + t0 * d[0]; e1 = A1[0] + t1 * d[0]
            else:
                _, j, Q0, Q1 = sg
                e0 = Q0[0] + S[j].n[0]; e1 = Q1[0]
            for xx in (xa, xb):
                dd = abs(e0 + e1 * xx - xx)
                if dd > drift:
                    drift = dd
        if drift > bound_exact:
            viol.append(('L2-drift', label, float(drift), float(bound_exact)))
        if bound_draft > 0:
            st['L2_ratio_max'] = max(st['L2_ratio_max'], mpq(drift) / bound_draft)
        if typ == 'W':
            for xx in (xa, xb):
                mn = min(xx, k - xx)
                if mpq(mn) > bound_draft * (1 + mp.mpf(10) ** -40):
                    viol.append(('L2-wall', label, float(mn), float(bound_draft)))
                st['L2W_ratio_max'] = max(st['L2W_ratio_max'], mpq(mn) / bound_draft)
        # ---------------- Lemma 3 at every entry candidate ----------------
        entries = []
        for idx, (_, j, Q0, Q1) in enumerate(passes):
            entries.append((idx, j, Q0, Q1))
        if typ in ('M', 'T'):
            j, Q0, Q1, gq0, gq1 = rec[6]
            entries.append((len(passes), j, Q0, Q1))
        for (m, j, Q0, Q1) in entries:
            st['nentry'] += 1
            D = sum((1 - S[passes[i][1]].c for i in range(m)), F0)
            suma = sum((S[passes[i][1]].a for i in range(m)), mp.mpf(0))
            # gap part: sum_{i <= m} t_i(x) * d_y(i)
            G0 = F0; G1 = F0; gc0 = F0; gc1 = F0
            for i in range(m + 1):
                _, A0, A1, d, t0, t1, src = gaps[i]
                G0 += t0 * d[1]; G1 += t1 * d[1]
                gc0 += t0; gc1 += t1
            # identity  q_y - m + D == sum g_i cos phi_i  (as affine functions of x)
            if (Q0[1] - m + D != G0) or (Q1[1] != G1):
                viol.append(('L3-identity', label, j))
            for xx in (xa, xb):
                val = G0 + G1 * xx
                if val < 0 or val > delta:
                    viol.append(('L3-range', label, j, float(val)))
                st['L3_gapfrac_max'] = max(st['L3_gapfrac_max'], val / delta)
                qy = Q0[1] + Q1[1] * xx
                dq = fdist(qy)
                st['L3_distq_max'] = max(st['L3_distq_max'], dq)
                if dq > max(D, delta):
                    viol.append(('L3-dist', label, j))
                # negative controls: drop D / drop the gap term
                if not (0 <= qy - m < delta):
                    st['L3_neg_noD'] += 1
                if not (-D <= qy - m <= 0):
                    st['L3_neg_nogap'] += 1
            if not (G0 + G1 * xm < delta):
                viol.append(('L3-strict', label, j))
            if rec[2] in ('M', 'T') and m == len(passes):
                # cumulative gap at the entry equals sum of gap lengths
                _, _, _, gq0, gq1 = rec[6]
                if gq0 != gc0 or gq1 != gc1:
                    viol.append(('L3-gapsum', label, j))
            if D > st['L3_Dmax']:
                st['L3_Dmax'] = D
            if D > 0:
                bnd = th / 2 * suma
                if mpq(D) > bnd:
                    viol.append(('L3-D>theta/2*sum a', label, float(D), float(bnd)))
                st['L3_ratio_D_bound_max'] = max(st['L3_ratio_D_bound_max'], mpq(D) / bnd)
                # negative control: halve the bound (theta/4) sum a
                if mpq(D) > bnd / 2:
                    st['L3_neg_halfbound'] = st.get('L3_neg_halfbound', 0) + 1
    # predicted bound of the chain: D <= (theta/2) * BL1 / cos(theta)
    st['L3_D_chain_bound'] = th / 2 * BL1 / costh
    return st, viol


# ---------------------------------------------------------------------------------------------
#  item 3: type-Z squares and entries on them
# ---------------------------------------------------------------------------------------------
def typeZ(P, sinB, tau0):
    out = []
    for s in P.S:
        if not (0 < s.sa < sinB):
            continue
        hit = False
        for (a, b) in s.phi_set():
            n0 = math.floor(a) - 1
            for n in range(n0, math.floor(b) + 2):
                c, d = n + tau0, n + 1 - tau0
                if a < d and b > c:
                    hit = True
                    break
            if hit:
                break
        if hit:
            out.append(s.i)
    return out


def lemma4_checks(flow, Zset, tau0, tau1, label=''):
    """entries (pass entries, M, T) on squares of Zset; geometric core of the proof of Lemma 4."""
    S = flow.S
    viol = []
    st = dict(nZ=len(Zset), n_entry_pieces=0, entry_measure=F0, minD_on_Z=None, min_dist_on_Z=None,
              n_entries_D_le_tau1=0, geo_min_margin=None)
    Zs = set(Zset)
    for i in Zset:
        s = S[i]
        md = mindist_range(s.v, s.v + s.sa)
        marg = mpq(md) - (mpq(tau0) - wfun(s.sa))
        if st['geo_min_margin'] is None or marg < st['geo_min_margin']:
            st['geo_min_margin'] = marg
        if marg < 0:
            viol.append(('L4-geometry', label, i, float(marg)))
    for rec in flow.recs:
        xa, xb = rec[0], rec[1]
        segs = path_segments(S, rec)
        passes = [sg for sg in segs if sg[0] == 'P']
        ents = [(idx, sg[1], sg[2], sg[3]) for idx, sg in enumerate(passes)]
        if rec[2] in ('M', 'T'):
            j, Q0, Q1, gq0, gq1 = rec[6]
            ents.append((len(passes), j, Q0, Q1))
        for (m, j, Q0, Q1) in ents:
            if j not in Zs:
                continue
            D = sum((1 - S[passes[i][1]].c for i in range(m)), F0)
            qa = Q0[1] + Q1[1] * xa; qb = Q0[1] + Q1[1] * xb
            dq = mindist_range(min(qa, qb), max(qa, qb))
            st['n_entry_pieces'] += 1
            st['entry_measure'] += xb - xa
            if st['minD_on_Z'] is None or D < st['minD_on_Z']:
                st['minD_on_Z'] = D
            if st['min_dist_on_Z'] is None or dq < st['min_dist_on_Z']:
                st['min_dist_on_Z'] = dq
            if D <= tau1:
                st['n_entries_D_le_tau1'] += 1
                viol.append(('L4-entry-with-small-deficit', label, j, float(D)))
    return st, viol


# ---------------------------------------------------------------------------------------------
#  item 4: Lemma 6 and losses
# ---------------------------------------------------------------------------------------------
def losses(flow):
    tot = {t: F0 for t in TYPES}
    for rec in flow.recs:
        tot[rec[2]] += rec[1] - rec[0]
    return tot


def passed_set(flow):
    out = set()
    for rec in flow.recs:
        for sg in path_segments(flow.S, rec):
            if sg[0] == 'P':
                out.add(sg[1])
    return out


def int_L(flow, h):
    """exact int_0^h L(y) dy = sum over loss records of int (h - tau(x)) dx."""
    tot = F0
    for rec in flow.recs:
        if rec[2] == 'H':
            continue
        xa, xb = rec[0], rec[1]
        tm = rec[4][1] + rec[5][1] * (xa + xb) / 2
        if tm > h:
            raise RuntimeError('loss above h')
        tot += (xb - xa) * (h - tm)
    return tot


def sample_heights(h, n):
    ys = []
    for i in range(n):
        y = h * Fr(2 * i + 1, 2 * n) + Fr(1, 7919 + 13 * i)
        if 0 < y < h:
            ys.append(y)
    return ys


def line_samples(flow, ys):
    out = []
    for y in ys:
        lf = line_full(flow, y, check=True)
        out.append(lf)
    return out


def lemma6_checks(flow, samples, Zsets, h, label=''):
    """Zsets: dict name -> list of squares (S below h assumed / filtered here) with mu == 0 required."""
    S = flow.S
    viol = []
    tot = losses(flow)
    LT = tot['T']
    supOv = max((s['Ovgap'] for s in samples), default=F0)
    Lam_b = tot['D'] + tot['M'] + tot['E'] + tot['W'] + supOv
    passed = passed_set(flow)
    iL = int_L(flow, h)
    # quadrature of Ovgap over (0, h) with the sample heights (midpoint rule; approximate)
    Ovq = sum((s['Ovgap'] for s in samples), F0) * h / max(1, len(samples))
    st = dict(LT=LT, D=tot['D'], M=tot['M'], E=tot['E'], W=tot['W'], H=tot['H'], supOvgap=supOv,
              Lam_b=Lam_b, intL=iL, intOvgap_quad=Ovq, nsamples=len(samples), exc_hits=0,
              tracer_bad=0)
    for s in samples:
        if s['Bd'] != 0:
            st['exc_hits'] += 1
        if s.get('bad'):
            st['tracer_bad'] += len(s['bad'])
            viol.append(('tracer-line-consistency', label, str(s['y']), s['bad'][:3]))
        # shadow inequality Sh <= L_T(y) + Lambda(y) = L(y) + Ovgap(y)
        if s['Bd'] == 0 and s['Sh'] > s['Lsum'] + s['Ovgap']:
            viol.append(('P4f-shadow', label, str(s['y'])))
    for name, Z in Zsets.items():
        Zb = [i for i in Z if S[i].top <= h]
        Zok = [i for i in Zb if i not in passed]
        nent = len(Zb) - len(Zok)
        rhs = h * (LT + Lam_b)
        r = dict(nZ=len(Zb), nZ_passed=nent, nZ_used=len(Zok), rhs_global=rhs,
                 ratio_global=(Fr(len(Zok)) / rhs) if rhs > 0 else None,
                 ratio_int=(Fr(len(Zok)) / (iL + Ovq)) if (iL + Ovq) > 0 else None,
                 ptw_min_margin=None, ptw_max_ratio=F0, neg_with_passed_fail=0, neg_LTonly_fail=0)
        if len(Zok) > rhs:
            viol.append(('L6-global', label, name, len(Zok), float(rhs)))
        for s in samples:
            if s['Bd'] != 0:
                continue
            y = s['y']
            lhs = sum((S[i].chord(y) for i in Zok), F0)
            lhs_all = sum((S[i].chord(y) for i in Zb), F0)
            rr = s['Lsum'] + s['Ovgap']
            marg = rr - lhs
            if r['ptw_min_margin'] is None or marg < r['ptw_min_margin']:
                r['ptw_min_margin'] = marg
            if rr > 0:
                r['ptw_max_ratio'] = max(r['ptw_max_ratio'], lhs / rr)
            elif lhs > 0:
                r['ptw_max_ratio'] = Fr(10**9)
            if lhs > s['Sh'] or lhs > rr:
                viol.append(('L6-pointwise', label, name, str(y)))
            # negative control 1: include squares that ARE passed (hypothesis mu == 0 dropped)
            if lhs_all > rr:
                r['neg_with_passed_fail'] += 1
            # negative control 2: only L_T(y) (the losses Lambda dropped)
            if lhs > s['L']['T']:
                r['neg_LTonly_fail'] += 1
        st['Z_' + name] = r
    return st, viol


# ---------------------------------------------------------------------------------------------
#  item 5: Theorem A with measured constants
# ---------------------------------------------------------------------------------------------
def theoremA(P, LD, tau0, tau1, sinB, nu0, h, LtotF, LtotC, flowF, flowC, label=''):
    """LtotF / LtotC: L_T + Lambda_b of the floor / ceiling flow (measured).
    flowF / flowC: flows (to check mu == 0 on Z_K)."""
    k = P.k; S = P.S
    beta = asin_q(sinB)
    viol = []
    hb_hi = k / 2 - 3; ht_lo = k / 2 + 3
    Hb = H_ivs(F0, hb_hi, tau0) if hb_hi > 0 else []
    Ht = H_ivs(ht_lo, k, tau0) if ht_lo < k else []
    calH = iv_measure(Hb) + iv_measure(Ht)
    mI = F0; mII = F0; mK = F0
    Kb = []; Kt = []
    intK = F0         # int_K (1 - omega - E)
    Emax = None
    ZKb = set(); ZKt = set()
    # maximal-K variant: K'' = {y in calH : all a < beta, 1 - omega - E > 0}
    intK2 = F0; mK2 = F0; ZK2b = set(); ZK2t = set()
    for e in LD.el:
        lo, hi = e['lo'], e['hi']
        for (Hs, side) in ((Hb, 'b'), (Ht, 't')):
            for (a, b) in iv_inter([(lo, hi)], Hs):
                om0, om1 = e['om']
                ivI = lin_set(a, b, om0, om1, '>=', nu0)
                ivL = lin_set(a, b, om0, om1, '<', nu0)
                if ivI:
                    mI += ivI[1] - ivI[0]
                if ivL:
                    if e['maxsin'] >= sinB:
                        mII += ivL[1] - ivL[0]
                    else:
                        c, d = ivL
                        mK += d - c
                        (Kb if side == 'b' else Kt).append((c, d))
                        E0, E1 = e['E']
                        intK += (d - c) * (1 - (om0 + om1 * (c + d) / 2) - (E0 + E1 * (c + d) / 2))
                        for yy in (c, d):
                            Ev = E0 + E1 * yy
                            if Emax is None or Ev > Emax:
                                Emax = Ev
                        for i in e['short']:
                            (ZKb if side == 'b' else ZKt).add(i)
                if e['maxsin'] < sinB:
                    # K'': 1 - omega - E > 0 :  (1 - om0 - E0) - (om1 + E1) y > 0
                    E0, E1 = e['E']
                    c0 = 1 - om0 - E0; c1 = -(om1 + E1)
                    iv2 = lin_set(a, b, c0, c1, '>', F0)
                    if iv2:
                        c, d = iv2
                        mK2 += d - c
                        intK2 += (d - c) * (c0 + c1 * (c + d) / 2)
                        for i in e['short']:
                            (ZK2b if side == 'b' else ZK2t).add(i)
    st = dict(calH=calH, calH_lower=(1 - 2 * tau0) * (k - 8), I=mI, II=mII, K=mK, Kb=iv_measure(Kb),
              Kt=iv_measure(Kt), Emax=Emax, intK=intK, nZKb=len(ZKb), nZKt=len(ZKt),
              K2=mK2, intK2=intK2, nZK2b=len(ZK2b), nZK2t=len(ZK2t))
    if calH < st['calH_lower']:
        viol.append(('TA-measure-calH', label))
    if mI + mII + mK != calH:
        viol.append(('TA-partition', label))
    W = P.W
    # (I)
    st['ratio_I'] = mI / (W / nu0) if W > 0 else None
    # (II): measured line cost
    AWb = LD.A_meas_W(sinB)
    st['AW_beta'] = AWb
    rhs_II = AWb / (2 * beta)
    st['ratio_II'] = mpq(mII) / rhs_II if rhs_II > 0 else (mp.mpf(0) if mII == 0 else mp.inf)
    if mpq(mII) > rhs_II * (1 + mp.mpf(10) ** -40):
        viol.append(('TA-II', label))
    st['neg_II_halfA'] = bool(mpq(mII) > rhs_II / 2)
    # Lemma 5 chain on K (exact up to the final beta)
    def chain(ZKset, Kivs_b, Kivs_t, intval):
        sum_int_c = F0
        sum_sc = F0
        for i in ZKset:
            s = S[i]
            for (pa, pb) in s.phi_set():
                for (a, b) in iv_inter([(pa, pb)], sorted(Kivs_b + Kivs_t)):
                    # c_S is affine on (pa, pb) (ramp part)
                    sum_int_c += (b - a) * s.chord((a + b) / 2)
            sum_sc += s.sa * s.c
        return sum_int_c, sum_sc
    sic, ssc = chain(ZKb | ZKt, Kb, Kt, intK)
    st['L5_int_K'] = intK; st['L5_sum_int_c'] = sic; st['L5_sum_sinacosa'] = ssc
    st['L5_beta_nZ'] = beta * (len(ZKb) + len(ZKt))
    if intK > sic:
        viol.append(('L5-3.5-integrated', label, float(intK), float(sic)))
    if sic > ssc:
        viol.append(('L5-sinacosa', label))
    if mpq(ssc) > beta * (len(ZKb) + len(ZKt)) and (len(ZKb) + len(ZKt)) > 0:
        viol.append(('L5-beta', label))
    # Lemma 4 + 6 on Z_K: mu == 0 for Z_{K_b} in the floor flow and Z_{K_t} in the ceiling flow
    passF = passed_set(flowF); passC = passed_set(flowC)
    st['ZKb_passed'] = len([i for i in ZKb if i in passF])
    st['ZKt_passed'] = len([i for i in ZKt if i in passC])
    st['ZKb_above_h'] = len([i for i in ZKb if S[i].top > h])
    st['ZKt_below_k-h'] = len([i for i in ZKt if S[i].v < k - h])
    st['L6b_ok'] = bool(len(ZKb) <= h * LtotF)
    st['L6t_ok'] = bool(len(ZKt) <= h * LtotC)
    # final |K| bound
    Ltot = LtotF + LtotC
    st['Ltot'] = Ltot
    den = (1 - nu0 - Emax) if Emax is not None else (1 - nu0)
    st['den'] = den
    if den > 0:
        boundK = beta * mpq(h * Ltot) / mpq(den)
        st['boundK'] = boundK
        st['ratio_K'] = mpq(mK) / boundK if boundK > 0 else (mp.mpf(0) if mK == 0 else mp.inf)
        if mpq(mK) > boundK:
            viol.append(('TA-K', label, float(mK), float(boundK)))
        # direct Lemma 5 inequality (1 - nu0 - Emax)|K| <= beta |Z_K|
        if mK > 0 and mpq(den * mK) > beta * (len(ZKb) + len(ZKt)):
            viol.append(('L5-final', label))
    else:
        st['boundK'] = None; st['ratio_K'] = None
    # K'' chain: int_{K''} (1 - omega - E) <= beta |Z_K''| <= beta h L_tot (per side)
    st['K2_chain_ratio'] = (mpq(intK2) / (beta * (len(ZK2b) + len(ZK2t)))) if (len(ZK2b) + len(ZK2t)) else None
    st['ZK2b_passed'] = len([i for i in ZK2b if i in passF])
    st['ZK2t_passed'] = len([i for i in ZK2t if i in passC])
    st['K2_L6_ratio_b'] = (Fr(len(ZK2b)) / (h * LtotF)) if LtotF > 0 else None
    st['K2_L6_ratio_t'] = (Fr(len(ZK2t)) / (h * LtotC)) if LtotC > 0 else None
    # total inequality of Theorem A with measured constants
    tot_rhs = mpq(W / nu0) + AWb / (2 * beta) + (st['boundK'] if st['boundK'] is not None else mp.inf)
    st['TA_rhs'] = tot_rhs
    st['TA_ratio'] = mpq(calH) / tot_rhs
    if mpq(calH) > tot_rhs:
        viol.append(('TA-total', label))
    return st, viol
