# -*- coding: utf-8 -*-
"""tbx_point.py -- independent exact point tracer (free paths) for cross-checking tbx_core.Flow.

Different algorithm from the piecewise tracer: for ONE rational floor point x it follows the free path
(every entry candidate treated as won, [Q, Definition 3.10]) by Cyrus-Beck clipping of the ray against
every square (no spatial grid, no edge envelope, no pieces).  All arithmetic is exact (Fraction).
"""
from fractions import Fraction as Fr

F0 = Fr(0); F1 = Fr(1); HALF = Fr(1, 2)


def first_hit(sq, p, d):
    """Cyrus-Beck: least t >= 0 with p + t d in the closed square, or None."""
    t_in = F0; t_out = None
    c = (sq.cx, sq.cy)
    for (A, B, m) in sq.E:
        # constraint (p + t d - A).m <= 0
        num = (p[0] - A[0]) * m[0] + (p[1] - A[1]) * m[1]
        den = d[0] * m[0] + d[1] * m[1]
        if den == 0:
            if num > 0:
                return None
            continue
        t = -num / den
        if den < 0:      # entering
            if t > t_in:
                t_in = t
        else:            # leaving
            if t_out is None or t < t_out:
                t_out = t
    if t_out is not None and t_in > t_out:
        return None
    if t_out is not None and t_out < 0:
        return None
    return t_in


def classify(sq, q):
    """feature of the boundary point q of sq: 'bot' (relint), 'top', 'left', 'right', 'vertex', or 'interior'."""
    xi = (q[0] - sq.cx) * sq.u[0] + (q[1] - sq.cy) * sq.u[1]
    ze = (q[0] - sq.cx) * sq.n[0] + (q[1] - sq.cy) * sq.n[1]
    axi = abs(xi); aze = abs(ze)
    if axi > HALF or aze > HALF:
        return 'outside'
    if axi == HALF and aze == HALF:
        return 'vertex'
    if aze == HALF:
        return 'bot' if ze < 0 else 'top'
    if axi == HALF:
        return 'right' if xi > 0 else 'left'
    return 'interior'


def free_path(P, x, sinA, delta, hF, stop_at=None):
    """free path of the floor point x.  Returns dict(events, term, passes) where
    passes = [(j, q, gcum)] (entry point and cumulative gap at the entry),
    term = (type, point, j or None, gcum)."""
    S = P.S; k = P.k
    p = (Fr(x), F0); d = (F0, F1); src = -1; g = F0
    passes = []; cands = []
    for _ in range(10000):
        # contact
        tS = None; jS = None
        for Y in S:
            if Y.i == src:
                continue
            t = first_hit(Y, p, d)
            if t is None:
                continue
            if tS is None or t < tS:
                tS = t; jS = Y.i
        tD = delta - g
        tW = None
        if d[0] > 0:
            tW = (k - p[0]) / d[0]
        elif d[0] < 0:
            tW = -p[0] / d[0]
        tH = (hF - p[1]) / d[1]
        cand = [(tD, 0)]
        if tW is not None:
            cand.append((tW, 1))
        if tS is not None:
            cand.append((tS, 2))
        cand.append((tH, 3))
        tmin = min(c[0] for c in cand)
        prio = min(c[1] for c in cand if c[0] == tmin)
        q = (p[0] + tmin * d[0], p[1] + tmin * d[1])
        if prio == 0:
            return dict(passes=passes, cands=cands, term=('D', q, None, g + tmin))
        if prio == 1:
            return dict(passes=passes, cands=cands, term=('W', q, None, g + tmin))
        if prio == 3:
            return dict(passes=passes, cands=cands, term=('H', q, None, g + tmin))
        Y = S[jS]
        feat = classify(Y, q)
        g2 = g + tmin
        if feat != 'bot':
            return dict(passes=passes, cands=cands, term=('E', q, jS, g2, feat))
        cands.append((jS, q, g2))
        if stop_at is not None and jS == stop_at:
            return dict(passes=passes, cands=cands, term=('STOP', q, jS, g2))
        if Y.sa >= sinA:
            return dict(passes=passes, cands=cands, term=('T', q, jS, g2))
        passes.append((jS, q, g2))
        e = (q[0] + Y.n[0], q[1] + Y.n[1])
        if e[1] >= hF:
            return dict(passes=passes, cands=cands, term=('H', e, None, g2))
        p = e; d = Y.n; src = jS; g = g2
    raise RuntimeError('too many steps')
