# -*- coding: utf-8 -*-
"""tbx_core.py -- exact rational core for the numerical test of the first draft of the argument (numbering: see README.md).

Written from scratch on 2026-10-07 for this test (it does not import any earlier tracer).

* Squares have rational centres and rational unit side vectors (c, s), c^2 + s^2 = 1
  (rotations by angles phi with tan(phi/2) rational).  Every vertex is rational, so every
  quantity of the flows (contact parameters, gaps, heights, pieces, R1 comparisons) is rational
  and is computed with fractions.Fraction, without tolerances.
* Angles (inclinations a(S) = asin|s|) are transcendental; they are only used in the
  inequalities of the lemmas and are evaluated with mpmath at 50 digits.
  Threshold comparisons a(S) >= theta are done exactly as |s| >= sin(theta) with sin(theta)
  a rational number (thresholds are specified by their sines).

Conventions (paper [Q], Definitions 3.9--3.10, 3.20): phase phi in (-pi/4, pi/4],
u = (cos phi, sin phi), n = (-sin phi, cos phi), Bot(S) = c - n/2 + xi u (|xi| <= 1/2),
priority D > W > contact > H, entry candidate = contact in relint Bot(Y),
R1 = lexicographic min of (g, x) among the paths that reach (Y, q) alive, then R2
(T if a(Y) >= alpha_F, else pass [q, q + n_Y]), squares processed by increasing centre height.
"""
from fractions import Fraction as Fr
import math, heapq, time, bisect, ctypes, os
import mpmath as mp

mp.mp.dps = 50
F0 = Fr(0)
F1 = Fr(1)
HALF = Fr(1, 2)


# ----------------------------------------------------------------------------------------
#  memory helper (Windows: working set of this process, in MB)
# ----------------------------------------------------------------------------------------
class _PMC(ctypes.Structure):
    _fields_ = [('cb', ctypes.c_ulong), ('PageFaultCount', ctypes.c_ulong),
                ('PeakWorkingSetSize', ctypes.c_size_t), ('WorkingSetSize', ctypes.c_size_t),
                ('QuotaPeakPagedPoolUsage', ctypes.c_size_t), ('QuotaPagedPoolUsage', ctypes.c_size_t),
                ('QuotaPeakNonPagedPoolUsage', ctypes.c_size_t), ('QuotaNonPagedPoolUsage', ctypes.c_size_t),
                ('PagefileUsage', ctypes.c_size_t), ('PeakPagefileUsage', ctypes.c_size_t)]


def rss_mb():
    """(current, peak) working set of this process in MB (Windows); (-1, -1) if unavailable."""
    try:
        from ctypes import wintypes
        k32 = ctypes.WinDLL('kernel32', use_last_error=True)
        k32.GetCurrentProcess.restype = wintypes.HANDLE
        f = k32.K32GetProcessMemoryInfo
        f.argtypes = [wintypes.HANDLE, ctypes.POINTER(_PMC), wintypes.DWORD]
        f.restype = wintypes.BOOL
        pmc = _PMC()
        pmc.cb = ctypes.sizeof(_PMC)
        if not f(k32.GetCurrentProcess(), ctypes.byref(pmc), pmc.cb):
            return -1.0, -1.0
        return pmc.WorkingSetSize / 2**20, pmc.PeakWorkingSetSize / 2**20
    except Exception:
        return -1.0, -1.0


class Budget(Exception):
    pass


def check_budget(t0, tlimit, mem_mb=380.0):
    if time.time() - t0 > tlimit:
        raise Budget('time limit %.0f s' % tlimit)
    cur, _ = rss_mb()
    if cur > mem_mb:
        raise Budget('memory limit %.0f MB (current %.0f MB)' % (mem_mb, cur))


# ----------------------------------------------------------------------------------------
#  rational rotations
# ----------------------------------------------------------------------------------------
def rot_from_t(t):
    """unit vector (cos phi, sin phi) with tan(phi/2) = t (rational)."""
    t = Fr(t)
    d = 1 + t * t
    return (1 - t * t) / d, 2 * t / d


def rot_near(phi, den=4000):
    """rational rotation whose angle is close to phi (|error| ~ 1/den)."""
    t = Fr(round(math.tan(phi / 2.0) * den), den)
    return rot_from_t(t)


def sin_rat(theta, den=10**6):
    """a rational number close to sin(theta) (used to SPECIFY a threshold angle exactly)."""
    return Fr(round(math.sin(theta) * den), den)


def asin_q(q):
    q = Fr(q)
    return mp.asin(mp.mpf(q.numerator) / q.denominator)


def mpq(q):
    q = Fr(q)
    return mp.mpf(q.numerator) / q.denominator


# ----------------------------------------------------------------------------------------
#  squares and packings
# ----------------------------------------------------------------------------------------
class Sq:
    __slots__ = ('i', 'cx', 'cy', 'c', 's', 'u', 'n', 'V', 'E', 'sa', 'v', 'ext', 'top',
                 'fx0', 'fx1', 'fy0', 'fy1', '_a', 'phi_f')

    def __init__(self, i, cx, cy, c, s):
        cx = Fr(cx); cy = Fr(cy); c = Fr(c); s = Fr(s)
        if c * c + s * s != 1:
            raise ValueError('not a unit vector')
        for (cc, ss) in ((c, s), (-s, c), (-c, -s), (s, -c)):
            if cc > abs(ss):
                break
        else:
            raise ValueError('45 degree square not supported')
        self.i = i; self.cx = cx; self.cy = cy; self.c = cc; self.s = ss
        self.u = (cc, ss); self.n = (-ss, cc)
        hu0, hu1 = cc / 2, ss / 2
        hn0, hn1 = -ss / 2, cc / 2
        V0 = (cx - hu0 - hn0, cy - hu1 - hn1)   # bottom-left  (xi=-1/2, zeta=-1/2)
        V1 = (cx + hu0 - hn0, cy + hu1 - hn1)   # bottom-right (xi=+1/2, zeta=-1/2)
        V2 = (cx + hu0 + hn0, cy + hu1 + hn1)   # top-right
        V3 = (cx - hu0 + hn0, cy - hu1 + hn1)   # top-left
        self.V = (V0, V1, V2, V3)
        # edges: 0 bottom (outward -n), 1 right (u), 2 top (n), 3 left (-u)
        self.E = ((V0, V1, (ss, -cc)), (V1, V2, (cc, ss)), (V2, V3, (-ss, cc)), (V3, V0, (-cc, -ss)))
        self.sa = abs(ss)
        self.v = cy - (cc + self.sa) / 2
        self.ext = cc + self.sa
        self.top = self.v + self.ext
        xs = [p[0] for p in self.V]; ys = [p[1] for p in self.V]
        self.fx0 = float(min(xs)); self.fx1 = float(max(xs))
        self.fy0 = float(min(ys)); self.fy1 = float(max(ys))
        self._a = None
        self.phi_f = math.atan2(float(ss), float(cc))

    @property
    def a(self):
        if self._a is None:
            self._a = asin_q(self.sa)
        return self._a

    def chord(self, y):
        """c_S(y): length of S cap l_y (exact)."""
        t = y - self.v
        if self.sa == 0:
            return F1 if (0 <= t <= 1) else F0
        if t <= 0 or t >= self.ext:
            return F0
        if t < self.sa:
            return t / (self.sa * self.c)
        if t <= self.c:
            return 1 / self.c
        return (self.ext - t) / (self.sa * self.c)

    def chord_lin(self, ym):
        """(slope, intercept) of c_S on the elementary interval containing ym (ym not a breakpoint)."""
        t = ym - self.v
        if self.sa == 0:
            return (F0, F1) if 0 < t < 1 else (F0, F0)
        if t <= 0 or t >= self.ext:
            return (F0, F0)
        sc = self.sa * self.c
        if t < self.sa:
            return (1 / sc, -self.v / sc)
        if t <= self.c:
            return (F0, 1 / self.c)
        return (-1 / sc, (self.v + self.ext) / sc)

    def chord_iv(self, y):
        """closed chord [xl, xr] of S on l_y, or None."""
        if y < self.v or y > self.top:
            return None
        xs = []
        for (A, B, _) in self.E:
            if A[1] == B[1]:
                if A[1] == y:
                    xs.append(A[0]); xs.append(B[0])
                continue
            lo, hi = (A[1], B[1]) if A[1] < B[1] else (B[1], A[1])
            if lo <= y <= hi:
                xs.append(A[0] + (y - A[1]) * (B[0] - A[0]) / (B[1] - A[1]))
        return (min(xs), max(xs))

    def phi_set(self):
        """Phi(S) = {y : 0 < c_S(y) < 1} as a list of open intervals."""
        if self.sa == 0:
            return []
        sc = self.sa * self.c
        return [(self.v, self.v + sc), (self.top - sc, self.top)]

    def breakpoints(self):
        if self.sa == 0:
            return [self.v, self.top]
        sc = self.sa * self.c
        return [self.v, self.v + sc, self.v + self.sa, self.v + self.c, self.top - sc, self.top]


def sat_separated(A, B):
    """exact: True iff the closed squares A, B are disjoint (separating axis among the 4 normals)."""
    for S in (A, B):
        for ax in (S.u, S.n):
            pa = [V[0] * ax[0] + V[1] * ax[1] for V in A.V]
            pb = [V[0] * ax[0] + V[1] * ax[1] for V in B.V]
            if max(pa) < min(pb) or max(pb) < min(pa):
                return True
    return False


class Packing:
    def __init__(self, k, sqs, name=''):
        self.k = Fr(k)
        self.name = name
        self.raw = [tuple(Fr(v) for v in t) for t in sqs]
        self.S = [Sq(i, *t) for i, t in enumerate(self.raw)]
        self.N = len(self.S)
        self.W = self.k * self.k - self.N
        self.grid = {}
        for s in self.S:
            for gx in range(int(math.floor(s.fx0)), int(math.floor(s.fx1)) + 1):
                for gy in range(int(math.floor(s.fy0)), int(math.floor(s.fy1)) + 1):
                    self.grid.setdefault((gx, gy), []).append(s.i)

    def query(self, x0, x1, y0, y1):
        out = set()
        for gx in range(int(math.floor(x0)) - 1, int(math.floor(x1)) + 1):
            for gy in range(int(math.floor(y0)) - 1, int(math.floor(y1)) + 1):
                l = self.grid.get((gx, gy))
                if l:
                    out.update(l)
        return out

    def verify(self):
        """exact check of a CLOSED packing: every square inside [0,k]^2, pairwise disjoint closed squares.
        Returns (ok, message)."""
        k = self.k
        for s in self.S:
            for (vx, vy) in s.V:
                if vx < 0 or vx > k or vy < 0 or vy > k:
                    return False, 'square %d outside the container' % s.i
        npairs = 0
        for s in self.S:
            for j in self.query(s.fx0, s.fx1, s.fy0, s.fy1):
                if j <= s.i:
                    continue
                t = self.S[j]
                if t.fx0 > s.fx1 + 1e-6 or t.fx1 < s.fx0 - 1e-6 or t.fy0 > s.fy1 + 1e-6 or t.fy1 < s.fy0 - 1e-6:
                    continue   # float boxes separated by more than 1e-6: certainly disjoint
                npairs += 1
                if not sat_separated(s, t):
                    return False, 'squares %d and %d intersect (closed sets)' % (s.i, j)
        return True, 'ok (%d close pairs checked exactly)' % npairs

    def reflect(self):
        """rho(x, y) = (x, k - y)."""
        return Packing(self.k, [(s.cx, self.k - s.cy, s.c, -s.s) for s in self.S], name=self.name + '/refl')


# ----------------------------------------------------------------------------------------
#  line structure (omega, E, short chords, max inclination) on elementary intervals
# ----------------------------------------------------------------------------------------
class Lines:
    """Elementary intervals between consecutive breakpoints of all chord functions.
    On each open elementary interval every c_S is affine and its class (0 / short / long) is constant."""

    def __init__(self, P):
        self.P = P
        k = P.k
        bps = {F0, k}
        for s in P.S:
            for b in s.breakpoints():
                if 0 <= b <= k:
                    bps.add(b)
        self.bps = sorted(bps)
        order = sorted(range(P.N), key=lambda i: P.S[i].v)
        vs = [P.S[i].v for i in order]
        self.el = []
        for a, b in zip(self.bps[:-1], self.bps[1:]):
            m = (a + b) / 2
            hiidx = bisect.bisect_left(vs, b)
            om1 = F0; om0 = k
            e1 = F0; e0 = F0
            sh1 = F0; sh0 = F0
            maxsin = F0
            short = []; meet = []; nlong = 0
            for ii in range(hiidx):
                s = P.S[order[ii]]
                if s.top <= a:
                    continue
                if not (s.v < m < s.top):
                    continue
                sl, ic = s.chord_lin(m)
                cm = sl * m + ic
                if cm <= 0:
                    continue
                meet.append(s.i)
                if s.sa > maxsin:
                    maxsin = s.sa
                om1 -= sl; om0 -= ic
                if cm >= 1:
                    nlong += 1
                    e1 += sl; e0 += ic - 1
                else:
                    sh1 += sl; sh0 += ic
                    short.append(s.i)
            self.el.append(dict(lo=a, hi=b, om=(om0, om1), E=(e0, e1), sh=(sh0, sh1), maxsin=maxsin,
                                short=short, meet=meet, nlong=nlong))
        self.lo_list = [e['lo'] for e in self.el]

    def find(self, y):
        i = bisect.bisect_right(self.lo_list, y) - 1
        return self.el[max(0, min(i, len(self.el) - 1))]


def lin_set(lo, hi, c0, c1, rel, thr):
    """sub-interval of (lo, hi) where c0 + c1*y (rel) thr, rel in {'<', '>=', '>'}; returns (a, b) or None."""
    # solve c0 + c1 y = thr
    if c1 == 0:
        val = c0
        ok = (val < thr) if rel == '<' else ((val >= thr) if rel == '>=' else (val > thr))
        return (lo, hi) if ok else None
    r = (thr - c0) / c1
    if rel == '<':
        a, b = ((lo, min(hi, r)) if c1 > 0 else (max(lo, r), hi))
    else:   # '>=' or '>' (same up to measure zero)
        a, b = ((max(lo, r), hi) if c1 > 0 else (lo, min(hi, r)))
    return (a, b) if b > a else None


def H_ivs(lo, hi, tau0):
    """(lo, hi) cap {y : dist(y, Z) >= tau0} as a list of intervals."""
    out = []
    n0 = math.floor(lo) - 1
    n1 = math.floor(hi) + 1
    for n in range(n0, n1 + 1):
        a = max(lo, n + tau0); b = min(hi, n + 1 - tau0)
        if b > a:
            out.append((a, b))
    return out


def iv_inter(A, B):
    """intersection of two sorted lists of disjoint intervals."""
    out = []
    i = j = 0
    while i < len(A) and j < len(B):
        a = max(A[i][0], B[j][0]); b = min(A[i][1], B[j][1])
        if b > a:
            out.append((a, b))
        if A[i][1] < B[j][1]:
            i += 1
        else:
            j += 1
    return out


def iv_measure(ivs):
    return sum((b - a for a, b in ivs), F0)


# ----------------------------------------------------------------------------------------
#  exact piecewise-affine flow tracer
# ----------------------------------------------------------------------------------------
TYPES = ('T', 'D', 'E', 'M', 'W', 'H')


class Flow:
    """F(alpha_F, delta, h_F) of [Q, Definitions 3.9-3.10], with alpha_F given by sinA = sin(alpha_F) (rational).

    recs: list of records (xa, xb, typ, node, T0, T1, info) on open floor intervals (xa, xb):
        termination point = T0 + x*T1 (affine); info = (j, Q0, Q1, gq0, gq1) for M/T (entry data), j for E.
    nodes: ('R', parent, (P0, P1, d, t0, t1, src, prio))  gap segment P(x) + t d, t in [0, t0 + t1 x]
           ('P', parent, (j, Q0, Q1, gq0, gq1))          pass of square j entered at Q(x)
    """

    def __init__(self, P, sinA, delta, hF, tlimit=500.0, mem_mb=380.0):
        self.P = P; self.S = P.S; self.k = P.k
        self.sinA = Fr(sinA); self.delta = Fr(delta); self.hF = Fr(hF)
        self.recs = []
        self.cands = {}
        self.heap = []
        self.inheap = set()
        self.done = set()
        self.anom = []
        self.nR1 = 0; self.nR1multi = 0
        self.t0 = time.time(); self.tlimit = tlimit; self.mem_mb = mem_mb
        self.nray = 0

    def run(self):
        k = self.k
        self._ray(F0, k, (F0, F0), (F1, F0), F0, F0, -1, None)
        while self.heap:
            _, j = heapq.heappop(self.heap)
            self.inheap.discard(j)
            self._square(j)
        # bookkeeping check: the record intervals partition (0, k) up to finitely many points
        tot = sum((r[1] - r[0] for r in self.recs), F0)
        if tot != k:
            self.anom.append(('measure', str(tot)))
        return self

    # --- one ray step on a piece ---
    def _ray(self, xa, xb, P0, P1, g0, g1, src, node):
        if xb <= xa:
            return
        self.nray += 1
        if self.nray % 200 == 0:
            check_budget(self.t0, self.tlimit, self.mem_mb)
        S = self.S; k = self.k; delta = self.delta; hF = self.hF
        if src < 0:
            ux, uy, dx, dy = F1, F0, F0, F1
        else:
            X = S[src]; ux, uy = X.u; dx, dy = X.n
        w0 = P0[0] * ux + P0[1] * uy; w1 = P1[0] * ux + P1[1] * uy
        hs = P0[0] * dx + P0[1] * dy
        if P1[0] * dx + P1[1] * dy != 0 or w1 <= 0:
            self.anom.append(('frame', src))
            return
        Wlo = w0 + w1 * xa; Whi = w0 + w1 * xb
        ga = g0 + g1 * xa; gb = g0 + g1 * xb
        reach = delta - min(ga, gb)
        pa = (P0[0] + P1[0] * xa, P0[1] + P1[1] * xa)
        pb = (P0[0] + P1[0] * xb, P0[1] + P1[1] * xb)
        rf = float(reach) + 1e-9
        xsf = (float(pa[0]), float(pb[0]), float(pa[0]) + rf * float(dx), float(pb[0]) + rf * float(dx))
        ysf = (float(pa[1]), float(pb[1]), float(pa[1]) + rf * float(dy), float(pb[1]) + rf * float(dy))
        cand = self.P.query(min(xsf) - 1e-6, max(xsf) + 1e-6, min(ysf) - 1e-6, max(ysf) + 1e-6)
        edges = []
        for j in cand:
            if j == src:
                continue
            Y = S[j]
            fw = [V[0] * ux + V[1] * uy for V in Y.V]
            if max(fw) <= Wlo or min(fw) >= Whi:
                continue
            fh = [V[0] * dx + V[1] * dy for V in Y.V]
            if min(fh) - hs > reach:
                continue
            for e in range(4):
                m = Y.E[e][2]
                if m[0] * dx + m[1] * dy >= 0:
                    continue
                ia, ib = e, (e + 1) % 4
                w1e, h1e, w2e, h2e = fw[ia], fh[ia], fw[ib], fh[ib]
                if w1e > w2e:
                    w1e, h1e, w2e, h2e = w2e, h2e, w1e, h1e
                if w2e == w1e or w2e <= Wlo or w1e >= Whi:
                    continue
                lo = max(w1e, Wlo); hi = min(w2e, Whi)
                hlo = h1e + (h2e - h1e) * (lo - w1e) / (w2e - w1e)
                hhi = h1e + (h2e - h1e) * (hi - w1e) / (w2e - w1e)
                if hlo < hs and hhi < hs:
                    continue
                if (hlo < hs) != (hhi < hs):
                    self.anom.append(('straddle', src, j))
                edges.append((w1e, h1e, w2e, h2e, j, e))
        bps = {Wlo, Whi}
        for ed in edges:
            if Wlo < ed[0] < Whi:
                bps.add(ed[0])
            if Wlo < ed[2] < Whi:
                bps.add(ed[2])
        bps = sorted(bps)
        segs = []
        for i in range(len(bps) - 1):
            b0, b1 = bps[i], bps[i + 1]
            wm = (b0 + b1) / 2
            best = None; bh = None
            for ed in edges:
                if ed[0] <= b0 and ed[2] >= b1:
                    he = ed[1] + (ed[3] - ed[1]) * (wm - ed[0]) / (ed[2] - ed[0])
                    if bh is None or he < bh:
                        bh = he; best = ed
                    elif he == bh:
                        self.anom.append(('tie-edges', ed[4], best[4]))
            x0 = (b0 - w0) / w1; x1 = (b1 - w0) / w1
            if segs and segs[-1][2] is best:
                segs[-1][1] = x1
            else:
                segs.append([x0, x1, best])
        d = (dx, dy)
        for (x0, x1, ed) in segs:
            Fs = [(0, delta - g0, -g1)]
            if dx > 0:
                Fs.append((1, (k - P0[0]) / dx, -P1[0] / dx))
            elif dx < 0:
                Fs.append((1, -P0[0] / dx, -P1[0] / dx))
            if ed is not None:
                w1e, h1e, w2e, h2e, j, e = ed
                sl = (h2e - h1e) / (w2e - w1e)
                Fs.append((2, h1e + sl * (w0 - w1e) - hs, sl * w1))
            Fs.append((3, (hF - P0[1]) / dy, -P1[1] / dy))
            xs = {x0, x1}
            for i in range(len(Fs)):
                for i2 in range(i + 1, len(Fs)):
                    d1 = Fs[i][2] - Fs[i2][2]
                    if d1 != 0:
                        xc = (Fs[i2][1] - Fs[i][1]) / d1
                        if x0 < xc < x1:
                            xs.add(xc)
            xs = sorted(xs)
            out = []
            for i in range(len(xs) - 1):
                u0, u1 = xs[i], xs[i + 1]
                xm = (u0 + u1) / 2
                vals = [f[1] + f[2] * xm for f in Fs]
                vmin = min(vals)
                ch = None
                for idx in range(len(Fs)):
                    if vals[idx] == vmin:
                        ch = idx
                        break
                if vmin < 0:
                    self.anom.append(('negative-t', src, Fs[ch][0]))
                if out and out[-1][2] == ch:
                    out[-1][1] = u1
                else:
                    out.append([u0, u1, ch])
            for (u0, u1, ch) in out:
                f = Fs[ch]
                self._emit(u0, u1, P0, P1, g0, g1, node, d, f[0], f[1], f[2], ed, src)

    def _emit(self, xa, xb, P0, P1, g0, g1, node, d, prio, t0, t1, ed, src):
        rn = ('R', node, (P0, P1, d, t0, t1, src, prio))
        T0 = (P0[0] + t0 * d[0], P0[1] + t0 * d[1])
        T1 = (P1[0] + t1 * d[0], P1[1] + t1 * d[1])
        if prio == 0:
            self.recs.append((xa, xb, 'D', rn, T0, T1, None))
        elif prio == 1:
            self.recs.append((xa, xb, 'W', rn, T0, T1, None))
        elif prio == 3:
            self.recs.append((xa, xb, 'H', rn, T0, T1, None))
        else:
            j, e = ed[4], ed[5]
            if e == 0:
                if j in self.done:
                    self.anom.append(('late-candidate', j))
                    self.recs.append((xa, xb, 'E', rn, T0, T1, j))
                    return
                self.cands.setdefault(j, []).append((xa, xb, T0, T1, g0 + t0, g1 + t1, rn))
                if j not in self.inheap:
                    heapq.heappush(self.heap, (self.S[j].cy, j))
                    self.inheap.add(j)
            else:
                self.recs.append((xa, xb, 'E', rn, T0, T1, j))

    # --- rule R1 and R2 at a square ---
    def _square(self, j):
        self.done.add(j)
        C = self.cands.pop(j, [])
        if not C:
            return
        Y = self.S[j]
        ux, uy = Y.u; nx, ny = Y.n
        self.nR1 += 1
        its = []
        for ci, (xa, xb, Q0, Q1, gq0, gq1, rn) in enumerate(C):
            xi0 = (Q0[0] - Y.cx) * ux + (Q0[1] - Y.cy) * uy
            xi1 = Q1[0] * ux + Q1[1] * uy
            if xi1 <= 0:
                self.anom.append(('xi1<=0', j))
                continue
            za = xi0 + xi1 * xa; zb = xi0 + xi1 * xb
            if za < -HALF or zb > HALF:
                self.anom.append(('not-relint', j))
            gA = gq0 - gq1 * xi0 / xi1; gB = gq1 / xi1
            xA = -xi0 / xi1; xB = 1 / xi1
            its.append((za, zb, gA, gB, xA, xB, ci, xi0, xi1))
        wins = []; loses = []
        if len(its) == 1:
            it = its[0]; c = C[it[6]]
            wins.append((c, c[0], c[1]))
        else:
            self.nR1multi += 1
            bps = set()
            for it in its:
                bps.add(it[0]); bps.add(it[1])
            for a in range(len(its)):
                for b in range(a + 1, len(its)):
                    A_, B_ = its[a], its[b]
                    lo = max(A_[0], B_[0]); hi = min(A_[1], B_[1])
                    if hi <= lo:
                        continue
                    if A_[3] != B_[3]:
                        zc = (B_[2] - A_[2]) / (A_[3] - B_[3])
                        if lo < zc < hi:
                            bps.add(zc)
                    elif A_[2] == B_[2] and A_[5] != B_[5]:
                        zc = (B_[4] - A_[4]) / (A_[5] - B_[5])
                        if lo < zc < hi:
                            bps.add(zc)
            bps = sorted(bps)
            lab = {it[6]: [] for it in its}
            for i in range(len(bps) - 1):
                b0, b1 = bps[i], bps[i + 1]
                cover = [it for it in its if it[0] <= b0 and it[1] >= b1]
                if not cover:
                    continue
                if len(cover) == 1:
                    lab[cover[0][6]].append((b0, b1, True))
                    continue
                zm = (b0 + b1) / 2
                keys = [((it[2] + it[3] * zm), (it[4] + it[5] * zm), it) for it in cover]
                win = min(keys, key=lambda t: (t[0], t[1]))[2]
                for it in cover:
                    lab[it[6]].append((b0, b1, it is win))
            for it in its:
                c = C[it[6]]
                xi0, xi1 = it[7], it[8]
                ivs = []
                for (b0, b1, wf) in lab[it[6]]:
                    x0 = (b0 - xi0) / xi1; x1 = (b1 - xi0) / xi1
                    if ivs and ivs[-1][2] == wf and ivs[-1][1] == x0:
                        ivs[-1][1] = x1
                    else:
                        ivs.append([x0, x1, wf])
                for (x0, x1, wf) in ivs:
                    if x1 > x0:
                        (wins if wf else loses).append((c, x0, x1))
        for (c, x0, x1) in loses:
            xa, xb, Q0, Q1, gq0, gq1, rn = c
            self.recs.append((x0, x1, 'M', rn, Q0, Q1, (j, Q0, Q1, gq0, gq1)))
        for (c, x0, x1) in wins:
            xa, xb, Q0, Q1, gq0, gq1, rn = c
            if Y.sa >= self.sinA:
                self.recs.append((x0, x1, 'T', rn, Q0, Q1, (j, Q0, Q1, gq0, gq1)))
                continue
            pn = ('P', rn, (j, Q0, Q1, gq0, gq1))
            E0 = (Q0[0] + nx, Q0[1] + ny)
            c0 = E0[1] - self.hF; c1 = Q1[1]
            parts = []
            if c1 == 0:
                parts.append((x0, x1, c0 >= 0))
            else:
                xr = -c0 / c1
                if xr <= x0 or xr >= x1:
                    xm = (x0 + x1) / 2
                    parts.append((x0, x1, c0 + c1 * xm >= 0))
                elif c1 > 0:
                    parts += [(x0, xr, False), (xr, x1, True)]
                else:
                    parts += [(x0, xr, True), (xr, x1, False)]
            for (u0, u1, isH) in parts:
                if u1 <= u0:
                    continue
                if isH:
                    self.recs.append((u0, u1, 'H', pn, E0, Q1, None))
                else:
                    self._ray(u0, u1, E0, Q1, gq0, gq1, j, pn)


def node_chain(node):
    out = []
    while node is not None:
        out.append(node)
        node = node[1]
    out.reverse()
    return out


def path_segments(S, rec):
    """segments of the trajectory of a record, in order:
    ('G', A0, A1, d, t0, t1, src)  gap segment from A(x) = A0 + x A1 along d, length t0 + t1 x
    ('P', j, Q0, Q1)              pass segment [Q(x), Q(x) + n_j]"""
    segs = []
    for nd in node_chain(rec[3]):
        if nd[0] == 'R':
            P0, P1, d, t0, t1, src, prio = nd[2]
            segs.append(('G', P0, P1, d, t0, t1, src))
        else:
            j, Q0, Q1, gq0, gq1 = nd[2]
            segs.append(('P', j, Q0, Q1))
    return segs


def tau_aff(rec, hF):
    """termination height as (c0, c1) (affine), or None for H (tau = h_F by convention)."""
    if rec[2] == 'H':
        return None
    return (rec[4][1], rec[5][1])


def measure_lt(xa, xb, c0, c1, y):
    """|{x in (xa, xb) : c0 + c1 x < y}|."""
    if c1 == 0:
        return (xb - xa) if c0 < y else F0
    r = (y - c0) / c1
    if c1 > 0:
        lo, hi = xa, min(xb, r)
    else:
        lo, hi = max(xa, r), xb
    return (hi - lo) if hi > lo else F0


def open_between(xa, xb, a0, a1, b0, b1, y):
    """open sub-interval of (xa, xb) where a0 + a1 x < y < b0 + b1 x; returns (lo, hi) or None."""
    lo, hi = xa, xb
    # a0 + a1 x < y
    if a1 == 0:
        if not (a0 < y):
            return None
    else:
        r = (y - a0) / a1
        if a1 > 0:
            hi = min(hi, r)
        else:
            lo = max(lo, r)
    # b0 + b1 x > y
    if b1 == 0:
        if not (b0 > y):
            return None
    else:
        r = (y - b0) / b1
        if b1 > 0:
            lo = max(lo, r)
        else:
            hi = min(hi, r)
    return (lo, hi) if hi > lo else None


def line_eval(flow, y, want_images=True):
    """exact line quantities of a flow at height y (0 < y < h_F).
    Returns dict with L (by type), mu (per square), Fgap, gap images [(l, r, dens)], pass images per square."""
    S = flow.S; k = flow.k
    L = {t: F0 for t in TYPES}
    mu = {}
    Fgap = F0
    gim = []
    pim = {}
    for rec in flow.recs:
        xa, xb, typ = rec[0], rec[1], rec[2]
        if typ != 'H':
            L[typ] += measure_lt(xa, xb, rec[4][1], rec[5][1], y)
        for seg in path_segments(S, rec):
            if seg[0] == 'G':
                _, A0, A1, d, t0, t1, src = seg
                a0, a1 = A0[1], A1[1]
                if a0 + a1 * xa >= y and a0 + a1 * xb >= y:
                    break
                b0 = a0 + t0 * d[1]; b1 = a1 + t1 * d[1]
                iv = open_between(xa, xb, a0, a1, b0, b1, y)
                if iv is None:
                    continue
                Fgap += iv[1] - iv[0]
                if want_images:
                    r = d[0] / d[1]
                    X0 = A0[0] + (y - A0[1]) * r
                    X1 = A1[0] - A1[1] * r
                    la = X0 + X1 * iv[0]; lb = X0 + X1 * iv[1]
                    gim.append((min(la, lb), max(la, lb), 1 / abs(X1)))
            else:
                _, j, Q0, Q1 = seg
                Z = S[j]
                a0, a1 = Q0[1], Q1[1]
                if a0 + a1 * xa >= y and a0 + a1 * xb >= y:
                    break
                b0 = a0 + Z.c; b1 = a1
                iv = open_between(xa, xb, a0, a1, b0, b1, y)
                if iv is None:
                    continue
                mu[j] = mu.get(j, F0) + (iv[1] - iv[0])
                if want_images:
                    nx, ny = Z.n
                    r = nx / ny
                    X0 = Q0[0] + (y - Q0[1]) * r
                    X1 = Q1[0] - Q1[1] * r
                    la = X0 + X1 * iv[0]; lb = X0 + X1 * iv[1]
                    pim.setdefault(j, []).append((min(la, lb), max(la, lb), 1 / abs(X1)))
    return dict(L=L, mu=mu, Fgap=Fgap, gim=gim, pim=pim)


def density_excess(ims):
    """for intervals with densities [(l, r, dens)]: returns (int (rho-1)_+, max rho, total mass)."""
    ev = []
    for (l, r, dd) in ims:
        if r > l:
            ev.append((l, dd)); ev.append((r, -dd))
    ev.sort(key=lambda t: t[0])
    rho = F0; prev = None; ov = F0; mx = F0; mass = F0
    i = 0
    while i < len(ev):
        x = ev[i][0]
        if prev is not None and x > prev:
            if rho > 1:
                ov += (rho - 1) * (x - prev)
            mass += rho * (x - prev)
        while i < len(ev) and ev[i][0] == x:
            rho += ev[i][1]; i += 1
        if rho > mx:
            mx = rho
        prev = x
    return ov, mx, mass


def line_full(flow, y, check=True):
    """line quantities + exact consistency checks (shadow identity, Bd = 0, Ov = 0, images in the waste)."""
    P = flow.P; S = flow.S; k = flow.k
    le = line_eval(flow, y, want_images=check)
    chords = {}
    cv = []
    for s in S:
        if s.v <= y <= s.top:
            c = s.chord(y)
            if c > 0:
                chords[s.i] = c
                if check:
                    cv.append(s.chord_iv(y) + (s.i,))
    sumc = sum(chords.values(), F0)
    omega = k - sumc
    mu = le['mu']
    Sh = F0; Ov = F0
    for i, c in chords.items():
        m = mu.get(i, F0)
        if c > m:
            Sh += c - m
        elif m > c:
            Ov += m - c
    for i, m in mu.items():
        if i not in chords and m > 0:
            Ov += m
    L = le['L']
    Lsum = sum((L[t] for t in ('T', 'D', 'E', 'M', 'W')), F0)
    musum = sum(mu.values(), F0)
    Bd = k - Lsum - musum - le['Fgap']
    Ovgap, rhomax, gmass = density_excess(le['gim']) if check else (None, None, None)
    out = dict(y=y, L=L, Lsum=Lsum, Sh=Sh, Ov=Ov, omega=omega, Fgap=le['Fgap'], Bd=Bd, Ovgap=Ovgap,
               mu=mu, chords=chords)
    if check:
        bad = []
        # gap crossers lie in the waste: their images avoid the interiors of all chords
        cv.sort()
        starts = [c[0] for c in cv]
        for (l, r, dd) in le['gim']:
            i0 = bisect.bisect_left(starts, r)
            for ii in range(max(0, i0 - 3), i0):
                cl, cr, ci = cv[ii]
                if cr > l and cl < r:
                    bad.append(('gap-in-square', ci))
                    break
        # pass crossers of Z lie in the chord of Z and have density <= 1
        for j, ims in le['pim'].items():
            civ = S[j].chord_iv(y)
            for (l, r, dd) in ims:
                if civ is None or l < civ[0] or r > civ[1]:
                    bad.append(('pass-outside-chord', j))
            ovp, mxp, _ = density_excess(ims)
            if mxp > 1:
                bad.append(('pass-density>1', j))
        if gmass != le['Fgap']:
            bad.append(('gap-mass', str(gmass - le['Fgap'])))
        # shadow identity (c) with Ov = 0 and Bd = 0 off Exc
        ident = (Sh - Ov) - (Lsum + le['Fgap'] - omega + Bd)
        if ident != 0:
            bad.append(('shadow-identity', str(ident)))
        out['bad'] = bad
        out['rhomax'] = rhomax
        out['U'] = Ovgap - le['Fgap'] + omega
    return out
