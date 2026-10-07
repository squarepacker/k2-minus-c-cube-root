# -*- coding: utf-8 -*-
"""tbx_gen.py -- generators of exact rational closed packings for the first-draft test.

Every proposed square is accepted only if it is (exactly) inside [0,k] x [0,ymax] and disjoint, as a closed
set, from every square accepted before (exact separating-axis test); the finished packing is verified again by
Packing.verify().  All families are built TWO-SIDED: a bottom half [0,k] x [0, k/2 - 1/100] is generated from
the floor upwards, a second independent bottom half is generated and reflected by (x,y) -> (x, k-y), so that the
floor flow and the ceiling flow are both non-trivial (the flows only reach height h = k/2 - 1).

Families (item 6 of the task):
  column   (a) tilted columns: squares at angle a stacked along n through which many paths pass
  phase    (b) stacks of tilted 'phase' squares shifting the squares above off the integers
  lshape   (c) grid of slightly tilted squares plus tilted 1 x b blocks in strips along two walls
  random   (d) random column-stacked packings (fractional heights, random small rotations, deletions)
  rows     (d') random row packings
  rotgrid  rotated grid (full-width low-waste bands; used for item 5 with large beta)
  stair    deficit staircase (negative control): a steep column feeding a slightly tilted square
"""
from fractions import Fraction as Fr
import math, random
from tbx_core import Sq, Packing, sat_separated, rot_near, rot_from_t, F0, F1, HALF


class Placer:
    def __init__(self, k, ymax=None):
        self.k = Fr(k)
        self.ymax = Fr(k) if ymax is None else Fr(ymax)
        self.sqs = []
        self.raw = []
        self.grid = {}

    def _cells(self, s):
        for gx in range(int(math.floor(s.fx0)), int(math.floor(s.fx1)) + 1):
            for gy in range(int(math.floor(s.fy0)), int(math.floor(s.fy1)) + 1):
                yield (gx, gy)

    def ok(self, cx, cy, c, s):
        k = self.k
        q = Sq(-1, cx, cy, c, s)
        for (vx, vy) in q.V:
            if vx < 0 or vx > k or vy < 0 or vy > self.ymax:
                return None
        seen = set()
        for cell in self._cells(q):
            for gx in (cell[0] - 1, cell[0], cell[0] + 1):
                for gy in (cell[1] - 1, cell[1], cell[1] + 1):
                    for i in self.grid.get((gx, gy), ()):
                        if i in seen:
                            continue
                        seen.add(i)
                        t = self.sqs[i]
                        if t.fx0 > q.fx1 + 1e-6 or t.fx1 < q.fx0 - 1e-6 or t.fy0 > q.fy1 + 1e-6 or t.fy1 < q.fy0 - 1e-6:
                            continue
                        if not sat_separated(q, t):
                            return None
        return q

    def add(self, cx, cy, c, s):
        q = self.ok(Fr(cx), Fr(cy), Fr(c), Fr(s))
        if q is None:
            return False
        i = len(self.sqs)
        q2 = Sq(i, q.cx, q.cy, q.c, q.s)
        self.sqs.append(q2)
        self.raw.append((q.cx, q.cy, q.c, q.s))
        for cell in self._cells(q2):
            self.grid.setdefault(cell, []).append(i)
        return True


def frq(x, den=10**4):
    return Fr(round(x * den), den)


def rot(phi, den=4000):
    if phi == 0:
        return (F1, F0)
    return rot_near(phi, den)


def bbox_half(c, s):
    return (abs(c) + abs(s)) / 2


def fill_grid(pl, x0, x1, y0, y1, rng, angles=(0,), pdel=0.0, gap=Fr(1, 500), jitter=0.0):
    """fill [x0,x1] x [y0,y1] with rows of squares in slots (proposals that collide are skipped)."""
    y = Fr(y0)
    while True:
        row = []
        x = Fr(x0)
        rowh = F0
        while True:
            phi = rng.choice(angles)
            c, s = rot(phi)
            hb = bbox_half(c, s)
            if x + 2 * hb > x1:
                break
            row.append((x + hb, hb, c, s))
            x += 2 * hb + gap
            rowh = max(rowh, 2 * hb)
        if not row or y + rowh > y1:
            break
        for (cx, hb, c, s) in row:
            if rng.random() < pdel:
                continue
            pl.add(cx, y + hb, c, s)
        y += rowh + gap


def fill_columns(pl, x0, x1, y0, y1, rng, angles=(0,), gap=(0.0005, 0.01), wslot=None):
    """fill [x0,x1] x [y0,y1] column-wise: squares stacked in unit-ish slots with small vertical gaps."""
    amax = max(abs(a) for a in angles)
    if wslot is None:
        wslot = Fr(math.cos(amax) + math.sin(amax)).limit_denominator(10**4) + Fr(1, 500)
    x = Fr(x0)
    while x + wslot <= x1:
        y = Fr(y0)
        while True:
            phi = rng.choice(angles)
            c, s = rot(phi)
            hb = bbox_half(c, s)
            if y + 2 * hb > y1:
                break
            pl.add(x + wslot / 2, y + hb, c, s)
            y += 2 * hb + frq(rng.uniform(*gap), 10**5)
        x += wslot + Fr(1, 300)


# ---------------------------------------------------------------------------------------------
#  one-sided generators (bottom part, squares below ymax); each returns (raw squares, meta)
# ---------------------------------------------------------------------------------------------
def half_random(rng, k, ymax, amax=0.06, pdel=0.08, pbig=0.15, ptilt=0.5):
    pl = Placer(k, ymax)
    wslot = Fr(math.cos(amax) + math.sin(amax)).limit_denominator(10**4) + Fr(1, 500)
    x = frq(rng.uniform(0, 0.05))
    while x + wslot <= k:
        y = frq(rng.uniform(0, 0.004), 10**5)
        while True:
            if rng.random() < ptilt:
                phi = rng.choice([rng.uniform(-amax, amax), rng.uniform(-amax / 3, amax / 3)])
            else:
                phi = 0
            c, s = rot(phi)
            hb = bbox_half(c, s)
            if y + 2 * hb > ymax:
                break
            jx = frq(rng.uniform(0, min(0.01, float(wslot - 2 * hb))))
            if rng.random() >= pdel:
                pl.add(x + hb + jx, y + hb, c, s)
            gap = rng.uniform(0.05, 0.35) if rng.random() < pbig else rng.uniform(0.0005, 0.015)
            y += 2 * hb + frq(gap, 10**5)
        x += wslot + frq(rng.uniform(0.001, 0.02))
    return pl.raw, dict(amax=amax, pdel=pdel, pbig=pbig, ptilt=ptilt)


def half_rows(rng, k, ymax, amax=0.04, pdel=0.08, pbig=0.2):
    pl = Placer(k, ymax)
    y = frq(rng.uniform(0, 0.004), 10**5)
    first = True
    while True:
        x = frq(rng.uniform(0, 0.4))
        row = []
        rowh = F0
        while True:
            phi = rng.choice([0, 0, 0, rng.uniform(-amax, amax), rng.uniform(-amax / 4, amax / 4)])
            c, s = rot(phi)
            hb = bbox_half(c, s)
            if x + 2 * hb > k:
                break
            row.append((x + hb, hb, c, s))
            x += 2 * hb + frq(rng.uniform(0.002, 0.05))
            rowh = max(rowh, 2 * hb)
        if not row or y + rowh > ymax:
            break
        for (cx, hb, c, s) in row:
            if rng.random() < pdel:
                continue
            jy = F0 if first else frq(rng.uniform(0, min(0.02, float(rowh - 2 * hb))), 10**5)
            pl.add(cx, y + hb + jy, c, s)
        first = False
        gap = rng.uniform(0.05, 0.35) if rng.random() < pbig else rng.uniform(0.0005, 0.015)
        y += rowh + frq(gap, 10**5)
    return pl.raw, dict(amax=amax, pdel=pdel)


def half_column(rng, k, ymax, a=0.035, ncols=1, g=Fr(1, 2000), base=Fr(1, 1000), fill_angles=(0,)):
    """(a) ncols columns of squares at angle a stacked along n with gap g, the lowest vertex at height base;
    the rest filled column-wise with squares at angles fill_angles."""
    pl = Placer(k, ymax)
    c, s = rot(a)
    sa = abs(s)
    xs = sorted(rng.sample(range(1, int(k) - 2), ncols))
    for xc0 in xs:
        cx = Fr(xc0) + HALF + frq(rng.uniform(0, 0.3))
        cy = base + (c + sa) / 2
        n = (-s, c)
        i = 0
        while pl.add(cx + i * (1 + g) * n[0], cy + i * (1 + g) * n[1], c, s):
            i += 1
    fill_columns(pl, F0, Fr(k), Fr(1, 2000), Fr(ymax), rng, angles=fill_angles)
    return pl.raw, dict(a=a, ncols=ncols)


def half_phase(rng, k, ymax, b=0.12, abeta=0.02, gap=Fr(1, 400)):
    """(b) in each column slot, p tilted 'phase' squares (angle +-b) stacked vertically, then slightly tilted
    squares (|angle| <= abeta) on top: their heights are shifted by about p (cos b + sin b - 1)."""
    pl = Placer(k, ymax)
    cb, sb = rot(b)
    hb_b = bbox_half(cb, sb)
    x = frq(rng.uniform(0, 0.1))
    while True:
        p = rng.randint(0, 4)
        wslot = 2 * hb_b
        if x + wslot > k:
            break
        y = frq(rng.uniform(0, 0.002), 10**5)
        for i in range(p):
            sgn = rng.choice([1, -1])
            pl.add(x + hb_b, y + hb_b, cb, sgn * sb)
            y += 2 * hb_b + gap
        while True:
            phi = rng.uniform(-abeta, abeta)
            c, s = rot(phi)
            hb = bbox_half(c, s)
            if y + 2 * hb > ymax:
                break
            pl.add(x + hb_b, y + hb, c, s)
            y += 2 * hb + gap
        x += wslot + frq(rng.uniform(0.002, 0.03))
    return pl.raw, dict(b=b, abeta=abeta)


def half_lshape(rng, k, ymax, w=1.25, a=0.05, abeta=0.015, nb=None):
    """(c) strips of width w along the floor and along the left wall filled with tilted 1 x b blocks
    (b squares in a row along u at angle a / in a column along n at angle -a); the rest a grid of slightly
    tilted squares (|angle| <= abeta) starting at (w, w): heights off the integers by w."""
    pl = Placer(k, ymax)
    w = frq(w)
    c, s = rot(a)
    sa = abs(s)
    g = Fr(1, 1000)
    if nb is None:
        nb = max(1, int(math.floor((float(w) - float(c) - 0.01) / float(sa) + 1e-9)))
    u = (c, s); n = (-s, c)
    x = w + Fr(1, 100)
    while True:
        base = Fr(1, 2000)
        cx = x + u[0] / 2 + n[0] / 2
        cy = base + u[1] / 2 + n[1] / 2
        placed = 0
        for i in range(nb):
            if pl.add(cx + i * (1 + g) * u[0], cy + i * (1 + g) * u[1], c, s):
                placed += 1
            else:
                break
        if placed == 0:
            break
        x += (nb * (1 + g)) * u[0] + sa + Fr(1, 100)
        if x >= k:
            break
    cl, sl = rot(-a)
    nl = (-sl, cl)
    ext = abs(cl) + abs(sl)
    y = w + Fr(1, 100)
    while True:
        cx = Fr(1, 2000) + ext / 2
        cy = y + ext / 2
        placed = 0
        for i in range(nb):
            if pl.add(cx + i * (1 + g) * nl[0], cy + i * (1 + g) * nl[1], cl, sl):
                placed += 1
            else:
                break
        if placed == 0:
            break
        y += (nb - 1) * (1 + g) * nl[1] + ext + Fr(1, 50)
        if y >= ymax:
            break
    fill_columns(pl, w + Fr(1, 50), Fr(k), w + Fr(1, 50), Fr(ymax), rng,
                 angles=(0, abeta, -abeta, abeta / 2), gap=(0.0005, 0.004))
    fill_grid(pl, F0, w, F0, w, rng, angles=(0,), gap=Fr(1, 400))
    return pl.raw, dict(w=float(w), a=a, abeta=abeta, nb=nb)


def half_rotgrid(rng, k, ymax, a=0.33, gap=Fr(1, 1000), oy=None):
    """rotated grid at angle a placed on the floor (lowest vertex at height ~0), as many squares as fit;
    rest filled column-wise axis-parallel.  oy shifts the grid vertically."""
    pl = Placer(k, ymax)
    c, s = rot(a)
    u = (c, s); n = (-s, c)
    sa = abs(s)
    if oy is None:
        oy = Fr(1, 2000)
    ncol = int(math.floor((float(k) - 0.01) / (float(c) * 1.0012))) + 2
    ox = Fr(1, 200) + sa * int(k)
    for r in range(int(k) + 2):
        for q in range(-int(k), ncol + 2):
            cx = ox + (q + HALF) * (1 + gap) * u[0] + (r + HALF) * (1 + gap) * n[0]
            cy = oy + (q + HALF) * (1 + gap) * u[1] + (r + HALF) * (1 + gap) * n[1]
            pl.add(cx, cy, c, s)
    fill_columns(pl, F0, Fr(k), Fr(1, 2000), Fr(ymax), rng, angles=(0,))
    return pl.raw, dict(a=a)


def half_stair(rng, k, ymax, a=0.34, m=6, ab=0.02, g=Fr(1, 3000), eps=Fr(1, 5000)):
    """deficit staircase (negative control): m squares at angle a > 0 stacked along n (lowest vertex at height
    ~0; the floor paths that enter near that vertex climb the whole column and leave its top side near the
    top-left vertex V3), then a 'target' square tilted by ab (0 < ab << a) whose bottom-right vertex lies at
    horizontal distance `off` right of V3, eps above the top side of the column; `off` in [0.02, 0.25] is chosen
    so that the target's bottom height is as far as possible from the integers.  The climbing paths enter the
    target with deficit D = m (1 - cos a)."""
    pl = Placer(k, ymax)
    c, s = rot(a)
    sa = abs(s)
    n = (-s, c)
    cx0 = Fr(int(k) // 2) + Fr(m) * sa / 2 + 1
    cy0 = Fr(1, 4000) + (c + sa) / 2
    for i in range(m):
        pl.add(cx0 + i * (1 + g) * n[0], cy0 + i * (1 + g) * n[1], c, s)
    tc = (cx0 + (m - 1) * (1 + g) * n[0], cy0 + (m - 1) * (1 + g) * n[1])
    V3 = (tc[0] - c / 2 + n[0] / 2, tc[1] - s / 2 + n[1] / 2)        # top-left vertex (lowest point of the top side)
    ta = s / c
    best = None
    for i in range(24):
        off = Fr(2, 100) + Fr(i, 100)
        yr = V3[1] + off * ta + eps
        f = yr - math.floor(yr)
        dz = min(f, 1 - f)
        if best is None or dz > best[0]:
            best = (dz, off)
    off = best[1]
    ct, st = rot(ab)
    ut = (ct, st); nt = (-st, ct)
    BR = (V3[0] + off, V3[1] + off * ta + eps)                          # bottom-right vertex of the target
    tcx = BR[0] - ut[0] / 2 + nt[0] / 2
    tcy = BR[1] - ut[1] / 2 + nt[1] / 2
    placed = pl.add(tcx, tcy, ct, st)
    fill_columns(pl, F0, Fr(k), Fr(1, 2000), Fr(ymax), rng, angles=(0,))
    return pl.raw, dict(a=a, m=m, ab=ab, off=float(off), target_dist=float(best[0]), target_placed=placed)


def half_converge(rng, k, ymax, a=0.08, gapx=Fr(1, 2000), eps=Fr(1, 4000), pskip=0.1):
    """merge stress: rows of V-pairs (left square at angle -a, right square at +a, horizontal gap gapx between
    their bounding boxes: the paths leaving them converge towards the notch), each notch bridged by an
    axis-parallel square just above; then the next row.  Paths from the two sides of a notch reach the same
    entry points of the bridge, and R1 merges them."""
    pl = Placer(k, ymax)
    cL, sL = rot(-a)
    cR, sR = rot(a)
    hb = bbox_half(cR, sR)
    ta = abs(sR) / cR
    y = frq(rng.uniform(0, 0.002), 10**5)
    while True:
        if y + 2 * hb + 1 + eps + HALF * ta + Fr(1, 100) > ymax:
            break
        x = frq(rng.uniform(0, 0.3))
        bridges = []
        while x + 4 * hb + gapx <= k:
            if rng.random() < pskip:
                x += 2 * hb + gapx
                continue
            cxL = x + hb; cxR = x + 3 * hb + gapx
            pl.add(cxL, y + hb, cL, sL)
            pl.add(cxR, y + hb, cR, sR)
            xn = x + 2 * hb + gapx / 2
            bridges.append(xn)
            x += 4 * hb + gapx + frq(rng.uniform(0.001, 0.01))
        # bridge squares: bottom at notch height + (1/2) tan a + eps clears both top sides
        yb = y + cR + HALF * ta + eps
        for xn in bridges:
            pl.add(xn, yb + HALF, F1, F0)
        y = yb + 1 + frq(rng.uniform(0.0005, 0.004), 10**5)
    return pl.raw, dict(a=a)


def half_wall(rng, k, ymax, a=0.1, eps=Fr(1, 2000), vgap=(0.0005, 0.01)):
    """wall stress: a column along each side wall of squares tilted TOWARDS the wall (left wall +a, right wall -a),
    each touching the wall up to eps, stacked vertically; paths leaving their top sides near the wall run into
    the wall (W terminations).  The interior is filled column-wise with squares of angle 0 or +-a/4."""
    pl = Placer(k, ymax)
    for side in (0, 1):
        phi = a if side == 0 else -a
        c, s = rot(phi)
        hb = bbox_half(c, s)
        cx = eps + hb if side == 0 else Fr(k) - eps - hb
        y = frq(rng.uniform(0, 0.002), 10**5)
        while y + 2 * hb <= ymax:
            pl.add(cx, y + hb, c, s)
            y += 2 * hb + frq(rng.uniform(*vgap), 10**5)
    c, s = rot(a)
    w0 = 2 * bbox_half(c, s) + 2 * eps + Fr(1, 200)
    fill_columns(pl, w0, Fr(k) - w0, Fr(1, 2000), Fr(ymax), rng, angles=(0, a / 4, -a / 4))
    return pl.raw, dict(a=a)


HALVES = dict(random=half_random, rows=half_rows, column=half_column, phase=half_phase,
              lshape=half_lshape, rotgrid=half_rotgrid, stair=half_stair, converge=half_converge,
              wall=half_wall)


def two_sided(fam, rng, k, kw_bot, kw_top=None):
    """bottom half from the floor + reflected independent bottom half (the top half)."""
    fn = HALVES[fam]
    ymax = Fr(k, 2) - Fr(1, 100)
    rb, mb = fn(rng, k, ymax, **kw_bot)
    rt, mt = fn(rng, k, ymax, **(kw_top if kw_top is not None else kw_bot))
    raw = list(rb) + [(cx, Fr(k) - cy, c, -s) for (cx, cy, c, s) in rt]
    P = Packing(k, raw, name=fam)
    ok, msg = P.verify()
    if not ok:
        raise RuntimeError('generator produced an invalid packing: ' + msg)
    P.meta = dict(bottom=mb, top=mt)
    return P
