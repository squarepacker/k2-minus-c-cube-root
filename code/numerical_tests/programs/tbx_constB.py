# -*- coding: utf-8 -*-
"""Side check (not one of the seven items): the constants quoted in Theorem B of the first draft of the argument, with interval
arithmetic (mpmath.iv, 60 digits).
  A   = 4 sqrt(K/Q*) / (2 - 3e-6),  K = 9 (Lemma 4.10) or 13 (Lemma 4.9), Q* = 0.32
  A'  = A / cos(beta_bar) + 1e-8,  beta_bar = 1.5e-6
  B1  = (A/2 + 2 beta_bar) / cos(beta_bar)
  c_inf^3 = sup over tau1 < tau0 < 1/2, nu0 > 0 of 2 tau1 (1-nu0)(1-2 tau0)^2 / (A A' B1)
          = (4/27) / (A A' B1)      (tau1 -> tau0 = 1/6, nu0 -> 0)
  draft:  'With B1 ~ A/2 ... c_inf ~ (0.2963/(A^2 A'))^(1/3) ~ 0.0627'; 'Without Lemma 4.10 it is ~ 0.0522'.
"""
import mpmath as mp
from mpmath import iv

iv.dps = 60
mp.mp.dps = 60


def consts(K):
    bb = iv.mpf('1.5e-6')
    A = 4 * iv.sqrt(iv.mpf(K) / iv.mpf('0.32')) / (2 - iv.mpf('3e-6'))
    Ap = A / iv.cos(bb) + iv.mpf('1e-8')
    B1 = (A / 2 + 2 * bb) / iv.cos(bb)
    c3 = iv.mpf(4) / 27 / (A * Ap * B1)
    c = c3 ** (iv.mpf(1) / 3)
    c3_approx = iv.mpf('0.2963') / (A * A * Ap)
    c_approx = c3_approx ** (iv.mpf(1) / 3)
    return A, Ap, B1, c, c_approx


# numerical check that the sup of 2 t (1-2 t)^2 over t in (0, 1/2) is 4/27 at t = 1/6
f = lambda t: 2 * t * (1 - 2 * t) ** 2
grid_max = max(f(mp.mpf(i) / 100000) for i in range(1, 50000))
print('max_t 2t(1-2t)^2 on a grid: %s   4/27 = %s' % (mp.nstr(grid_max, 12), mp.nstr(mp.mpf(4) / 27, 12)))
for K, name in ((9, 'K_w*=9 (with [R, Lemma 4.10])'), (13, 'K_w=13 (analytic [R, Lemma 4.9])')):
    A, Ap, B1, c, ca = consts(K)
    print(name)
    print('  A   in', A)
    print('  A\'  in', Ap)
    print('  B1  in', B1, '   (A/2 =', A / 2, ')')
    print('  c_inf = ((4/27)/(A A\' B1))^(1/3) in', c)
    print('  draft formula (0.2963/(A^2 A\'))^(1/3) in', ca)
