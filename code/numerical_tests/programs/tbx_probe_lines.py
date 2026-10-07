# -*- coding: utf-8 -*-
"""probe: line structure of a generated packing (measure of low-waste lines, E there, inclinations)."""
import sys, random
from fractions import Fraction as Fr
import tbx_core as C
import tbx_gen as G
import tbx_checks as K

fam = sys.argv[1]; k = int(sys.argv[2]); seed = int(sys.argv[3])
kw = eval(sys.argv[4]) if len(sys.argv) > 4 else {}
rng = random.Random(seed)
P = G.two_sided(fam, rng, k, kw)
print(fam, 'k', k, 'N', P.N, 'W', P.W, P.verify())
LD = K.LineData(P)
for nu in (Fr(3, 10), Fr(49, 100), Fr(3, 4), Fr(1)):
    m = Fr(0); Emax = None; Emin = None; msmax = Fr(0); segs = []
    for e in LD.el:
        iv = C.lin_set(e['lo'], e['hi'], e['om'][0], e['om'][1], '<', nu)
        if iv:
            m += iv[1] - iv[0]
            E0, E1 = e['E']
            for yy in iv:
                Ev = E0 + E1 * yy
                Emax = Ev if Emax is None or Ev > Emax else Emax
                Emin = Ev if Emin is None or Ev < Emin else Emin
            msmax = max(msmax, e['maxsin'])
            segs.append((float(iv[0]), float(iv[1])))
    print('omega <', float(nu), ': measure %.4f' % float(m), 'E range', None if Emax is None else (float(Emin), float(Emax)),
          'max sin a %.4f' % float(msmax), 'pieces', len(segs), segs[:6])
