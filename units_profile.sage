import sys
import nprofile
import units
import goodprime
import mult
import div
import subsetprod
import powerprod
import sqrt

n = 5
if len(sys.argv) > 1: n = ZZ(sys.argv[1])

start = 0
if len(sys.argv) > 2: start = ZZ(sys.argv[2])

d = ()
for p in primes(start,infinity):
  if len(d) >= n: break
  d += (p,)

G = units.generators(d)
for u in G: print u

nprofile.output([goodprime,mult,div,subsetprod,powerprod,sqrt,units])
