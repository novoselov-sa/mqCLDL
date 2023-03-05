import sys
import nprofile
import units
import goodprime
import mult
import div
import subsetprod
import powerprod
import sqrt
import field
import classgroup
import cryptosystem
import ideal

n = 8
if len(sys.argv) > 1: n = ZZ(sys.argv[1])
N = 2^n

start = 0
if len(sys.argv) > 2: start = ZZ(sys.argv[2])

loops = 1
if len(sys.argv) > 3: loops = ZZ(sys.argv[3])

d = ()
for p in primes(start,infinity):
  if len(d) >= n: break
  d += (p,)

d = (5,13,17,29,37,41,53,61)
params = cryptosystem.parameters(d)

for loop in range(loops):
  K = field.field(d)
  test = cryptosystem.semi_random_element(d,params)
  print test
  print K(test,1).conj()
  print classgroup.mainwork(n,N,d,12)
  sys.stdout.flush()


nprofile.output([goodprime,mult,div,subsetprod,powerprod,sqrt,cryptosystem,units,ideal])
