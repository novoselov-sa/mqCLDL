import nprofile
import goodprime
import mult
import div
import sqrt

n = 3

d = (2,3,5,7,11,13,17,19)[:n]
N = 2^n

bits = 1000
f = [randrange(1-2^bits,2^bits) for j in range(N)]
h = mult.mult(n,N,d,f,f)
for loop in range(100):
  g,S = sqrt.squareroot(n,d,n,h)

nprofile.output([goodprime,mult,div,sqrt])
