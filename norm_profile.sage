import nprofile
import goodprime
import mult
import norm

for n in range(8):
  d = (2,3,5,7,11,13,17,19,23)[:n]
  N = 2^n

  bits = 1000
  f = [randrange(1-2^bits,2^bits) for j in range(N)]
  for loop in range(10):
    for t in range(n):
      norm.normonestep(n,N,d,f,t)

nprofile.output([goodprime,mult,norm])
