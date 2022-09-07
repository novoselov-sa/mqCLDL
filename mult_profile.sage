import nprofile
import goodprime
import mult

for n in range(10):
  d = (2,3,5,7,11,13,17,19,23)[:n]
  N = 2^n

  bits = 10
  f = [randrange(1-2^bits,2^bits) for j in range(N)]
  g = [randrange(1-2^bits,2^bits) for j in range(N)]
  for loop in range(100):
    h = mult.mult(n,N,d,f,g)

nprofile.output([goodprime,mult])
