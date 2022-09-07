import nprofile
import goodprime
import mult
import div

for n in range(10):
  d = (2,3,5,7,11,13,17,19,23)[:n]
  N = 2^n

  bits = 1000
  f = [randrange(1-2^bits,2^bits) for j in range(N)]
  g = [randrange(1-2^bits,2^bits) for j in range(N)]
  h = mult.mult(n,N,d,f,g)
  for loop in range(20):
    f = div.div(n,N,d,h,g)

nprofile.output([goodprime,mult,div])
