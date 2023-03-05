import mult
import sqrt

for n in range(9):
  print(n)
  d = (-2,3,-5,7,-11,13,-17,19)[:n]
  N = 2^n

  for loop in range(10):

    for bits in range(200):
      S = set(j for j in range(n) if randrange(2))
      g = randrange(1-2^bits,2^bits)
      h = g^2 * prod(d[j] for j in S)
      if h == 0:
        assert sqrt.twisted(n,d,h) == (0,set())
      else:
        assert sqrt.twisted(n,d,h) == (abs(g),S)
  
    for bits in range(200):
      h = randrange(1-2^bits,2^bits)
      try:
        g,S = sqrt.twisted(n,d,h)
        assert all(j in range(n) for j in S)
        assert h == g^2 * prod(d[j] for j in S)
      except:
        pass
  
    for bits in range(200):
      for m in range(n+1):
        M = 2^m
        f = [randrange(1-2^bits,2^bits) for j in range(M)]
        ff = mult.square(m,M,d[:m],f)
        R = set(j for j in range(m,n) if randrange(2))
        pi = prod(d[j] for j in R)
        h = tuple(ffj*pi for ffj in ff)
        g,S = sqrt.squareroot(n,d,m,h)
        assert(len(g) == M)
        if all(gj == 0 for gj in g):
          assert(all(hj == 0) for hj in h)
          assert(S == set())
        else:
          assert(S == R)
          assert [g[j] == f[j] for j in range(M)] or [g[j] == -f[j] for j in range(M)]
