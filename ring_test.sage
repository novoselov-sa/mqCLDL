from ring import *

for n in range(0,7):
  print n

  N = 2^n
  d = (-2,3,-5,7,-11,13,-17,19,-23,29)[:n]
  R = ring(d)

  P = PolynomialRing(ZZ,'x',n)
  x = P.gens()
  Q = P.quotient([x[i]^2 - d[i] for i in range(n)])

  for bits in [0,1,2,4,8,16,32,64,128,256]:
    f = R(ZZ.random_element(1-2^bits,2^bits) for j in range(N))
    g = R(ZZ.random_element(1-2^bits,2^bits) for j in range(N))
    assert f - g == f + (-g)
    assert (f + g) - g == f
    h = f * g
    assert (h,) == subsetprod(d,(f,g),((1,1),))
    assert (h,h) == subsetprod(d,(f,f,g,g,h,h),((1,0,0,1,0,0),(0,1,1,0,0,0)))
    if g.is_nonzero():
      assert h.divexact(g) == f
      assert h / g == f
      scale = R((N,) + (0,)*(N-1))
      assert (f,g,h) == scaledpowerprod(d,(scale,g,h),(1,1,1),((-1,-1,1),(-1,1,0),(-1,0,1)),1)

    for t in range(10):
      assert h.symbol(t) == f.symbol(t) * g.symbol(t)
    assert h.symbols(2,20) == tuple(h.symbol(t) for t in range(2,20))

    h = f.square()
    assert h == f * f

    s = h.sqrt()
    assert s == f or s == -f

    Pf = P(sum(f.c[j] * prod(x[i] for i in range(n) if j & (1 << i)) for j in range(N)))
    Qf = Q(Pf)
    Qg = Q(sum(g.c[j] * prod(x[i] for i in range(n) if j & (1 << i)) for j in range(N)))

    h = f + g
    Qh = Q(sum(h.c[j] * prod(x[i] for i in range(n) if j & (1 << i)) for j in range(N)))
    assert Qh == Qf + Qg

    h = f - g
    Qh = Q(sum(h.c[j] * prod(x[i] for i in range(n) if j & (1 << i)) for j in range(N)))
    assert Qh == Qf - Qg

    h = -f
    Qh = Q(sum(h.c[j] * prod(x[i] for i in range(n) if j & (1 << i)) for j in range(N)))
    assert Qh == -Qf

    h = f * g
    Qh = Q(sum(h.c[j] * prod(x[i] for i in range(n) if j & (1 << i)) for j in range(N)))
    assert Qh == Qf * Qg

    h = f.square()
    Qh = Q(sum(h.c[j] * prod(x[i] for i in range(n) if j & (1 << i)) for j in range(N)))
    assert Qh == Qf^2

    for t in range(n):
      h = f.normonestep(t)
      z = [0 if (j&(1<<t)) else h.c[(j % (1<<t)) + ((j >> (t+1)) << t)] for j in range(N)]
      Qh = Q(sum(z[j] * prod(x[i] for i in range(n) if j & (1 << i)) for j in range(N)))
      assert Qh == Qf*Pf.substitute({x[t]:-x[t]})

    h = f.quadnorms()
    for t in range(n):
      z = [0 for j in range(N)]
      z[0] = h[t].c[0]
      z[2^t] = h[t].c[1]
      Qh = Q(sum(z[j] * prod(x[i] for i in range(n) if j & (1 << i)) for j in range(N)))
      y = Pf
      for i in range(n):
        if i != t:
	  y = y*y.substitute({x[i]:-x[i]})
      assert Qh == Q(y)
