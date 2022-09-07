from field import *

for n in range(0,7):
  print n
  N = 2^n
  d = (-2,3,-5,7,-11,13,-17,19,-23,29)[:n]
  K = field(d)

  P = PolynomialRing(QQ,'x',n)
  x = P.gens()
  Q = P.quotient([x[i]^2 - d[i] for i in range(n)])

  for bits in range(20):
    f = K(tuple(ZZ.random_element(1-2^bits,2^bits) for j in range(N)),ZZ.random_element(1,1+2^bits))
    h = f.square()
    assert h == f * f
    s = h.sqrt()
    assert s == f or s == -f
    g = K(tuple(ZZ.random_element(1-2^bits,2^bits) for j in range(N)),ZZ.random_element(1,1+2^bits))
    assert f - g == f + (-g)
    assert (f + g) - g == f
    q = K(tuple(ZZ.random_element(1-2^bits,2^bits) for j in range(N)),1)
    h = f * q
    if f.is_nonzero():
      assert h.divexact(f) == q
      assert (q,) == powerprod(d,(f,h),((-1,1),))

    Pf = P(sum(f.numer.c[j] * prod(x[i] for i in range(n) if j & (1 << i)) for j in range(N)))/f.denom
    Qf = Q(Pf)
    Qg = Q(sum(g.numer.c[j] * prod(x[i] for i in range(n) if j & (1 << i)) for j in range(N)))/g.denom

    h = f + g
    Qh = Q(sum(h.numer.c[j] * prod(x[i] for i in range(n) if j & (1 << i)) for j in range(N)))/h.denom
    assert Qh == Qf + Qg

    h = f - g
    Qh = Q(sum(h.numer.c[j] * prod(x[i] for i in range(n) if j & (1 << i)) for j in range(N)))/h.denom
    assert Qh == Qf - Qg

    h = -f
    Qh = Q(sum(h.numer.c[j] * prod(x[i] for i in range(n) if j & (1 << i)) for j in range(N)))/h.denom
    assert Qh == -Qf

    h = f * g
    Qh = Q(sum(h.numer.c[j] * prod(x[i] for i in range(n) if j & (1 << i)) for j in range(N)))/h.denom
    assert Qh == Qf * Qg
    assert (h,) == subsetprod(d,(f,g),((1,1),))
    assert (h,h) == subsetprod(d,(f,f,g,g,h,h),((1,0,0,1,0,0),(0,1,1,0,0,0)))
    for t in range(10):
      assert h.symbol(t) == f.symbol(t) * g.symbol(t)
    assert h.symbols(2,20) == tuple(h.symbol(t) for t in range(2,20))

    h = f.square()
    Qh = Q(sum(h.numer.c[j] * prod(x[i] for i in range(n) if j & (1 << i)) for j in range(N)))/h.denom
    assert Qh == Qf^2

    for t in range(n):
      h = f.normonestep(t)
      z = [0 if (j&(1<<t)) else h.numer.c[(j % (1<<t)) + ((j >> (t+1)) << t)] for j in range(N)]
      Qh = Q(sum(z[j] * prod(x[i] for i in range(n) if j & (1 << i)) for j in range(N)))/h.denom
      assert Qh == Qf*Pf.substitute({x[t]:-x[t]})
