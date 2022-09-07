import mult
import powerprod

for n in range(5):
  d = (-2,3,-5,7,-11,13,-17,19)[:n]
  N = 2^n

  for bits in [1,2,4,8,16,32]:
    print n,bits
    for inputs in range(8):
      for outputs in range(8):
	denom = [ZZ.random_element(10,20) for j in range(inputs)]
        g = [[denom[j]*ZZ.random_element(1-2^bits,2^bits) for k in range(N)] for j in range(inputs)]
	e = [[randrange(5) for j in range(inputs)] for i in range(outputs)]
	h = powerprod.scaledpowerprod(n,N,d,g,denom,e,1)
	assert len(h) == outputs
	for hi in h: assert len(hi) == N
	for hi,ei in zip(h,e):
	  check = (N,) + (0,)*(N-1)
	  for j in range(inputs):
	    for k in range(ei[j]):
	      check = mult.mult(n,N,d,check,g[j])
	      check = tuple(ZZ(v/denom[j]) for v in check)
	  assert check == hi

	denom = [ZZ.random_element(10,20) for j in range(inputs)]
        g = [[denom[j]*ZZ.random_element(1,2^bits) for k in range(N)] for j in range(inputs)]
	e1 = [[randrange(0,5) for j in range(inputs)] for i in range(outputs)]
	e2 = [[randrange(0,5) for j in range(inputs)] for i in range(outputs)]
	h1 = powerprod.scaledpowerprod(n,N,d,g,denom,e1,2^100)
	h2 = powerprod.scaledpowerprod(n,N,d,g,denom,e2,2^100)
	e = [[e2[i][j] - e1[i][j] for j in range(inputs)]+[ZZ(i==j) for j in range(outputs)] for i in range(outputs)]
	h3 = powerprod.scaledpowerprod(n,N,d,g + list(h1),denom + [1]*outputs,e,2^300)
	h4 = powerprod.scaledpowerprod(n,N,d,h2,(1,)*outputs,[[i == j for j in range(outputs)] for i in range(outputs)],2^300)
	assert h3 == h4
