import mult
import subsetprod

for n in range(5):
  print(n)
  d = (-2,3,-5,7,-11,13,-17,19)[:n]
  N = 2^n

  for bits in [1,2,4,8,16,32,64,128,256]:
    for inputs in range(8):
      for outputs in range(8):
        g = [[ZZ.random_element(1-2^bits,2^bits) for k in range(N)] for j in range(inputs)]
    e = [[randrange(2) for j in range(inputs)] for i in range(outputs)]
    h = subsetprod.subsetprod(n,N,d,g,e)
    assert len(h) == outputs
    for hi in h: assert len(hi) == N
    for hi,ei in zip(h,e):
      check = (1,) + (0,)*(N-1)
      for j in range(inputs):
        if ei[j]:
          check = mult.mult(n,N,d,check,g[j])
      assert check == hi
