import goodprime
import mult
import div

for n in range(0,8):
  N = 2^n

  for loop in range(10):
    goodprime.sqrtprime_internal.cache = {}
    goodprime.sqrtprime_product_internal.cache = {}
    goodprime.sqrtprime_suminv_internal.cache = {}
    goodprime.sqrtprime_scale_cache = {}

    while True:
      d = tuple(uniq(sorted([random_prime(100) for i in range(n)])))
      if len(d) == n:
        break

    print(n,d)
    for bits1 in range(0,64):
      for bits2 in range(0,64):
        f = tuple(ZZ.random_element(1-2^bits1,2^bits1) for j in range(N))
        g = tuple(ZZ.random_element(1-2^bits2,2^bits2) for j in range(N))
        h = mult.mult(n,N,d,f,g)
        if any(fj != 0 for fj in f):
          assert g == div.div(n,N,d,h,f)
        if any(gj != 0 for gj in g):
          assert f == div.div(n,N,d,h,g)
