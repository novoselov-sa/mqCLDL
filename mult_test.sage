import goodprime
import mult
import div

n = 2
N = 4
d = (67,97)
goodprime.sqrtprime_internal.cache[d,0] = (9350342201084777911, (4662121590309669551, 3975690820250549520), (8443024681423826033, 9198538143539221145))
f = (-25672053, 44811798, 49776505, -59708499)
g = (-2580735, 2843271, 8780279, 10322305) 
h = mult.mult(n,N,d,f,g)

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

    print n,d
    for bits1 in range(0,64):
      for bits2 in range(0,64):
        f = tuple(ZZ.random_element(1-2^bits1,2^bits1) for j in range(N))
        g = tuple(ZZ.random_element(1-2^bits2,2^bits2) for j in range(N))
        h = mult.mult(n,N,d,f,g)
