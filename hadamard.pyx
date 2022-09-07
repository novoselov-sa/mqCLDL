from sage.rings.integer cimport Integer
from libc.stdlib cimport malloc, free
from sage.libs.gmp.mpz cimport *

def hadamard(n_py,N_py,q_py,s,f):
  cdef int n = n_py
  cdef int N = N_py
  cdef Integer q_Int = q_py
  cdef int i
  cdef int j
  cdef int k
  cdef int x
  cdef int twoi
  cdef Integer v
  cdef Integer r
  cdef mpz_t si
  cdef mpz_t sg1
  cdef mpz_t q
  cdef mpz_t *g = <mpz_t *> malloc(N * sizeof(mpz_t))
  if not g: raise MemoryError

  for j in range(N):
    v = Integer(f[j])
    mpz_init_set(g[j],v.value)

  mpz_init(si)
  mpz_init(sg1)
  mpz_init_set(q,q_Int.value)

  for i in range(n):
    v = s[i]
    mpz_set(si,v.value)
    twoi = 1 << i
    for j in range(0,N >> (i + 1)):
      x = j << (i + 1)
      for k in range(x,x + twoi):
        mpz_mul(sg1,si,g[k + twoi])
        mpz_mod(sg1,sg1,q)
        mpz_sub(g[k + twoi],g[k],sg1)
        mpz_add(g[k],g[k],sg1)

  mpz_clear(si)
  mpz_clear(sg1)
  mpz_clear(q)

  result = []
  for j in range(N):
    r = Integer(0)
    mpz_set(r.value,g[j])
    mpz_clear(g[j])
    result.append(r)

  free(g)
  return tuple(result)

def scaledinvhadamard(n_py,N_py,q_py,s,f):
  cdef int n = n_py
  cdef int N = N_py
  cdef Integer q_Int = q_py
  cdef int i
  cdef int j
  cdef int k
  cdef int x
  cdef int twoi
  cdef Integer v
  cdef Integer r
  cdef mpz_t si
  cdef mpz_t sg1
  cdef mpz_t q
  cdef mpz_t *g = <mpz_t *> malloc(N * sizeof(mpz_t))
  if not g: raise MemoryError

  for j in range(N):
    v = Integer(f[j])
    mpz_init_set(g[j],v.value)

  mpz_init(si)
  mpz_init(sg1)
  mpz_init_set(q,q_Int.value)

  for i in range(n):
    v = s[i]
    mpz_set(si,v.value)
    twoi = 1 << i
    for j in range(0,N >> (i + 1)):
      x = j << (i + 1)
      for k in range(x,x + twoi):
        mpz_sub(sg1,g[k],g[k + twoi])
        mpz_mul(sg1,sg1,si)
        mpz_add(g[k],g[k],g[k + twoi])
        mpz_mod(g[k + twoi],sg1,q)

  for i in range(n):
    mpz_mod(g[i],g[i],q)

  mpz_clear(si)
  mpz_clear(sg1)
  mpz_clear(q)

  result = []
  for j in range(N):
    r = Integer(0)
    mpz_set(r.value,g[j])
    mpz_clear(g[j])
    result.append(r)

  free(g)
  return tuple(result)

def invhadamard(n_py,N_py,q_py,s,f):
  cdef int n = n_py
  cdef int N = N_py
  cdef Integer q_Int = q_py
  cdef int i
  cdef int j
  cdef int k
  cdef int x
  cdef int twoi
  cdef Integer v
  cdef Integer r
  cdef mpz_t si
  cdef mpz_t sg1
  cdef mpz_t q
  cdef mpz_t *g = <mpz_t *> malloc(N * sizeof(mpz_t))
  if not g: raise MemoryError

  for j in range(N):
    v = Integer(f[j])
    mpz_init_set(g[j],v.value)

  mpz_init(si)
  mpz_init(sg1)
  mpz_init_set(q,q_Int.value)

  for i in range(n):
    v = s[i]
    mpz_set(si,v.value)
    twoi = 1 << i
    for j in range(0,N >> (i + 1)):
      x = j << (i + 1)
      for k in range(x,x + twoi):
        mpz_sub(sg1,g[k],g[k + twoi])
        mpz_mul(sg1,sg1,si)
        mpz_add(g[k],g[k],g[k + twoi])
        mpz_mod(g[k + twoi],sg1,q)

  for i in range(n):
    mpz_mod(g[i],g[i],q)

  result = []
  for j in range(N):
    for i in range(n):
      if mpz_tstbit(g[j],0):
        mpz_add(g[j],g[j],q)
      mpz_fdiv_q_2exp(g[j],g[j],1)
    r = Integer(0)
    mpz_set(r.value,g[j])
    mpz_clear(g[j])
    result.append(r)

  mpz_clear(si)
  mpz_clear(sg1)
  mpz_clear(q)
  free(g)
  return tuple(result)
