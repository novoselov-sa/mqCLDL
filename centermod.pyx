from sage.rings.integer cimport Integer
from libc.stdlib cimport malloc, free
from sage.libs.gmp.mpz cimport *

def vector(f,q_py):
  cdef int N = len(f)
  cdef Integer q_Int = q_py
  cdef Integer v
  cdef int j
  cdef mpz_t q
  cdef mpz_t x
  cdef mpz_t y

  mpz_init_set(q,q_Int.value)
  mpz_init(x)
  mpz_init(y)

  result = []
  for j in range(N):
    v = Integer(f[j])
    mpz_set(x,v.value)
    mpz_fdiv_r(x,x,q)
    mpz_sub(y,x,q)
    if mpz_cmpabs(x,y) > 0:
      mpz_set(v.value,y)
    else:
      mpz_set(v.value,x)
    result.append(v)

  mpz_clear(x)
  mpz_clear(y)
  mpz_clear(q)
  return tuple(result)

def vector_mult(f,g,q_py):
  cdef int N = len(f)
  cdef Integer q_Int = q_py
  cdef Integer v
  cdef Integer w
  cdef int j
  cdef mpz_t q
  cdef mpz_t x
  cdef mpz_t y

  mpz_init_set(q,q_Int.value)
  mpz_init(x)
  mpz_init(y)

  result = []
  for j in range(N):
    v = Integer(f[j])
    w = Integer(g[j])
    mpz_set(x,v.value)
    mpz_set(y,w.value)
    mpz_mul(x,x,y)
    mpz_fdiv_r(x,x,q)
    mpz_sub(y,x,q)
    if mpz_cmpabs(x,y) > 0:
      mpz_set(v.value,y)
    else:
      mpz_set(v.value,x)
    result.append(v)

  mpz_clear(x)
  mpz_clear(y)
  mpz_clear(q)
  return tuple(result)
