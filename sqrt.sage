import nprofile
profiling = ['twisted','squareroot']
import goodprime
import norm
import div
from memoized import memoized

@memoized
# assumes d is tuple of ZZ
# assumes not has_any_nontrivial_square_product(d)
# assumes k in range(len(d))
# returns q,s,sinv
# each r gives a new random choice, distinct from previous
def twistedprime(d,k,r):
  n = len(d)
  mu = tuple(-1 if i == k else 1 for i in range(n))
  while True:
    q = goodprime.goodprime(d,mu)
    if all(q != twistedprime(d,k,s) for s in range(r)):
      return q

# TwistedSquareRoot from paper
# assumes d is tuple of ZZ
# assumes not has_any_nontrivial_square_product(d)
# assumes h is ZZ
@nprofile.profile
def twisted(n,d,h):
  if h == 0: return 0,set()
  r = 0
  while True:
    q = [twistedprime(d,j,r) for j in range(n)]
    # XXX: use remainder trees
    hq4 = [h % (4*q[j]) for j in range(n)]
    hq = [h % qj for qj in q]
    if not 0 in hq: break
    r += 1
  S = set(j for j in range(n) if hq4[j].jacobi(q[j]) == -1)
  pi = prod(d[j] for j in range(n) if j in S)
  if h % pi != 0:
    raise Exception('%s does not divide %s' % (pi,h))
  h /= pi
  if not h.is_square():
    raise Exception('%s/%s is not square' % (h,pi))
  return sqrt(h),S

# SquareRoot from paper
# assumes d is tuple of ZZ with len(d) == n
# assumes not has_any_nontrivial_square_product(d)
# assumes m in {0,1,...,n}
# assumes h is tuple of ZZ with len(h) == 2^m
def squareroot_internal(n,d,m,h):
  if all(hj == 0 for hj in h):
    return h,set()
  if m == 0:
    g,S = twisted(n,d,h[0])
    return (g,),S
  hnorm = norm.normonestep(m,2^m,d[:m],h,m-1)
  m -= 1
  M = 2^m
  h0,h1 = h[:M],h[M:]
  s,S0 = squareroot_internal(m,d[:m],m,hnorm)
  if s == h0:
    t,S = squareroot_internal(n,d,m,h0)
  else:
    t,S = squareroot_internal(n,d,m,tuple(ZZ((h0[j]-s[j])/2) for j in range(M)))
  pi = prod(d[j] for j in range(n+1) if j in S and j != m)
  numer = tuple(ZZ(h1[j]/(2*pi)) for j in range(M))
  u = div.div(m,M,d[:m],numer,t)
  if m in S:
    S.discard(m)
    return u+t,S
  else:
    return t+u,S

@nprofile.profile
def squareroot(n,d,m,h):
  return squareroot_internal(n,d,m,h)
