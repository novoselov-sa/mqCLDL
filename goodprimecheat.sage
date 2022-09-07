merge = 8
wheel = 30030
units = tuple(j for j in range(wheel) if gcd(j,wheel)==1)

def setlambda(b):
  global bits
  global A
  global B
  bits = b
  A = ZZ(floor(2^(bits-1)/wheel))
  B = ZZ(ceil(2^bits/wheel))

setlambda(64)

proof.arithmetic(False)

from memoized import memoized
import nprofile
profiling = ['goodprime','sqrtprime','sqrtprime_product','sqrtprime_suminv','sqrtprime_scale_cachefill']

# 1/s[0],1/s[1],... assuming all invertible
def recip(s):
  if 0: # montgomery approach
    n = len(s)
    if n == 0: return s
    result = [s[0]]
    for j in range(1,n):
      result += [result[-1]*s[j]]
    recip = 1/result[-1]
    for j in range(n-1,0,-1):
      result[j],recip = recip*result[j-1],recip*s[j]
    result[0] = recip
    return tuple(result)
  return tuple(1/sj for sj in s)

@memoized
def fastjacobi(d,q):
  return d.jacobi(q)

def jacobi(d,q):
  return fastjacobi(d,q % (4 * d))

# assumes not has_any_nontrivial_square_product(d)
# assumes mu is a tuple of {-1,1} with len(mu) == len(d)
@nprofile.profile
def goodprime(d,mu):
  while True:
    q = ZZ.random_element(A,B)*wheel + units[randrange(len(units))]
    if all(q % dj for dj in d if abs(dj) >= 2):
      if all(jacobi(dj,q) == muj for dj,muj in zip(d,mu)):
          if q.is_prime():
            if q >= 2^(bits-1):
              if q < 2^bits:
                return q

@memoized
# assumes d is tuple of ZZ
# assumes not has_any_nontrivial_square_product(d)
# returns q,s,sinv
# each j>=0 gives a new random choice, distinct from previous
# negative j's merge ranges of nonnegative j's
def sqrtprime_internal(d,j):
  if j < 0:
    n = len(d)
    # skip 0..merge-1
    high = merge*(1-j)
    low = high - merge
    x = [sqrtprime_internal(d,i) for i in range(low,high)]
    q = sqrtprime_product(d,low,high)
    s = tuple(sum(((x[i-low][1][j]*sqrtprime_scale(d,low,high,i)) % x[i-low][0]) * (q / x[i-low][0]) for i in range(low,high)) % q for j in range(n))
    sinv = tuple(sum(((x[i-low][2][j]*sqrtprime_scale(d,low,high,i)) % x[i-low][0]) * (q / x[i-low][0]) for i in range(low,high)) % q for j in range(n))
    return q,s,sinv
  while True:
    q = goodprime(d,(1,)*len(d))
    s = tuple(lift(sqrt(Mod(dj,q))) for dj in d)
    sinv = recip([Mod(sj,q) for sj in s])
    sinv = tuple(lift(sj) for sj in sinv)
    if all((q,s,sinv) != sqrtprime_internal(d,i) for i in range(j)):
      return q,s,sinv

def sqrtprime(d,j):
  return sqrtprime_internal(d,j)

def lowhigh(B):
  numprimes = ceil(B / (bits - 1))
  notmerged = numprimes % merge
  if notmerged < 3: # XXX: measure to find optimum
    merged = numprimes // merge
    return -merged,notmerged
  return -ceil(numprimes / merge),0

def nextlowhigh(B,low,high):
  numprimes = ceil(B / (bits - 1))
  if low < 0:
    return low - numprimes,low
  return high,high + numprimes

@memoized
def sqrtprime_product_internal(d,low,high):
  if high == low + 1:
    q,s,sinv = sqrtprime(d,low)
    return q
  mid = low + 1
  while mid + (mid - low) < high:
    mid = mid + (mid - low)
  return sqrtprime_product_internal(d,low,mid) * sqrtprime_product_internal(d,mid,high)

# assumes d is tuple of ZZ
# assumes not has_any_nontrivial_square_product(d)
# returns product of q for j in range(low,high)
def sqrtprime_product(d,low,high):
  low = ZZ(low)
  high = ZZ(high)
  if high <= low:
    raise Exception('%s must be larger than %s',high,low)
  return sqrtprime_product_internal(d,low,high)

@memoized
def sqrtprime_suminv_internal(d,low,high):
  if high == low + 1:
    return 1
  mid = low + 1
  while mid + (mid - low) < high:
    mid = mid + (mid - low)
  return (sqrtprime_suminv_internal(d,low,mid) * sqrtprime_product_internal(d,mid,high)
        + sqrtprime_suminv_internal(d,mid,high) * sqrtprime_product_internal(d,low,mid))

@nprofile.profile
# assumes d is tuple of ZZ
# assumes not has_any_nontrivial_square_product(d)
# returns numerator of sum 1/qj for j in range(low,high)
def sqrtprime_suminv(d,low,high):
  low = ZZ(low)
  high = ZZ(high)
  if high <= low:
    raise Exception('%s must be larger than %s',high,low)
  return sqrtprime_suminv_internal(d,low,high)

sqrtprime_scale_cache = {}

# X is congruent to suminv(low,high) mod product(low2,high2)
def sqrtprime_scale_cachefill_internal(d,low,high,X,low2,high2):
  q = sqrtprime_product_internal(d,low2,high2)
  if high2 == low2 + 1:
    sqrtprime_scale_cache[d,low,high,low2] = ZZ(1/Mod(X,q))
    return
  X %= q
  mid = low2 + 1
  while mid + (mid - low2) < high2:
    mid = mid + (mid - low2)
  sqrtprime_scale_cachefill_internal(d,low,high,X,low2,mid)
  sqrtprime_scale_cachefill_internal(d,low,high,X,mid,high2)

@nprofile.profile
def sqrtprime_scale_cachefill(d,low,high):
  X = sqrtprime_suminv(d,low,high)
  sqrtprime_scale_cachefill_internal(d,low,high,X,low,high)

# returns reciprocal of product(low,high) mod qj
def sqrtprime_scale(d,low,high,j):
  if not (d,low,high,j) in sqrtprime_scale_cache:
    sqrtprime_scale_cachefill(d,low,high)
  return sqrtprime_scale_cache[d,low,high,j]
