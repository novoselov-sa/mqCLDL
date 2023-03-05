import nprofile
#profiling = ['hadamard','invhadamard','transform','transform_multiply','invtransform','mult_mainwork','mult','square']
import goodprime
import hadamard as fasthadamard
import centermod

def hadamard(n,N,q,s,f):
  return fasthadamard.hadamard(n,N,q,s,f)

# hadamard(n,N,q,s,f) where q,s,sinv = sqrtprime(d,j) for j in range(low,high)
def transform_internal(n,N,d,f,low,high):
  if high == low + 1:
    q,s,sinv = goodprime.sqrtprime(d,low)
    result = hadamard(n,N,q,s,f)
    return (result,)
  mid = low + 1
  while mid + (mid - low) < high:
    mid = mid + (mid - low)
  q0 = goodprime.sqrtprime_product(d,low,mid)
  q1 = goodprime.sqrtprime_product(d,mid,high)
  f0 = centermod.vector(f,q0)
  f1 = centermod.vector(f,q1)
  return transform_internal(n,N,d,f0,low,mid) + transform_internal(n,N,d,f1,mid,high)

#@nprofile.profile
def transform(n,N,d,f,low,high):
  return transform_internal(n,N,d,f,low,high)

def invhadamard(n,N,q,sinv,v):
  return fasthadamard.invhadamard(n,N,q,sinv,v)

def invtransform_internal(n,N,d,v,low,high,outerlow,outerhigh):
  if high == low + 1:
    q,s,sinv = goodprime.sqrtprime(d,low)
    V = invhadamard(n,N,q,sinv,list(v[0]))
    x = goodprime.sqrtprime_scale(d,outerlow,outerhigh,low)
    if x == 1: return V
    return centermod.vector([x*v for v in V],q)
  mid = low + 1
  while mid + (mid - low) < high:
    mid = mid + (mid - low)
  q0 = goodprime.sqrtprime_product(d,low,mid)
  q1 = goodprime.sqrtprime_product(d,mid,high)
  f0 = invtransform_internal(n,N,d,v[:mid - low],low,mid,outerlow,outerhigh)
  f1 = invtransform_internal(n,N,d,v[mid - low:],mid,high,outerlow,outerhigh)
  return [f0j * q1 + f1j * q0 for f0j,f1j in zip(f0,f1)]

#@nprofile.profile
def invtransform(n,N,d,v,low,high):
  h = invtransform_internal(n,N,d,v,low,high,low,high)
  q = goodprime.sqrtprime_product(d,low,high)
  return centermod.vector(h,q)

#@nprofile.profile
def transform_multiply(n,N,d,v,w,low,high):
  result = []
  for j in range(low,high):
    vj = v[j - low]
    wj = w[j - low]
    q,s,sinv = goodprime.sqrtprime(d,j)
    result += [centermod.vector_mult(vj,wj,q)]
  return result

#@nprofile.profile
def mult_mainwork(n,N,d,f,g,low,high):
  ftransform = transform(n,N,d,f,low,high)
  gtransform = transform(n,N,d,g,low,high)
  htransform = transform_multiply(n,N,d,ftransform,gtransform,low,high)
  return invtransform(n,N,d,htransform,low,high)

# assumes N = 2^n
# assumes d is tuple of ZZ with len(d) == n
# assumes not has_any_nontrivial_square_product(d)
# assumes f is tuple of ZZ with len(f) == N
# assumes g is tuple of ZZ with len(g) == N
#@nprofile.profile
def mult(n,N,d,f,g):
  fmax = max(abs(fj) for fj in f)
  gmax = max(abs(gj) for gj in g)
  hbound = fmax * gmax * prod(1 + abs(dj) for dj in d)
  B = hbound.nbits() + 1 # so hbound <= 2^(B-1)-1
  low,high = goodprime.lowhigh(B)
  h = mult_mainwork(n,N,d,f,g,low,high)
  assert(all(abs(hj) <= hbound for hj in h))
  return h

#@nprofile.profile
def square(n,N,d,f):
  fmax = max(abs(fj) for fj in f)
  hbound = fmax^2 * prod(1 + abs(dj) for dj in d)
  B = hbound.nbits() + 1 # so hbound <= 2^(B-1)-1
  low,high = goodprime.lowhigh(B)
  ftransform = transform(n,N,d,f,low,high)
  htransform = transform_multiply(n,N,d,ftransform,ftransform,low,high)
  h = invtransform(n,N,d,htransform,low,high)
  assert(all(abs(hj) <= hbound for hj in h))
  return h
