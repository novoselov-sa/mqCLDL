# XXX: reduce redundancy with div and powerprod

import nprofile
profiling = ['bits','mainwork','subsetprod_internal','subsetprod']
import goodprime
import mult

def vector_scalar(n,s,f):
  return [s * fj for fj in f]

def hadamard(n,N,s,f):
  if n == 0: return f
  n -= 1
  N /= 2
  twistf = vector_scalar(n,s[n],f[N:])
  v0 = [f[i] + twistf[i] for i in range(N)]
  v1 = [f[i] - twistf[i] for i in range(N)]
  return hadamard(n,N,s,v0) + hadamard(n,N,s,v1)

def scaledinvhadamard(n,N,sinv,v):
  if n == 0: return v
  n -= 1
  N /= 2
  v0 = scaledinvhadamard(n,N,sinv,v[:N])
  v1 = scaledinvhadamard(n,N,sinv,v[N:])
  f0 = [v0[i] + v1[i] for i in range(N)]
  f1 = [sinv[n] * (v0[i] - v1[i]) for i in range(N)]
  return f0 + f1

@nprofile.profile
def bits(n,N,d,g,e):
  m = len(g)

  prec = 26
  while True:
    if len(d) == 0 or min(d) > 0:
      I = RealIntervalField(prec)
    else:
      I = ComplexIntervalField(prec)
    gI = [[I(gjk) for gjk in gj] for gj in g]
    s = [sqrt(I(di)) for di in d]
    sinv = [1/si for si in s]
    gtransform = [hadamard(n,N,s,gIj) for gIj in gI]
    htransform = [
      [I.prod(gtransform[j][k] for j in range(m) if ei[j]) for k in range(N)]
      for ei in e]
    h = [scaledinvhadamard(n,N,sinv,htransformi) for htransformi in htransform]
    h = [[hij/N for hij in hi] for hi in h]
    hbound = max(abs(hij).upper() for hi in h for hij in hi)
    if not hbound.is_infinity():
      B = ceil(hbound).nbits() + 1
      if B <= 16 or max(abs(hij).lower() for hi in h for hij in hi) >= 2^(B-5):
        return B,hbound

    prec *= 2

def transform_subsetprod_single(n,N,d,v,ei,low,high):
  m = len(v)
  result = []
  for t in range(low,high):
    q,s,sinv = goodprime.sqrtprime(d,t)
    K = Integers(q)
    result += [[lift(K.prod(K(v[j][t - low][k]) for j in range(m) if ei[j])) for k in range(N)]]
  return result

def transform_subsetprod(n,N,d,v,e,low,high):
  # XXX: use lupanov etc.
  return [transform_subsetprod_single(n,N,d,v,ei,low,high) for ei in e]

@nprofile.profile
def mainwork(n,N,d,g,e,low,high):
  gtransform = [mult.transform(n,N,d,gj,low,high) for gj in g]
  htransform = transform_subsetprod(n,N,d,gtransform,e,low,high)
  return [mult.invtransform(n,N,d,htransformi,low,high) for htransformi in htransform]

# assumes N = 2^n
# assumes d is length-n tuple of ZZ
# assumes not has_any_nontrivial_square_product(d)
# assumes g is length-m tuple; each g[j] is length-N tuple of ZZ
# assumes e is tuple, each e[i] being length-m tuple of 0 or 1
# returns tuple of:
#   g[0]^e[0][0]*g[1]^e[0][1]*...
#   g[0]^e[1][0]*g[1]^e[1][1]*...
#   g[0]^e[2][0]*g[1]^e[2][1]*...
#   etc.
@nprofile.profile
def subsetprod_internal(n,N,d,g,e):
  m = len(g)
  used = [False] * m
  for ei in e:
    if len(ei) != m: raise Exception('mismatched length')
    for j in range(m):
      if ei[j] == 1:
        used[j] = True
      elif ei[j] != 0:
        raise Exception('subsetprod works only with exponents 0 and 1')
  
  g = tuple(g[j] for j in range(m) if used[j])
  e = tuple(tuple(ei[j] for j in range(m) if used[j]) for ei in e)
  m = len(g)

  if m == 0: # optimization: no g used, so all outputs are 1
    return tuple((1,) + (0,)*(N-1) for ei in e)

  B,hbound = bits(n,N,d,g,e)
  low,high = goodprime.lowhigh(B)
  h = mainwork(n,N,d,g,e,low,high)
  assert(max(abs(hij) for hi in h for hij in hi) <= hbound)
  return h

# same semantics as subsetprod_internal
# extra optimization: checks for e[i] weight<2
@nprofile.profile
def subsetprod(n,N,d,g,e):
  esum = [sum([ZZ(eij) for eij in ei]) for ei in e]
  e2 = [e[i] for i in range(len(e)) if esum[i]>=2]
  h2 = subsetprod_internal(n,N,d,g,e2)
  k = 0
  h = ()
  for i in range(len(e)):
    if esum[i]>=2:
      h += (h2[k],)
      k += 1
    elif esum[i]==1:
      for j in range(len(g)):
        if e[i][j] == 1:
          h += (tuple(g[j]),)
          break
    elif esum[i]==0:
      h += ((1,) + (0,)*(N-1),)
  return h
