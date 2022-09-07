# XXX: reduce redundancy with div and subsetprod

import nprofile
profiling = ['bits','transform_invertible','transform_recip','mainwork','scaledpowerprod']
import goodprime
import centermod
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
def bits(n,N,d,g,denom,e,multiplier):
  m = len(g)
  inverting = [any(ei[j] < 0 for ei in e) for j in range(m)]

  for j in range(m):
    if inverting[j]:
      if all(gjk == 0 for gjk in g[j]):
        raise Exception('division by zero')

  prec = 26
  while True:
    if len(d) == 0 or min(d) > 0:
      I = RealIntervalField(prec)
    else:
      I = ComplexIntervalField(prec)
    gI = [[I(gjk)/denom[j] for gjk in g[j]] for j in range(m)]
    s = [sqrt(I(di)) for di in d]
    sinv = [1/si for si in s]
    gtransform = [None] * m
    ok = True
    for j in range(m):
      if inverting[j]:
        gtransform[j] = hadamard(n,N,s,gI[j])
        if any(0 in v for v in gtransform[j]):
          ok = False
          break
    if ok:
      for j in range(m):
        if not inverting[j]:
          gtransform[j] = hadamard(n,N,s,gI[j])
      htransform = [
        [I.prod(gtransform[j][k]^ei[j] for j in range(m) if ei[j]) for k in range(N)]
        for ei in e]
      h = [scaledinvhadamard(n,N,sinv,htransformi) for htransformi in htransform]
      # skip division by N since we want scaled result
      hbound = abs(multiplier) * max(abs(hij).upper() for hi in h for hij in hi)
      if not hbound.is_infinity() and not hbound.is_NaN():
        B = ceil(hbound).nbits() + 1
        if B <= 16+n or max(abs(multiplier) * abs(hij).lower() for hi in h for hij in hi) >= 2^(B-5):
          return B,hbound

    prec *= 2

@nprofile.profile
def transform_invertible(n,N,d,w,low,high):
  for j in range(low,high):
    q,s,sinv = goodprime.sqrtprime(d,j)
    wj = w[j - low]
    for i in range(N):
      if wj[i] % q == 0:
        return False
  return True

@nprofile.profile
def transform_recip(n,N,d,w,low,high):
  result = []
  for j in range(low,high):
    q,s,sinv = goodprime.sqrtprime(d,j)
    wj = w[j - low]
    mont = [wj[0]]
    for i in range(1,N):
      mont.append(mont[-1] * wj[i] % q)
    recip = lift(1/Mod(mont[-1],q))
    for i in range(N-1,0,-1):
      mont[i] = recip*mont[i - 1] % q
      recip = recip*wj[i] % q
    mont[0] = recip
    result += [mont]
  return result

def transform_powerprod_single(n,N,d,v,vinv,denom,ei,low,high):
  m = len(v)
  result = []
  for t in range(low,high):
    q,s,sinv = goodprime.sqrtprime(d,t)
    K = Integers(q)
    result += [
      [lift(K.prod(K(v[j][t - low][k]/denom[j])^ei[j] if ei[j]>0 else 
                   K(vinv[j][t - low][k]*denom[j])^(-ei[j])
                   for j in range(m)
                   if ei[j])
           ) for k in range(N)
      ]
    ]
  return result

def transform_powerprod(n,N,d,v,vinv,denom,e,low,high):
  m = len(v)
  # XXX: use pippenger etc.
  return [transform_powerprod_single(n,N,d,v,vinv,denom,ei,low,high) for ei in e]

@nprofile.profile
def mainwork(n,N,d,g,denom,e,low,high):
  m = len(g)
  inverting = [any(ei[j] < 0 for ei in e) for j in range(m)]

  gtransform = [None] * m
  gtransforminv = [None] * m
  for j in range(m):
    if inverting[j]:
      gtransform[j] = mult.transform(n,N,d,g[j],low,high)
      assert transform_invertible(n,N,d,gtransform[j],low,high)
  for j in range(m):
    if inverting[j]:
      gtransforminv[j] = transform_recip(n,N,d,gtransform[j],low,high)
    else:
      gtransform[j] = mult.transform(n,N,d,g[j],low,high)

  htransform = transform_powerprod(n,N,d,gtransform,gtransforminv,denom,e,low,high)
  return [mult.invtransform(n,N,d,htransformi,low,high) for htransformi in htransform]

# assumes N = 2^n
# assumes d is length-n tuple of ZZ
# assumes not has_any_nontrivial_square_product(d)
# assumes g is length-m tuple; each g[j] is length-N tuple of ZZ
# assumes denom is length-m tuple of ZZ, all nonzero
# assumes e is tuple; each e[i] is length-m tuple of ZZ
# returns tuple of:
#   multiplier*N*(g[0]/denom[0])^e[0][0]*(g[1]/denom[1])^e[0][1]*...
#   multiplier*N*(g[0]/denom[0])^e[1][0]*(g[1]/denom[1])^e[1][1]*...
#   multiplier*N*(g[0]/denom[0])^e[2][0]*(g[1]/denom[1])^e[2][1]*...
#   etc.
# note the multiplier*N factor!
# assumes that each of these results has integer coefficients
@nprofile.profile
def scaledpowerprod(n,N,d,g,denom,e,multiplier):
  m = len(g)
  used = [False] * m
  for ei in e:
    if len(ei) != m: raise Exception('mismatched length')
    for j in range(m):
      if ei[j]:
        used[j] = True

  g = tuple(g[j] for j in range(m) if used[j])
  denom = tuple(denom[j] for j in range(m) if used[j])
  e = tuple(tuple(ZZ(ei[j]) for j in range(m) if used[j]) for ei in e)
  m = len(g)

  if m == 0: # optimization: no g used, so all outputs are multiplier*N
    return tuple((multiplier*N,) + (0,)*(N-1) for ei in e)

  B,hbound = bits(n,N,d,g,denom,e,multiplier)
  low,high = goodprime.lowhigh(B)
  while True:
    try:
      q = goodprime.sqrtprime_product(d,low,high)
      h = mainwork(n,N,d,g,denom,e,low,high)
      h = tuple(centermod.vector(tuple(multiplier*N*hj for hj in hi),q) for hi in h)
      break
    except:
      low,high = goodprime.nextlowhigh(B,low,high)
  assert(max(abs(hij) for hi in h for hij in hi) <= hbound)
  return h
