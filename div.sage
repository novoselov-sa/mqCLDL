import nprofile
profiling = ['transform_invertible','transform_divide','div_bits','div_mainwork','div']
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
def transform_invertible(n,N,d,w,low,high):
  for j in range(low,high):
    q,s,sinv = goodprime.sqrtprime(d,j)
    wj = w[j - low]
    for i in range(N):
      if wj[i] % q == 0:
        return False
  return True

@nprofile.profile
def transform_divide(n,N,d,v,w,low,high):
  result = []
  for j in range(low,high):
    q,s,sinv = goodprime.sqrtprime(d,j)
    vj = v[j - low]
    wj = w[j - low]
    mont = [wj[0]]
    for i in range(1,N):
      mont.append(mont[-1] * wj[i] % q)
    recip = lift(1/Mod(mont[-1],q))
    for i in range(N-1,0,-1):
      mont[i] = recip*mont[i - 1] % q
      recip = recip*wj[i] % q
    mont[0] = recip
    result += [[vj[i] * mont[i] % q for i in range(N)]]
  return result

@nprofile.profile
# assumes g nonzero (otherwise loops forever!)
def div_bits(n,N,d,h,g):
  prec = 26
  while True:
    if len(d) == 0 or min(d) > 0:
      I = RealIntervalField(prec)
    else:
      I = ComplexIntervalField(prec)
    hI = [I(hj) for hj in h]
    gI = [I(gj) for gj in g]
    s = [sqrt(I(di)) for di in d]
    sinv = [1/si for si in s]
    gtransform = hadamard(n,N,s,gI)
    if all(not 0 in v for v in gtransform):
      htransform = hadamard(n,N,s,hI)
      ftransform = [htransform[j]/gtransform[j] for j in range(N)]
      fI = scaledinvhadamard(n,N,sinv,ftransform)
      fI = [fj/N for fj in fI]
      fbound = max(abs(fj).upper() for fj in fI)
      if not fbound.is_infinity():
        B = ceil(fbound).nbits() + 1 # so fbound <= 2^(B-1)-1
        if max(abs(fj).lower() for fj in fI) >= 2^(B-5):
          return B,fbound
    prec *= 2
  
@nprofile.profile
def div_mainwork(n,N,d,h,g,low,high):
  gtransform = mult.transform(n,N,d,g,low,high)
  assert transform_invertible(n,N,d,gtransform,low,high)
  htransform = mult.transform(n,N,d,h,low,high)
  ftransform = transform_divide(n,N,d,htransform,gtransform,low,high)
  return mult.invtransform(n,N,d,ftransform,low,high)

# assumes N = 2^n
# assumes d is tuple of ZZ with len(d) == n
# assumes not has_any_nontrivial_square_product(d)
# assumes h is tuple of ZZ with len(h) == N
# assumes g is tuple of ZZ with len(g) == N
@nprofile.profile
def div(n,N,d,h,g):
  if all(gj == 0 for gj in g): raise Exception('division by zero')
  if all(hj == 0 for hj in h): return h
  B,fbound = div_bits(n,N,d,h,g)
  low,high = goodprime.lowhigh(B)
  while True:
    try:
      f = div_mainwork(n,N,d,h,g,low,high)
      break
    except:
      low,high = goodprime.nextlowhigh(B,low,high)
  assert(max(abs(fj) for fj in f) <= fbound)
  return f
