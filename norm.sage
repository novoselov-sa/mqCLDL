import nprofile
profiling = ['normonestep','absnorm']
import goodprime
import mult

def transform_normonestep(n,N,d,v0,v1,dt,low,high):
  result = []
  for j in range(low,high):
    v0j = v0[j - low]
    v1j = v1[j - low]
    q,s,sinv = goodprime.sqrtprime(d,j)
    result += [[(v0j[i]^2 - dt * v1j[i]^2) % q for i in range(N)]]
  return result

@nprofile.profile
def normonestep(n,N,d,f,t):
  if not t in range(n):
    raise Exception('normonestep %s not in range(%s)' % (t,n))
  fmax = max(abs(fj) for fj in f)
  hbound = fmax^2 * prod(1 + abs(dj) for dj in d)
  B = hbound.nbits() + 1 # so hbound <= 2^(B-1)-1
  low,high = goodprime.lowhigh(B)
  f0 = [f[j] for j in range(N) if not (j & (1 << t))]
  f1 = [f[j] for j in range(N) if j & (1 << t)]
  dt = d[t]
  n -= 1
  N /= 2
  d = d[:t] + d[t+1:]
  v0 = mult.transform(n,N,d,f0,low,high)
  v1 = mult.transform(n,N,d,f1,low,high)
  v = transform_normonestep(n,N,d,v0,v1,dt,low,high)
  h = mult.invtransform(n,N,d,v,low,high)
  assert(all(abs(hj) <= hbound for hj in h))
  return h

@nprofile.profile
def absnorm(n,N,d,f):
  while n > 0:
    f = normonestep(n,N,d,f,n-1)
    n -= 1
    N /= 2
    d = d[:n]
  return f[0]
