import char
import div
import ring
import nprofile
profiling = ['parameters','keygen_try','keygen','randommessage','encrypt','expandseckey','decrypt']

@nprofile.profile
def parameters(d,scale=0.5):
  n = len(d)
  N = 2^n
  D = [prod(d[i] for i in range(n) if j & (1 << i)) for j in range(N)]

  badprimeslimit = ZZ(ceil(0.5*scale*N/log(1.0+N)))
  badprimes = prod(primes(2,badprimeslimit))
  gmax = tuple(ZZ(ceil(2.0^(scale*N/log(1.0+N))*sqrt(abs(D[N-1]/D[j])))) for j in range(N))
  return badprimeslimit,badprimes,gmax

# k is finite field
# R = k[x[0],...,x[n-1]]/(x[0]^2-d[0],...,x[n-1]^2-d[n-1])
# returns uniform random element of R
def random_element(d,k):
  n = len(d)
  N = 2^n
  return [k.random_element() for j in range(N)]

# k is finite field; XXX: see degree restriction below
# R = k[x[0],...,x[n-1]]/(x[0]^2-d[0],...,x[n-1]^2-d[n-1])
# returns uniform random invertible element of R
def random_invertible_element(d,k):
  n = len(d)
  N = 2^n
  if n == 0:
    while True:
      r = k.random_element()
      if r != 0: return [r]
  if not k(d[-1]).is_square():
    P.<x> = k[]
    K = k.extension(x^2-d[-1],'s')
    z = random_invertible_element(d[:-1],K)
    # XXX: using polynomial() here is ok only for deg K <= 2
    return [zj.polynomial()[0] for zj in z] + [zj.polynomial()[1] for zj in z]
  if k(2*d[-1]) == 0:
    # x[n-1]^2-d[n-1] = (x[n-1]-d[n-1])^2
    # R = R'[x[n-1]]/(x[n-1]-d[n-1])^2
    # z + (x[n-1]-d[n-1])y is invertible iff z is invertible
    z = random_invertible_element(d[:-1],k)
    y = random_element(d[:-1],k)
    return [z[j]-d[-1]*y[j] for j in range(N/2)] + y
  s = sqrt(k(d[-1]))
  # s is nonzero and 2 is invertible
  # R = R'[x[n-1]]/(x[n-1]-s)(x[n-1]+s)
  z = random_invertible_element(d[:-1],k)
  y = random_invertible_element(d[:-1],k)
  return [(z[j]+y[j])/2 for j in range(N/2)] + [(z[j]-y[j])/(2*s) for j in range(N/2)]

def keygen_g(d,params):
  badprimeslimit,badprimes,gmax = params
  n = len(d)
  N = 2^n
  R = ring.ring(d)

  u = [Mod(0,1) for j in range(N)]
  for p in primes(2,badprimeslimit):
    v = random_invertible_element(d,GF(p))
    u = [u[j].crt(v[j]) for j in range(N)]
  u = [lift(u[j] + badprimes//2) - badprimes//2 for j in range(N)]

  s = [ZZ.random_element(2*gmax[j]+1)-gmax[j] for j in range(N)]
  g = [u[j] + badprimes*s[j] for j in range(N)]
  return R(g)

@nprofile.profile
def keygen_try(d,params):
  n = len(d)

  g = keygen_g(d,params)
  q = abs(g.absnorm())

  if q % 2 == 0: raise Exception('norm is even')
  for dj in d:
    if q % dj == 0:
      raise Exception('norm is divisible by %s' % dj)

  qn = g.quadnorms()
  if n > 0:
    mont = [Mod(qn[i].c[1],q) for i in range(n)]
    for i in range(1,n):
      mont[i] *= mont[i-1]
    recip = mont[-1]
    if not recip.is_unit(): raise Exception('non-invertible')
    recip = 1/recip
    for i in range(n-1,0,-1):
      mont[i] = mont[i-1]*recip
      recip *= qn[i].c[1]
    mont[0] = recip

  qs = tuple(lift(-qn[i].c[0]*mont[i]) for i in range(n))

  pubkey = q,qs
  seckey = g,q
  return pubkey,seckey

@nprofile.profile
def keygen(d,params):
  while True:
    try:
      return keygen_try(d,params)
    except:
      pass

@nprofile.profile
def randommessage(d,params,gap=32.0):
  badprimeslimit,badprimes,gmax = params
  n = len(d)
  N = 2^n
  mmax = tuple(ZZ(floor(badprimes*gmax[j]/gap)) for j in range(N))
  return tuple(ZZ(randrange(2*mmax[i]+1)-mmax[i]) for i in range(len(mmax)))

@nprofile.profile
def encrypt(d,message,pubkey):
  q,qs = pubkey
  n = len(d)
  N = 2^n
  return char.evaluate(n,N,q,qs,message)

@nprofile.profile
def expandseckey(d,seckey):
  g,q = seckey
  n = len(d)
  N = 2^n
  R = ring.ring(d)
  # XXX: speed up division using fact that numerator is integer
  qoverg = R((q,) + (0,)*(N-1)).divexact(g)
  return g,q,qoverg

@nprofile.profile
def decrypt(d,c,expandedseckey):
  g,q,qoverg = expandedseckey
  n = len(d)
  N = 2^n
  R = ring.ring(d)

  r = R((c*x + q//2)//q for x in qoverg.c)

  # XXX: limit precision to -max(mmax)...max(mmax)
  rg = r * g
  decrypt = (c - rg.c[0],) + tuple(-rg.c[j] for j in range(1,N))

  return decrypt
