logbits = 300
#logbits = 30

import nprofile
profiling = ['fundamentalunit','symbols_log','kernel','squares','adjoin_sqrt','shortening_matrix','shorten','generators_mod_torsion','generators']
import field
from memoized import memoized

@memoized
# XXX: need to analyze memory cost and time benefit of this memoization
def quadfield(D, name='a'):
  return QuadraticField(D, name)

@nprofile.profile
@memoized
def fundamentalunit(D):
  assert D > 0, "Imaginary quadratic fields doesn't have fundamental units"
  K = quadfield(D)
  u = K.units(proof=False)[0]
  u0,u1 = list(u)
  u0,u1 = abs(u0),abs(u1)
  s = lcm(QQ(u0).denom(),QQ(u1).denom())
  R = RealField(2*logbits)
  approxlog = R(log(R(u0)) + log(R(1) + (R(u1)/R(u0))*sqrt(R(D)))) # same sign so no cancellation
  #print "aplog1", approxlog
  #print R(2^logbits)*approxlog
  #print R(2^logbits)
  #print "round", floor(R(2^logbits*approxlog))
  #print ZZ(floor(R(2^logbits)*approxlog))
  approxlog = ZZ(ceil(R(2^logbits)*approxlog))/R(2^logbits)
  return ZZ(u0*s),ZZ(u1*s),s,approxlog

@memoized
def units(d):
  n = len(d)
  N = 2^n
  K = field.field(d)
  class U:
    gens = d
    def __init__(f,*args):
      #print 'args', args
      if len(args) == 0:
        if n == 1: # U() for quadratic is fundamental unit
          u0s,u1s,s,approxlog = fundamentalunit(d[0])
          f.element = K(((u0s,u1s)),s)
          f.powers = (0,1)
          f.approxlog = (approxlog,-approxlog)
          return
      if len(args) == 1:
        g = args[0]
        if g in [-1,1]: # U(1) is 1; U(-1) is -1
          f.element = K((g,)+(0,)*(N-1),1)
          f.powers = (0,)*N
          f.approxlog = (0,)*N
          return
        if n >= 1 and  g.__class__.gens == units(d[:-1]).gens:
          f.element = K(g.element)
          f.powers = g.powers + (0,)*(N//2)
          f.approxlog = g.approxlog + g.approxlog
          return
        if n >= 2 and g.__class__.gens == units(d[:-2] + d[-1:]).gens:
          f.element = K(g.element)
          f.powers = g.powers[:N//4] + (0,)*(N//4) + g.powers[N//4:] + (0,)*(N//4)
          f.approxlog = g.approxlog[:N//4]*2 + g.approxlog[N//4:]*2
          return
        if n >= 2 and g.__class__.gens == units(d[:-2] + (d[-2]*d[-1],)).gens:
          f.element = K(g.element)
          f.powers = g.powers[:N//4] + (0,)*(N//2) + g.powers[N//4:]
          f.approxlog = g.approxlog + g.approxlog[N//4:] + g.approxlog[:N//4]
          return
      if len(args) == 3: # XXX: trusting caller to provide suitable values
        f.element,f.powers,f.approxlog = args
        return
      raise Exception('not known how to initialize units(%s)(%s)' % (str(d),args))
    def __repr__(f):
      return 'units.units(%s)(%s,%s,%s)' % (d,f.element,f.powers,f.approxlog)
    def __mul__(f,g):
      helement = f.element * g.element
      hpowers = tuple(f.powers[j] + g.powers[j] for j in range(N))
      happroxlog = tuple(f.approxlog[j] + g.approxlog[j] for j in range(N))
      return U(helement,hpowers,happroxlog)
    #def __div__(f,g): # Python 2
    def __truediv__(f,g): # Python 3
      helement = f.element.divexact(g.element)
      hpowers = tuple(f.powers[j] - g.powers[j] for j in range(N))
      happroxlog = tuple(f.approxlog[j] - g.approxlog[j] for j in range(N))
      return U(helement,hpowers,happroxlog)
    def __pow__(f,e):
      if e == 0: return U(1)
      if e == 1: return f
      if e < 0: return U(1)/(f^(-e))
      g = f^(e // 2)
      g *= g
      if e & 1: g *= f
      return g
    def sqrt(f):
      helement = f.element.sqrt()
      hpowers = tuple(f.powers[j]/2 for j in range(N))
      happroxlog = tuple(f.approxlog[j]/2 for j in range(N))
      return U(helement,hpowers,happroxlog)
    def symbols2_log(f,low,high):
      return f.element.symbols2_log(low,high)
    def conj(f): # conjugate last entry of d
      helement = f.element.conj()
      hpowers = f.powers[:N//2] + tuple(-x for x in f.powers[N//2:])
      happroxlog = f.approxlog[N//2:] + f.approxlog[:N//2]
      return U(helement,hpowers,happroxlog)
  return U

def subsetprod(d,g,e):
  U = units(d)
  for gj in g:
    if gj.__class__.gens != U.gens:
      raise Exception('%s not in units(%s)' % (gj,d))
  n = len(d)
  N = 2^n
  helement = field.subsetprod(d,[gj.element for gj in g],e)
  hpowers = [tuple(sum(gj.powers[k] for gj,eij in zip(g,ei) if eij) for k in range(N)) for ei in e]
  happroxlog = [tuple(sum(gj.approxlog[k] for gj,eij in zip(g,ei) if eij) for k in range(N)) for ei in e]
  return tuple(U(e,p,a) for e,p,a in zip(helement,hpowers,happroxlog))

def powerprod(d,g,e):
  U = units(d)
  for gj in g:
    if gj.__class__.gens != U.gens:
      raise Exception('%s not in units(%s)' % (gj,d))
  n = len(d)
  N = 2^n
  helement = field.powerprod(d,[gj.element for gj in g],e)
  hpowers = [tuple(sum(eij*gj.powers[k] for gj,eij in zip(g,ei) if eij) for k in range(N)) for ei in e]
  happroxlog = [tuple(sum(eij*gj.approxlog[k] for gj,eij in zip(g,ei) if eij) for k in range(N)) for ei in e]
  return tuple(U(e,p,a) for e,p,a in zip(helement,hpowers,happroxlog))

@nprofile.profile
def kernel(n,M):
  return M.left_kernel().basis_matrix().rows()

@nprofile.profile
def squares(n,d,S,e):
  # XXX: figure out bottlenecks
  return subsetprod(d,S,e)
  # return [prod(Sj for Sj,eij in zip(S,ei) if eij == 1) for ei in e]

@nprofile.profile
def adjoin_sqrt(n,S,E):
  return S + [sqrt(Ei) for Ei in E]

# return matrix that turns vectors into short nonzero vectors
@nprofile.profile
def shortening_matrix(n,N,S):
  M = matrix(ZZ,[[ZZ(i==j) for i in range(len(S))]
                 + [ZZ(v*2^(logbits+n-1)) for v in S[j].approxlog]
                 for j in range(len(S))])
  M = M.LLL()
  return [Mi[:len(S)] for Mi in M if not Mi[len(S):].is_zero()]

@nprofile.profile
def shorten(n,d,S,L):
  return powerprod(d,S,L)
  # return [prod(Sj^ZZ(ej) for Sj,ej in zip(S,e) if ej != 0) for e in L]

@nprofile.profile
def symbols2_log(n,S,rank):
  low,high = 0,rank + 64 # XXX: 64 is overkill
  while True:
    try:
      return [u.symbols2_log(low,high) for u in S],low,high
    except:
      low,high = high,high + (high - low)

@memoized
# number of roots of unity, and generator
def torsion(d):
  n = len(d)
  N = 2^n

  if n == 1 and (-d[0]).is_square():
    K = field.field(d)
    felement = K((0,1),sqrt(-d[0]))
    fpowers = (0,0)
    fapproxlog = (0,0)
    return 4,units(d)(felement,fpowers,fapproxlog)

  if n == 1 and (-3*d[0]).is_square():
    K = field.field(d)
    r = sqrt(-3*d[0])
    felement = K((r,3),2*r)
    fpowers = (0,0)
    fapproxlog = (0,0)
    return 6,units(d)(felement,fpowers,fapproxlog)

  if n <= 1 or min(d) >= 0:
    return 2,units(d)(-1)

  U = units(d)
  w1,t1 = torsion(d[:-1])
  w2,t2 = torsion(d[:-2] + d[-1:])
  w3,t3 = torsion(d[:-2] + (d[-2]*d[-1],))
  t1 = U(t1)
  t2 = U(t2)
  t3 = U(t3) # no need for conj here
  # now <t1,t2,t3> contains squares of all roots of unity

  w = lcm([w1,w2,w3])
  if w == w1:
    t = t1
  elif w == w2:
    t = t2
  elif w == w3:
    t = t3
  else:
    t = U(1)
    # can replace this by factoring into coprimes
    # but w is always small, so factor() is fine
    for p,e in factor(w):
      if w1 % p^e == 0:
        t *= t1^(w1//p^e)
      elif w2 % p^e == 0:
        t *= t2^(w2//p^e)
      elif w3 % p^e == 0:
        t *= t3^(w3//p^e)
      else:
        raise Exception('internal error in torsion')
  
  # now <t> contains squares of all roots of unity
  if all(s == 0 for s in symbols2_log(n,[t],N)[0][0]):
    w,t = 2*w,sqrt(t)

  return w,t

# generator of cyclic group of roots of unity
@memoized
def generator_of_torsion(d):
  return torsion(d)[1]

# generators of full unit group mod torsion
@memoized
def generators_internal(d):
  n = len(d)

  if n == 0:
    return ()

  if n == 1:
    D = d[0]
    if D < 0: return ()
    return (units(d)(),)

  N = 2^n
  freerank = N - 1
  if min(d) < 0: freerank = N//2 - 1
  #print(f'freerank: {freerank}')
  
  rank = freerank + 1

  U = units(d)
  U1 = generators_internal(d[:-1])
  U2 = generators_internal(d[:-2] + d[-1:])
  U3 = generators_internal(d[:-2] + (d[-2]*d[-1],))
  U1 = [U(u) for u in U1]
  U2 = [U(u) for u in U2]
  U3 = [U(u).conj() for u in U3]

  root = [generator_of_torsion(d)]
  #S = uniq(sorted(root + U1 + U2 + U3))
  S = list(set(root + U1 + U2 + U3))

  M = matrix(GF(2),symbols2_log(n,S,rank)[0])
  e = kernel(n,M)
  E = squares(n,d,S,e)
  S = adjoin_sqrt(n,S,E)
  L = shortening_matrix(n,N,S)
  if len(L) != freerank:
    raise Exception('internal error; maybe logbits not high enough?')
  T = shorten(n,d,S,L)

  #print "T", T

  return tuple(T)

@nprofile.profile
def generators_mod_torsion(d):
  return generators_internal(d)

@nprofile.profile
def generators(d):
  return (generator_of_torsion(d),) + generators_mod_torsion(d)

def compresslog(d,approxlog):
  n = len(d)
  N = 2^n
  for i in range(n):
    if d[i] < 0:
      result = []
      for j in range(1,N):
        if j & (1 << i) == 0:
          result += [2*approxlog[j]]
      return tuple(result)
  return approxlog[1:]

def approxregulator(d, prec=500):
  U = generators_mod_torsion(d)
  RP = RealField(prec)
  M = matrix(RP,[compresslog(d,u.approxlog) for u in U])
  #M = matrix(RP,[u.approxlog[1:] for u in U])
  return abs(det(M))

def mqgenerators_mod_torsion(d):
  n = len(d)
  N = 2^n
  K = field.field(d)
  U = units(d)
  result = []
  for j in range(1,N):
    D = prod(d[i] for i in range(n) if j & (1 << i))
    if D > 0:
      u0s,u1s,s,approxlog = fundamentalunit(D)
      numer = [0]*N
      numer[0] = u0s
      numer[j] = u1s
      element = K(tuple(numer),s)
      print(D, element)
      powers = tuple(ZZ(k == j) for k in range(N))
      approxlog = tuple(approxlog*(-1)^ZZ(k & j).popcount() for k in range(N))
      result += [U(element,powers,approxlog)]
  return tuple(result)

def approxmqregulator(d):
  U = mqgenerators_mod_torsion(d)
  M = matrix(RR,[compresslog(d,u.approxlog) for u in U])
  return abs(det(M))

def compresspowers(d,powers):
  if all(di > 0 for di in d):
    return powers[1:]
  n = len(d)
  N = 2^n
  return tuple(powers[j] for j in range(1,N) if prod(d[i] for i in range(n) if j & (1 << i)) > 0)

def mqindex(d):
  U = generators_mod_torsion(d)
  M = matrix(QQ,[compresspowers(d,u.powers) for u in U])
  return ZZ(1/abs(det(M)))
