import nprofile
profiling = ['unitsymbols','unitloginverse','mainwork_quad','hsymbols','exponents_square','buildsquare','squareroot','idealsqrtshorten','exponents_rounded','times_units','shorten','mainwork']
import field
import units
import div
from memoized import memoized

@nprofile.profile
@memoized
def unitloginverse(d):
  n = len(d)
  N = 2^n
  S = units.generators_mod_torsion(d)
  M = matrix(QQ,[(1,)*N] + [u.approxlog for u in S])
  return 1/M

def approxlog(d,g):
  n = len(d)
  N = 2^n

  prec = 2*units.logbits
  while True:
    if min(d) > 0:
      K = RealIntervalField(prec)
    else:
      K = ComplexIntervalField(prec)
    gsize = sqrt(K(sum(gj*gj for gj in g.numer.c)))
    gK = [K(gj)/gsize for gj in g.numer.c]
    s = [sqrt(K(di)) for di in d]
  
    v = div.hadamard(n,N,s,gK)
    v = [abs(vj) for vj in v]
    try:
      if any(0 in vj for vj in v):
        raise Exception('computing log of 0, need more precision')
      v = [log(vj) for vj in v]
      v = [vj*2^units.logbits for vj in v]
      if any(vj.absolute_diameter() > 1 for vj in v):
        raise Exception('not happy with log precision')
      v = [ZZ(round(vj.center())) for vj in v]
      v = [vj/2^units.logbits for vj in v]
      return vector(QQ,v)
    except:
      prec *= 2

@nprofile.profile
def exponents_rounded(d,g,scale):
  v = approxlog(d,g)
  e = v * unitloginverse(d) / scale
  return [-ZZ(round(ej)) for ej in e[1:]]

@nprofile.profile
def times_units(d,g,e):
  S = units.generators_mod_torsion(d)
  # better if h is small:
  #   field.powerprod(d,[g] + [u.element for u in S],[[1] + e])[0]
  return g * field.powerprod(d,[u.element for u in S],[e])[0]

@nprofile.profile
def shorten(d,g):
  e = exponents_rounded(d,g,1)
  return times_units(d,g,e)

@nprofile.profile
def mainwork_quad(d,q,s):
  D = d[0]
  K = units.quadfield(D)
  u0,u1,denom,approxlog = units.fundamentalunit(D)
  u0,u1 = u0/denom,u1/denom
  O = K.ring_of_integers()
  I = O.ideal([q,K.gens()[0]-s[0]])
  g = I.gens_reduced()
  if len(g) > 1: raise Exception('ideal is not principal')
  g = g[0]
  g0,g1 = list(g)
  s = lcm(QQ(g0).denom(),QQ(g1).denom())
  return field.field(d)((ZZ(g0*s),ZZ(g1*s)),s)

@nprofile.profile
@memoized
def unitsymbols(n,S,low,high):
  return matrix(GF(2),[u.symbols_log(low,high) for u in S])

@nprofile.profile
def hsymbols(n,h,rank):
  s,low,high = units.symbols_log(n,[h],rank)
  return s[0],low,high

@nprofile.profile
def exponents_square(n,d,M,v):
  return M.solve_left(vector(GF(2),v))

@nprofile.profile
def buildsquare(n,d,S,h,e):
  return h * field.powerprod(d,[u.element for u in S],[e])[0]

@nprofile.profile
def squareroot(n,d,uh):
  return uh.sqrt()

@nprofile.profile
def idealsqrtshorten(n,N,d,h):
  S = units.generators(d)
  v,low,high = hsymbols(n,h,len(S))
  M = unitsymbols(n,S,low,high)

  e1 = exponents_square(n,d,M,v)
  e1 = [ZZ(ej) for ej in e1]

  v = approxlog(d,h)
  v += vector(QQ,e1) * matrix(QQ,[u.approxlog for u in S])

  e2 = v * unitloginverse(d) / 2
  e2 = [-ZZ(round(ej)) for ej in e2]

  gg = buildsquare(n,d,S,h,[e1[j]+e2[j]*2 for j in range(len(S))])
  return squareroot(n,d,gg)

@memoized
def mainwork_internal(n,N,d,q,s):
  if n == 0:
    return field.field(d)((q,),1)
  if n == 1:
    return shorten(d,mainwork_quad(d,q,s))

  K = field.field(d)
  d1 = d[:-1]
  d2 = d[:-2] + d[-1:]
  d3 = d[:-2] + (d[-2]*d[-1],)
  s1 = s[:-1]
  s2 = s[:-2] + s[-1:]
  s3 = s[:-2] + (s[-2]*s[-1],)
  g1 = K(mainwork_internal(n-1,N//2,d1,q,s1))
  g2 = K(mainwork_internal(n-1,N//2,d2,q,s2))
  g3 = K(mainwork_internal(n-1,N//2,d3,q,s3))
  g3 = g3.conj()

  h = (g1*g2).divexact(g3)

  return idealsqrtshorten(n,N,d,h)

@nprofile.profile
def mainwork(n,N,d,q,s):
  return mainwork_internal(n,N,d,q,s)

def ideals(d):
  n = len(d)
  N = 2^n
  class I:
    def __init__(f,q,s):
      q = ZZ(q)
      if q % 2 == 0:
        raise Exception('only odd ideals are supported')
      for j in range(n):
        if gcd(d[j],q) != 1:
          raise Exception('only ideals coprime to d are supported')
      s = tuple(ZZ(sj) % q for sj in s)
      if len(s) != n:
        raise Exception('%s does not have length %s' % (s,n))
      for j in range(n):
        if Mod(s[j]^2-d[j],q) != 0:
          raise Exception('invalid image of square root of %s' % d[j])
      f.q = q
      f.s = s
    def __repr__(f):
      return 'ideal.ideals(%s)(%s,%s)' % (d,f.q,f.s)
    def generator(f):
      return mainwork(n,N,d,f.q,f.s)
  return I
