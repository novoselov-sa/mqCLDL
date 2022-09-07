logbits = 100

import nprofile
import numpy as np
profiling = ['fundamentalunit','symbols_log','kernel','squares','adjoin_sqrt','shortening_matrix','shorten','generators_mod_torsion','generators']
import field
from memoized import memoized

def quadfield(D, name='a'):
  return QuadraticField(D, name)

def belabas(d):
        return prod(d)^(2^(len(d)-1))

def matrixmaker(kd,roots,p):
        dI = [1]
        unit = kd*[0]
        unit[0] = 1
        dMat = [unit]
        dBasis = [1]
        for i in range(len(roots)):
                for j in range(len(dI)):
                        dBasis += [(1+basis[i])/2*dBasis[j]]
                        dI += [roots[i]*dI[j]]
                        rl = dI[-1].polynomial().list()
                        while len(rl) < kd: rl += [0]
                        dMat += [rl]
        return matrix(ZZ,dMat)

#turns sage primes into this code framework
def sage_prime(d, p, elt, name='a'):
  p = ZZ(p)
  K = quadfield(d, name)
  a = K.gen()
  root = (p- ZZ(elt - 1/2 - a/2))%p
  return [p,[p,elt],1,[root],matrix(GF(p),[[1,0],[root,0]])]

def load_primes(d, file, food = None):
  if food == None:
    food = {"a": 5, "b": 13, "c": 17, "d": 29, "e": 37, "f": 41, "g": 53, "h": 61}
  #foodinv = {5: "a", 13: "b", 17: "c" , 29: "d" , 37: "e", 41: "f", 53: "g", 61: "h" }
  foodinv = dict((k,v) for v,k in food.iteritems())
  P = factorsb(d)
  #print d
  #if len(d) ==1:
  #  file = "a_"
  #  if d == (5,):
  #    file += "b"
  #  elif d[0]%5 == 0:
  #    file += ''.join([foodinv[i] for i in zip(*list(ZZ(d[0]/5).factor()))[0]])
  #  else:
  #    file += ''.join([foodinv[i] for i in zip(*list(d[0].factor()))[0]])
  #print file
  f = open("trees/" + file + ".christine","r")
  #if d[0]> 0:
  #  f = open("trees/" + file + ".christine","r")
  #else:
  #  f = open("trees/imaginary/" + file + ".christine","r")
  myList = []
  printing = False
  iti = 0
  currentprime = []
  for line in f:
    if line[0] == 'F':
      line = line[:-1]
      spl = line.split(' ')
      lst = spl[2]
      field = []
      for elt in lst.split('_'):
        field += [prod(int(food[k]^elt.count(k)) for k,_ in food.iteritems())]
      if tuple(field) == d: printing = True
      else: printing = False
    else:
      if printing:
        #very f-in ugly
        if iti < 2:
          currentprime += [line]
          iti += 1 
        else:
          iti = 0
          myList += [P(int(currentprime[0]),currentprime[1][:-1],[int(l) for l in line[1:-2].split(', ')])]
          currentprime = []
  f.close()
  return myList

def load_primes2(d):
     P = factorsb(d) 
     f = open("primes.txt","r")
     print("loading", P)
     myList = []
     printing = False
     for line in f:
       if line[0] == 'F':
         line = line[:-1]
         spl = line.split(' ')
         list = spl[1]
         field = [int(l) for l in list.split(',')]
         print("field: ", field)
         print(d)
         if tuple(field) == d: printing = True
         else: printing = False
       else:
         if printing:
            ln = line.splitlines()
            myList += [P(int(ln[0]), ln[1][:-1], [int(l) for l in ln[2][1:-1].split(', ')])]
     f.close()
     return myList

def elts_from_kernel(M,d):
  n = len(d)
  N = 2^n
  K = field.field(d)
  form = field.formmaker(d)
  elts = []
  for h in M.kernel().basis():
    h1 = np.array(h).nonzero()[0]
    elt = K(N*[0],1)
    for h2 in h1:
       elt += (K([h[h2]] + (N-1)*[0],1)*form[h2])
    elts += [elt]
  return elts

def primes_internal(p, d):
  n = len(d)
  N = 2^n
  primes = []
  if n == 1:
    print("d", d)
    P = primesid(tuple(d))
    primes = [P(fp) for fp in fundamentalprimes(p,d[0])]
    print("d done:", d[0])
    return primes
  else:
    K = field.field(tuple(d))
    print("d", d)
    form = field.formmaker(d)
    P = primesid(tuple(d))
    F = GF(p^2)
    f = F.gen()
    R.<x> = F[]
    poly = R(x^2 -x - ZZ((d[-1]-1)/4))
    print(poly.roots())
    primecur = primes_internal(p,d[:-1])
    for curp in primecur:
      pm = curp.matrix
      for roots in poly.roots():
        rl = roots[0].polynomial().list()
        while len(rl) < 2: rl += [0]
        rl2 = (f*roots[0]).polynomial().list()
        while len(rl2) < 2: rl2 += [0]
        rm = matrix(GF(p),[rl,rl2])
        newm = pm.stack(pm*rm)
        ib = elts_from_kernel(newm,tuple(d))
        primes.append(P(p,ib,roots[1],curp.roots + [roots[0]],newm))
    print("d done:", d)
  return primes

#there is some switching of variables that is wrong
def primes_right(p, d):
  n = len(d)
  N = 2^n
  primes = []
  if n == 1:
    print "d", d
    P = primesid(tuple(d))
    primes = [P(fp) for fp in fundamentalprimes(p,d[0])]
    print "d done:", d[0]
    return primes
  else:
    K = field.field(tuple(d))
    print "d", d
    form = field.formmaker(d)
    P = primesid(tuple(d))
    F = GF(p^2)
    f = F.gen()
    R.<x> = F[]
    poly = R(x^2 -x - ZZ((d[-2]-1)/4))
    primecur = primes_right(p,d[:-2] + d[-1:])
    for curp in primecur:
      pm = curp.matrix
      for roots in poly.roots():
        rl = roots[0].polynomial().list()
        while len(rl) < 2: rl += [0]
        rl2 = (f*roots[0]).polynomial().list()
        while len(rl2) < 2: rl2 += [0]
        rm = matrix(GF(p),[rl,rl2])
        newm = matrix(list(sum(map(list, zip(pm, pm*rm)), [])))
        ib = elts_from_kernel(newm,tuple(d))
        primes.append(P(p,ib,roots[1],curp.roots + [roots[0]],newm))
    print "d done:", d
  return primes

def factorsb(d):
  n = len(d)
  N = 2^n
  K = field.field(d)
  class P:
    def __init__(f,*args):
      if len(args) == 1: # list instead of seperates
        c = args[0]
        [f.prime,f.elts,f.powers] = c
        return
      if len(args) == 3: # XXX: trusting caller to provide suitable values
        f.prime,f.elts,f.powers = args
        return
      raise Exception('not known how to initialize units(%s)(%s)' % (str(d),args))
    def __repr__(f):
      return 'fb.prime(%s)(%s,%s,%s)' % (d,f.prime,f.elts,f.powers)
    def __mul__(f,g):
      if f == 0: return ZZ(0)
      if f == 1: return g
      #helement = f.element * g.element
      #hpowers = tuple(f.powers[j] + g.powers[j] for j in range(N))
      #happroxlog = tuple(f.approxlog[j] + g.approxlog[j] for j in range(N))
      #return U(helement,hpowers,happroxlog)
    def __div__(f,g):
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
    def symbols_log(f,low,high):
      return f.element.symbols_log(low,high)
    def conj(f): # conjugate last entry of d
      helement = f.element.conj()
      hpowers = f.powers[:N//2] + tuple(-x for x in f.powers[N//2:])
      happroxlog = f.approxlog[N//2:] + f.approxlog[:N//2]
      return U(helement,hpowers,happroxlog)
  return P

def subsetprod(d,g,e):
  U = units(d)
  for gj in g:
    if gj.__class__ != U:
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
    if gj.__class__ != U:
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
def symbols_log(n,S,rank):
  low,high = 0,rank + 64 # XXX: 64 is overkill
  while True:
    try:
      return [u.symbols_log(low,high) for u in S],low,high
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
  if all(s == 0 for s in symbols_log(n,[t],N)[0][0]):
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
  rank = freerank + 1

  U = units(d)
  U1 = generators_internal(d[:-1])
  U2 = generators_internal(d[:-2] + d[-1:])
  U3 = generators_internal(d[:-2] + (d[-2]*d[-1],))
  print U3[0]
  U1 = [U(u) for u in U1]
  U2 = [U(u) for u in U2]
  U3 = [U(u).conj() for u in U3]

  root = [generator_of_torsion(d)]
  S = uniq(sorted(root + U1 + U2 + U3))

  M = matrix(GF(2),symbols_log(n,S,rank)[0])
  e = kernel(n,M)
  E = squares(n,d,S,e)
  S = adjoin_sqrt(n,S,E)
  L = shortening_matrix(n,N,S)
  if len(L) != freerank:
    raise Exception('internal error; maybe logbits not high enough?')
  T = shorten(n,d,S,L)

  print("T", T)

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

def approxregulator(d):
  U = generators_mod_torsion(d)
  M = matrix(RR,[compresslog(d,u.approxlog) for u in U])
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
