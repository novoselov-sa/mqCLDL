logbits = 100

import nprofile
import numpy as np
profiling = ['fundamentalunit','symbols_log','kernel','squares','adjoin_sqrt','shortening_matrix','shorten','generators_mod_torsion','generators']
import field
import trees
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
  myList = []
  printing = False
  iti = 0
  currentprime = []
  for line in f:
    if line[0] == 'F':
      line = line[:-1]
      spl = line.split(' ')
      field_id = spl[2]
      field = []
      for elt in field_id.split('_'):
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
          # Prime ideal is given in two gens representation (p, elts).
          # The list pows contains pointers to prime ideals the field extension (for root field it contains zeroes).
          p = int(currentprime[0])
          elts = currentprime[1][:-1]
          pows = [int(l) for l in line[1:-2].split(', ')]
          names = field_id.split("_")
          myList += [P(p, elts, pows, names=names)]
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
    print("d done:", d)
  return primes

def factorsb(d):
  n = len(d)
  N = 2^n
  K = field.field(d)
  class P:
    def __init__(f,*args, **kwargs):
      f.names = kwargs.get("names") # names of variables in elts

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
    def to_sage(f, food = None):
      if food == None:
        food = trees.get_food()
      K_sage = K.sage(food)
      return K_sage.ideal(f.prime, K_sage(f.elts))
  return P

# Checks whether an integer b is smooth with respect to the factor base FB. 
def is_smooth(b, FB):
  r = b
  for i in range(len(FB)):
    #if Mod(r, FB[i].prime) == 0:
    r = r / (QQ(FB[i].prime) ^ valuation(QQ(r), QQ(FB[i].prime)))
    if r == 1:
      return True
  return False

# Factor Sage's ideal of the field K over given factor base.
# Returns vector of exponents.
def factor(K, I, FB):
  N_I = pari.idealnorm(K.pari_nf(), I)
  e = []
  for i in range(len(FB)):
    if Mod(N_I, FB[i].prime) != 0:
      e.append(0)
    else:
      P = K.ideal(FB[i].prime, K(FB[i].elts))
      v = ZZ(pari.idealval(K.pari_nf(), I, P.pari_prime()))
      e.append(v)
  return e

# Given a vector of exponets for ideals in factor base computes explicit product of these ideals.
# Very slow, use for testing purposes only.
def prod_ideal(K, FB, e):
  assert len(FB) == len(e), f"Wrong length of factor base or input vector, {len(FB)} != {len(e)}"
  return prod(K.ideal(FB[i].prime, K(FB[i].elts))^int(e[i]) for i in range(len(FB)))
