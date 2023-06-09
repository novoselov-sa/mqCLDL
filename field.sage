import ring
import char
import polynomial_ring
import ideal
import fb
import relations
import units
import sunits
import clgp
import field_compact
import trees

from memoized import memoized

@memoized
def field_knownniceness(d,nice):
  n = len(d)
  N = 2^n
  R = ring.ring(d)
  if nice:
    multiplier = N
  else:
    multiplier = N * abs(prod(d))
  class K:
    gens = d
    multiplieroverN = multiplier / N
    def __init__(f,*args):
      if len(args) == 2:
        numer,denom = args
        if numer.__class__ == R: numer = numer.c
        numer = tuple(ZZ(cj) for cj in numer)
        denom = ZZ(denom)
        g = gcd(numer + (denom,))
        numer = tuple(ZZ(cj/g) for cj in numer)
        denom = ZZ(denom/g)
        f.numer = R(numer)
        f.denom = denom
        return
      if len(args) == 1:
        g = args[0]
        if n >= 1 and g.__class__ == field(d[:-1]):
          f.numer = R(g.numer.c + (0,)*(N//2))
          f.denom = g.denom
          return
        if n >= 2 and g.__class__ == field(d[:-2] + d[-1:]):
          gc = g.numer.c
          f.numer = R(gc[:N//4] + (0,)*(N//4) + gc[N//4:] + (0,)*(N//4))
          f.denom = g.denom
          return
        if n >= 2 and g.__class__ == field(d[:-2] + (d[-2]*d[-1],)):
          gc = g.numer.c
          f.numer = R(gc[:N//4] + (0,)*(N//2) + gc[N//4:])
          f.denom = g.denom
          return
    def __repr__(f):
      return 'field.field(%s)(%s,%s)' % (d,f.numer,f.denom)
    def __eq__(f,g):
      if g.__class__ != f.__class__: raise Exception('do not know how to compare %s to %s' % (f,g))
      return f.numer == g.numer and f.denom == g.denom
    def __ne__(f,g):
      return not f == g
    def is_nonzero(f):
      return f.numer.is_nonzero()
    def __add__(f,g):
      if g.__class__ == Integer:
        return f + K.from_ZZ(g)
      if g.__class__ != f.__class__: raise Exception('do not know how to compute %s + %s' % (f,g))
      return K(f.numer * g.denom + g.numer * f.denom,f.denom * g.denom)
    def __neg__(f):
      return K(-f.numer,f.denom)
    def __sub__(f,g):
      if g.__class__ == Integer:
        return f - K.from_ZZ(g)
      if g.__class__ != f.__class__: raise Exception('do not know how to compute %s - %s' % (f,g))
      return K(f.numer * g.denom - g.numer * f.denom,f.denom * g.denom)
    def __mul__(f,g):
      if g.__class__ == Integer:
        return f * K.from_ZZ(g)
      if g.__class__ != f.__class__: raise Exception('do not know how to compute %s * %s' % (f,g))
      return K(f.numer * g.numer,f.denom * g.denom)
    def __truediv__(f,g):
      #FIXME: slow, we use methods of Sage since we have only exact division.
      f0 = f.to_sage(K.sage())
      g0 = g.to_sage(K.sage())
      return K.from_sage(f0/g0)
    def square(f):
      return K(f.numer.square(),f.denom^2)
    def sqrt(f):
      return K((f.numer * (multiplier * multiplier * f.denom)).sqrt(),multiplier * f.denom)
    def divexact(f,g):
      hdenom = multiplier
      top = f.numer * (hdenom * g.denom)
      bot = g.numer * f.denom
      hnumer = top.divexact(bot)
      return K(hnumer,hdenom)
    def conj(f): # conjugate last entry of d
      fc = f.numer.c
      return K(fc[:N//2] + tuple(-x for x in fc[N//2:]),f.denom)
    def normonestep(f,t):
      h = f.numer.normonestep(t)
      K1 = field_internal(d[:t] + d[t+1:],nice)
      return K1(h,f.denom^2)
    def absnorm(f):
      return f.numer.absnorm() / f.denom^N
    def norm_size(f):
      return abs(f.absnorm().numerator()) + abs(f.absnorm().denominator())
    def getd(f):
      return d
    def get_qs(f,low,high):
      return char.altsymbols(d,f.numer.c,f.denom,low,high)
    def const(f):
      return f.numer[0] / f.denom
    def symbol(f,t):
      return char.symbol(d,f.numer.c,f.denom,t)
    def symbols(f,low,high):
      return char.symbols(d,f.numer.c,f.denom,low,high)
    def symbols_log(f,low,high):
      s = char.symbols(d,f.numer.c,f.denom,low,high)
      if 0 in s: raise Exception('zero symbol')
      return tuple(0 if si == 1 else 1 for si in s)
    def symbols2_log(f,low,high):
      s = char.symbols2(d,f.numer.c,f.denom,low,high)
      if 0 in s: raise Exception('zero symbol')
      return tuple(0 if si == 1 else 1 for si in s)
    def valuation(f, P):
      return f.to_sage(K.sage()).valuation(P)
    def approxlog(f, prec=None):
      return ideal.approxlog(d, f, prec=prec)
    def lognorm(f, p=2, prec=None):
      if prec == None:
        prec = units.logbits
      R = RealField(prec)
      v = vector(f.approxlog(prec=prec))
      return R(v.norm(p=2))
    def __pow__(f, a):
      assert type(a) == Integer
      if a < 0:
        return K.one() / f ^ abs(a)
      if a == 0:
        return K.one()
      return sage.arith.power.generic_power(f, a)

    def reduce_mod_units(f, rollback=False):
      r = ideal.shorten(d, f)
      if rollback and r.bit_size() >= f.bit_size():
        return f
      #print(f"\n\nf = {f.to_sage(K.sage())} reduced to r = {r.to_sage(K.sage())}")
      #print(f"-> f.bit_size() = {f.bit_size()} reduced to r.bit_size() = {r.bit_size()}")
      #print(f"-> f.lognorm() = {f.lognorm()} reduced to r.lognorm() = {r.lognorm()}")
      return r

    def reduce_mod_sunits(f, SU, mod_units=True):
      '''
      Given element f of the field K returns pair (f/s, s) where s is an S-unit in power-product representation and ||Log_S(f/s)||_2 <= ||Log_S(f)||_2.
      '''
      A = SU.relations_matrix()
      FB = SU.factor_base(sage=True)
      FB0 = SU.factor_base(sage=False)

      v = [f.valuation(FB[j]) for j in range(len(FB))]
      if fb.is_smooth(f.absnorm(), FB0):
        return (K.one(), SU.from_prime_prod(v).compute_roots(ideal_sqrt=mod_units))

      w = clgp.approx_CVP(d, v)
      assert len(w) == len(SU.gens())

      wA = vector(w) * A

      # reduction failed, resulting vector has larger norm
      if RR((vector(v)-wA).norm()) >= RR(vector(v).norm()):
        return (f, SU.one())

      s = prod([SU.gens()[i]^w[i] for i in range(len(w))], SU.one()).compute_roots(ideal_sqrt=mod_units)

      if mod_units:
        s = s.reduce_mod_units(rollback=True)

      KC = field_compact.field_compact(d)

      #h,r,s = KC(f).rdiv(s, l=2, bits=True, log_sunits=True, FB = FB)
      h,r,s = KC(f).rdiv(s, l=2, bits=True, mod_units=mod_units)

      assert len(h.elements) == 1 and len(h.powers) == 1 and h.powers[0] == 1
      h = h.elements[0]

      if r != KC.one() and (not vector(r.powers).is_zero()):
        print(f"[red_mod_sunits_incomplete({len(r.elements)})]")

      return (h, s)

    def mul_mod_sunits(f, g, SU, mod_units=True):
      return (f*g).reduce_mod_sunits(SU, mod_units=mod_units)

    def pow_mod_sunits(f, a, SU, mod_units=True):
      assert type(a) == Integer
      if a < 0:
        #return (K.one() / f) ^ abs(a)
        return (K.one() / f).pow_mod_sunits(abs(a), SU, mod_units=mod_units)
      if a == 0:
        return (K.one(), SU.one())
      if a == 1:
        #return (f, SU.one())
        return f.reduce_mod_sunits(SU, mod_units=mod_units)
      #g = f^(a // 2)
      g,s = f.pow_mod_sunits(a // 2, SU, mod_units=mod_units)
      #g *= g
      g,s0 = g.mul_mod_sunits(g, SU, mod_units=mod_units)
      s = s^2 * s0
      if a & 1:
        #g *= f
        g,s0 = g.mul_mod_sunits(f, SU, mod_units=mod_units)
        s *= s0
      return (g,s)

    def apply_aut(f, mu):
      """
      Applying an automorphism mu to element of the field.
      
      Automorphism mu is defined by a vector with elements in {0,1}.
      It is acts as d[i] -> (-1)^mu[i]*d[i].
      """
      assert len(mu) == n, f"Wrong size of vector representation of automorphism, {n} != {len(mu)}"

      J = [1]
      for i in range(n):
        J += [(-1)^mu[i]*g for g in J]
      nm = tuple(J[i]*f.numer.c[i] for i in range(N))
      return K(nm, f.denom)

    def to_sage(f, K = None, names='a', quotient=True):
      """
      Conversion the field element to Sage's polynomial ring QQ[x_1,...,x_n] (quotient = False),
      to the quotient ring QQ[x_1,...,x_n]/(x_1^2 - d_1, ..., x_n^2 - d_n) (quotient = True)), or
      given field K.
      """
      num = f.numer.to_sage(K = K, names=names, quotient=quotient)
      if K == None and quotient==True:
        rng = PolynomialRing(QQ, n, names=names[:n])
        I = num.parent().defining_ideal()
        rng = rng.quotient(I, names=names[:n])
        num = rng(num.lift())
      return num / QQ(f.denom)

    def minimal_polynomial(f):
      return polynomial_ring.minimal_polynomial(d, f)

    # Computation of field element mod p assuming that all d_i a squares in F_p.
    def evaluate_mod(f, p):
      d_sq = [Mod(di,p).sqrt() for di in d]
      I_sq = [1]
      for di in d_sq:
        I_sq += [di*i for i in I_sq]
      val = 0
      for i in range(len(I_sq)):
        if f.numer.c[i] != 0:
          val += f.numer.c[i] * I_sq[i]
      return Mod(val, p)/Mod(f.denom, p)

    # Evaluate element of the field mod p in the case when there are d_i s.t. d_i is not square in F_p.
    # This is only possible when element of the field is of the form sum_J (a_J*prod_{i in J}(d_i)),
    # where we have in each prod_{i in J}(d_i) the even number of d_i s.t. d_i is not square in F_p
    # (product of two non-squares is a square in F_p).
    def evaluate_mod_ext(f, p):
      F = GF(p)
      # Since square roots are defined up to sign, we fix for one of the variants each prod_{i in J}(d_i).
      # This guaranties that evaluation procedure will return the same result for each call.
      I_sq = field_sq_basis_mod(d, p)
      val = 0
      for i in range(len(I_sq)):
        if f.numer.c[i] != 0:
          val += f.numer.c[i] * I_sq[i]
      return F(val)/F(f.denom)
    
    def relnorm(f, d_M):
      M = field(d_M)
      if d[-1]*d[-2] == d_M[-1]:
        #raise Exception('unimplemented')
        # gc = g.numer.c
        # f.numer = R(gc[:N//4] + (0,)*(N//2) + gc[N//4:])
        # f.denom = g.denom
        f_s = f.apply_aut([0]*(n-2) + [1,1])
        g = f * f_s
        assert g.numer.c[N//4:N//4+N//2] == (0,)*(N//2), "error: element does not belongs to subfield [secondary subfield]"
        g_numer = g.numer.c[:N//4] + g.numer.c[N//4+N//2:]
        assert len(g_numer) == N//2
        return M(g_numer, g.denom)
      elif d[:-1] == d_M:
        #raise Exception('unimplemented')
        # f.numer = R(g.numer.c + (0,)*(N//2))
        # f.denom = g.denom
        f_s = f.apply_aut([0]*(n-1) + [1])
        g = f * f_s
        assert g.numer.c[N//2:] == (0,)*(N//2), "error: element does not belongs to subfield [primary subfield 1]"
        g_numer = g.numer.c[:N//2]
        assert len(g_numer) == N//2
        return M(g_numer, g.denom)
      elif d[:-2] + (d[-1],) == d_M:
        #raise Exception('unimplemented')
        #gc = g.numer.c
        #f.numer = R(gc[:N//4] + (0,)*(N//4) + gc[N//4:] + (0,)*(N//4))
        #f.denom = g.denom
        f_s = f.apply_aut([0]*(n-2) + [1, 0])
        g = f * f_s
        assert (g.numer.c[N//4:N//4+N//4] == (0,)*(N//4)) and (g.numer.c[-N//4:] == (0,)*(N//4)), "error: element does not belongs to subfield [primary subfield 2]"
        g_numer = g.numer.c[:N//4] + g.numer.c[N//4 + N//4:-N//4]
        assert len(g_numer) == N//2
        return M(g_numer, g.denom)
      else:
         raise Exception('unimplemented')
    def bit_size(f):
      return f.numer.bit_size() + RR(f.denom).abs().log(2).floor() + 1

    @staticmethod
    def abs_gen():
      # Return sqrt(d1) + ... + sqrt(d1).
      return field_abs_gen(d)
    
    @staticmethod
    def zero():
      return K(tuple([0]*N), ZZ(1))
    
    @staticmethod
    def one():
      return K(tuple([1] + [0]*(N-1)), ZZ(1))

    @staticmethod
    def from_ZZ(a):
      return K(tuple([a] + [0]*(N-1)), ZZ(1))
    
    @staticmethod
    def from_QQ(a):
      return K(tuple([a.numerator()] + [0]*(N-1)), a.denominator())
    
    @staticmethod
    def from_sage(a):
      return from_sage_internal(d, a, 0)

    @staticmethod
    def absolute_polynomial():
      return polynomial_ring.absolute_polynomial(d)

    @staticmethod
    def gens():
      # Return [sqrt(d1), ..., sqrt(d1)]
      return field_gens(d)

    @staticmethod
    def sq_basis_mod(p):
      # Return all possible products prod_{j in J}{sqrt(d_j)} as elements of Fp^2
      return field_sq_basis_mod(d, p)

    @staticmethod
    def discriminant():
      return mquad_disc(d)
    
    @staticmethod
    def degree():
      return N
    
    @staticmethod
    def idx():
      # Computes [O_K:Z[theta_K]]
      return mquad_idx(d)
    
    @staticmethod
    def sage(food = None, names = "a"):
      return field_sage(d, food = food, names = names)
    
    @staticmethod
    def random(bits = 16):
      return K(tuple(ZZ.random_element(1-2^bits,2^bits) for j in range(N)),ZZ.random_element(1,1+2^bits))

  return K

def field_internal(d,definitelynice):
  nice = True
  if not definitelynice:
    if not all(dj.is_squarefree() for dj in d):
      nice = False
    for j in range(len(d)):
      for i in range(j):
        if gcd(d[i],d[j]) != 1:
          nice = False
  return field_knownniceness(d,nice)

def field(d):
  return field_internal(d,False)

def subsetprod(d,g,e):
  K = field(d)
  for gj in g:
    if gj.__class__ != K:
      raise Exception('%s not in field(%s)',gj,d)
  hnumer = ring.subsetprod(d,[gj.numer for gj in g],e)
  hdenom = [ZZ.prod(gj.denom for gj,eij in zip(g,ei) if eij) for ei in e]
  return tuple(K(n,d) for n,d in zip(hnumer,hdenom))

def formmaker(d):
  K = field(d)
  n = len(d)
  N = 2^n
  list = [K([1] + (N-1)*[0],1)]
  for i in range(len(d)):
    Di = K([1] + (2**i - 1)*[0] + [1] + (N - 1 - 2**i)*[0],2)
    list2 = [Di*li for li in list]
    list.extend(list2)
  return list

def powerprod(d,g,e):
  K = field(d)
  for gj in g:
    if gj.__class__ != K:
      raise Exception('%s not in field(%s)',gj,d)
  hnumer = ring.scaledpowerprod(d,[gj.numer for gj in g],[gj.denom for gj in g],e,K.multiplieroverN)
  N = 2^len(d)
  return tuple(K(n,N * K.multiplieroverN) for n in hnumer)

@memoized
def field_abs_gen(d):
  n = len(d)
  N = 2^n
  K = field(d)
  a = [0] * N
  for i in range(n):
    a[2^i] = 1
  a = K(tuple(a), ZZ(1))
  return a

@memoized
def field_gens(d):
  n = len(d)
  N = 2^n
  K = field(d)
  res = []
  for i in range(n):
    a = [0] * N
    a[2^i] = 1
    a = K(tuple(a), ZZ(1))
    res.append(a)
  return res

# variant without global dictionary
@memoized # Warning: removing memoized will cause errors in evaluate_mod because the sqrt's are defined up to sign.
def field_sq_basis_local_mod(d, p):
  F = GF(p)
  F2 = F.extension(2)
  d_sq = [F2(di).sqrt() for di in d]
  I_sq = [1]
  for di in d_sq:
    I_sq += [di*i for i in I_sq]
  return I_sq

# Since sqrt is defined up to sign, we have to fix one variant for correctness.
# In SQROOTS we store fixed square roots globally for all fields.
SQROOTS = {}
# This procedure fix sqrt consistently with multiplication of the roots.
# We ensure that sqrt(a)*sqrt(b) == sqrt(a*b) over F_p^2
@memoized
def fixed_sqrt(a, p):
  assert a != -1, "not implemented"
  global SQROOTS
  if not (p in SQROOTS):
    SQROOTS[p] = {}

  if a in SQROOTS[p]:
    #print(f"fixed_sqrt(a = {a}, p = {p}) => {SQROOTS[p][a]} [stored]")
    return SQROOTS[p][a]

  F = GF(p)
  F2 = F.extension(2)

  if F2(a) == 0:
    return F2(0)
  
  #if F2(a) == 1:
  #  return F2(1)

  for di in SQROOTS[p]:
    if abs(a) == abs(di):
      continue
    if Mod(a, di) == 0:
      SQROOTS[p][a] = fixed_sqrt(di, p) * fixed_sqrt(a / di, p)
      return SQROOTS[p][a]
    if Mod(di, a) == 0:
     assert abs(di).valuation(abs(a)) == 1, f"not implemented, p = {p}, di = {di}, a = {a}, SQROOTS[p] = {SQROOTS[p]}"
     if di / a == a:
       SQROOTS[p][a] = F2(a).sqrt()
     elif di / a in SQROOTS[p]:
       assert SQROOTS[p][di/a] != 0, f"p = {p}, di = {di}, a = {a}, SQROOTS[p] = {SQROOTS[p]}"
       SQROOTS[p][a] = SQROOTS[p][di] / SQROOTS[p][di/a]
     #elif a < di / a:
     #  SQROOTS[p][a] = F2(a).sqrt()
     #  SQROOTS[p][di/a] = SQROOTS[p][di] / SQROOTS[p][a]
     #else:
     #  SQROOTS[p][a] = SQROOTS[p][di] / fixed_sqrt(di / a, p)
     else:
       SQROOTS[p][a] = F2(a).sqrt()
       SQROOTS[p][di/a] = SQROOTS[p][di] / SQROOTS[p][a]
     return SQROOTS[p][a]

  SQROOTS[p][a] = F2(a).sqrt()
  return SQROOTS[p][a]

@memoized # Warning: removing memoized will cause errors in evaluate_mod/evaluate_mod_ext, because the sqrt's are defined up to sign.
def field_sq_basis_mod(d, p):
  d_sq = []
  for di in d:
    d_sq.append(fixed_sqrt(di, p))
  I_sq = [1]
  for di in d_sq:
    I_sq += [di*i for i in I_sq] #TODO: fix all products (add them to SQROOTS) and require calling sqrt fix on top field before descending to subfields.
  return I_sq

@memoized
def mquad_disc(d):
    r""" Computes the discriminant of the multiquadratic number field $L = Q[sqrt(d_1), ..., sqrt(d_n)]$.
    
    Ref: Schmal B. Diskriminanten, ℤ-Ganzheitsbasen und relative Ganzheitsbasen bei multiquadratischen Zahlkörpern // Archiv der Mathematik. 1989. vol. 52, no. 3, p. 245-257.
    """
    n = len(d)
    m4 = True
    m3 = False
    m2 = False
    P = []
    res = 1
    for di in d:
        if Mod(di, 4) != 1:
            m4 = False
        if Mod(di, 4) == 3:
            m3 = True
        for p,e in factor(di):
            if p == -1 or p in P:
                continue
            if p == 2:
                m2 = True
            res = res * p
    if m4:
        r = 0
    elif ((not m2) and m3) or (m2 and m4):
        r = 2
    else:
        r = 3
    res = (2^r * res)^(2^(n-1))
    return res

@memoized
def mquad_idx(d): # Computes [O_K:Z[theta_K]] for K = Q(theta_K)
  K = field(d)
  disc_pol = K.absolute_polynomial().to_sage().change_ring(QQ).discriminant()
  disc_K = K.discriminant()
  idx_K = disc_pol / disc_K
  return idx_K

@memoized
def field_sage(d, food = None, names = "a"):
  if food == None:
    food = trees.get_food()

  if food != None:
    names = list(sorted(food.keys()))[:len(d)]
  K = NumberField([x^2 - d[i] for i in range(len(d))], names=names)
  return K

def from_sage_internal(d, a, i=0):
  K = field(d)
  if a.parent() == QQ:
    return K.from_QQ(a)
  if a.parent() == ZZ:
    return K.from_ZZ(a)

  #print(f"i = {i}, {a.base_ring().gens()[0]^2}, {d[i]}, {d}, {a.base_ring().gens()}")
  assert a.parent().gens()[0]^2 == d[i], f"{a.parent().gens()[0]^2} != {d[i]}"
  g = K.gens()[i]

  a0, a1 = a.list()
  r = from_sage_internal(d, a0, i+1) + from_sage_internal(d, a1, i+1) * g
  return r
