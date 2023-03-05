# Arithmetic using representation of field elements as power products.

import field
import ideal
from memoized import memoized

@memoized
def field_compact(d):
  n = len(d)
  N = 2^n
  K = field.field(d)
  class KC:
    gens = d
    def __init__(f,*args):
      if len(args) == 1:
        g = args[0]
        if g.__class__ == K:
          f.elements = (g,)
          f.powers = (1,)
          return
        if n >= 1 and  g.__class__.gens == field_compact(d[:-1]).gens:
          f.elements = tuple(K(g.elements[i]) for i in range(len(g.elements))) 
          f.powers = tuple(g.powers)
          return
        if n >= 2 and g.__class__.gens == field_compact(d[:-2] + d[-1:]).gens:
          f.elements = tuple(K(g.elements[i]) for i in range(len(g.elements)))
          f.powers = tuple(g.powers)
          return
        if n >= 2 and g.__class__.gens == field_compact(d[:-2] + (d[-2]*d[-1],)).gens:
          f.elements = tuple(K(g.elements[i]) for i in range(len(g.elements)))
          f.powers = tuple(g.powers)
          return
      if len(args) == 2: # XXX: trusting caller to provide suitable values
        f.elements,f.powers = args
        return
      raise Exception('not known how to initialize field_compact(%s)(%s)' % (str(d),args))
    def __repr__(f):
      return 'field_compact.field_compact(%s)(%s,%s)' % (d,f.elements,f.powers)
    def __eq__(f,g):
      if g.__class__ != f.__class__: raise Exception('do not know how to compare %s to %s' % (f,g))
      f.compactify()
      g.compactify()
      if len(f.elements) != len(g.elements): # We assume that compactify() joined duplicate entries.
        return False
      if f.absnorm() != g.absnorm():
        return False
      for i in range(len(f.elements)):
        if i >= len(g.elements):
          return False
        if f.elements[i] != g.elements[i]:
          try:
            j = g.elements.index(f.elements[i])
            if g.powers[j] != f.powers[i]:
              return False
          except ValueError:
            return False
        elif f.powers[i] != g.powers[i]:
          return False
      return True
    def __mul__(f,g):
      if f.__class__ != g.__class__:
        g = KC(g)
      helements = list(f.elements)
      hpowers = list(f.powers)
      for i in range(len(g.elements)):
        if i < len(f.elements) and g.elements[i] == f.elements[i]:
          hpowers[i] += g.powers[i]
        else:
          try:
            j = f.elements.index(g.elements[i])
            hpowers[j] += g.powers[i]
          except ValueError:
            helements.append(g.elements[i])
            hpowers.append(g.powers[i])
      h = KC(tuple(helements),tuple(hpowers))
      return h.trim()
    def __truediv__(f,g):
      g_inv_powers = [-g.powers[i] for i in range(len(g.powers))]
      g_inv = KC(g.elements, tuple(g_inv_powers))
      h = f * g_inv
      # TODO: make reduction/shortening
      return h
    def __pow__(f,e):
      if e == 0: return KC(K.one())
      if e == 1: return f
      g = vector(f.powers) * e
      return KC(tuple(f.elements), tuple(g))
    def sqrt(f, ideal_sqrt=False):
      helements = [0]*len(f.elements)
      hpowers = [0]*len(f.elements)
      el = K.one()
      for i in range(len(f.elements)):
        e1,e0 = divmod(f.powers[i], 2) # f.powers[i] = e0 + 2*e1
        hpowers[i] = e1
        helements[i] = f.elements[i]
        if e0 == 1:
          el *= f.elements[i]
      h = KC(tuple(helements),tuple(hpowers))
      if ideal_sqrt:
        h = KC(ideal.idealsqrtshorten(len(d), 2^len(d), d, el)) * h
      else:
        h = KC(el.sqrt()) * h
      return h.trim()
    def conj(f): # conjugate last entry of d
      helements = [f.elements[i].conj() for i in range(len(f.elements))]
      hpowers = f.powers
      return KC(tuple(helements),tuple(hpowers))
    def symbols2_log(f,low,high):
      if len(f.elements) == 0:
        return 0
      r = vector(GF(2), f.elements[0].symbols2_log(low,high)) * f.powers[0]
      for i in range(1,len(f.elements)):
        r += vector(GF(2), f.elements[i].symbols2_log(low,high)) *  f.powers[i]
      return vector(ZZ, r)
    def symbols_log(f, low, high):
      if len(f.elements) == 0:
        return 0
      r = vector(GF(2), f.elements[0].symbols_log(low,high)) * f.powers[0]
      for i in range(1,len(f.elements)):
        r += vector(GF(2), f.elements[i].symbols_log(low,high)) *  f.powers[i]
      return vector(ZZ, r)
    def absnorm(f):
      r = 1
      s_r = 1
      for i in range(len(f.elements)):
        if f.powers[i].denominator() == 2:
          s_r *= f.elements[i].absnorm() ^ f.powers[i].numerator()
        else:
          r *= f.elements[i].absnorm() ^ f.powers[i]
      return r * s_r ^ (1/2)
    def valuation(f, P):
      v = 0
      for el,e in zip(f.elements, f.powers):
        v += el.valuation(P) * e
      return v
    def factor(f, S): # S is a list of prime ideals
      h = KC(K.one())
      r = vector([0]*len(S))
      for el,e in zip(f.elements, f.powers):
        nrm = el.absnorm()
        el_v = vector([0]*len(S))
        for i in range(len(S)):
          el_v[i] = el.valuation(S[i]) * e
          nrm /= S[i].absolute_norm() ^ el_v[i]
        if abs(nrm) != 1:
          h *= el
        else:
          r += el_v
      return (h,tuple(r))
    def evaluate(f, ideal_sqrt=True): # evaluates product. Very slow, works only for small elements
      el = K.one()
      el_r = K.one()
      for i in range(len(f.elements)):
        if f.powers[i].denominator() == 1:
          el *= f.elements[i] ^ ZZ(f.powers[i])
        elif f.powers[i].denominator() == 2:
          el_r *= f.elements[i] ^ f.powers[i].numerator()
        else:
          raise NotImplementedError(f"Extraction of {f.powers[i].denominator()}-roots is not implemented!")
      if not ideal_sqrt:
        el_r = el_r.sqrt()
      else:
        el_r = ideal.idealsqrtshorten(len(d), 2^len(d), d, el_r)
      return el * el_r
    def compute_roots(f, ideal_sqrt=False): # expand rational exponents
      el = KC.one()
      el_r = KC.one()
      for i in range(len(f.elements)):
        if f.powers[i].denominator() == 1:
          el *= KC(f.elements[i]) ^ ZZ(f.powers[i])
        elif f.powers[i].denominator() == 2:
          if ideal_sqrt:
            try:
              rt = ideal.idealsqrtshorten(len(d), 2^len(d), d, f.elements[i])
              el *= KC(rt) ^ f.powers[i].numerator()
            except:
              el_r *= KC(f.elements[i]) ^ f.powers[i].numerator()
          else:
            try:
              rt = f.elements[i].sqrt()
              el *= KC(rt) ^ f.powers[i].numerator()
            except:
              el_r *= KC(f.elements[i]) ^ f.powers[i].numerator()
        else:
          raise NotImplementedError(f"Extraction of {f.powers[i].denominator()}-roots is not implemented!")
      el_r = el_r.sqrt(ideal_sqrt=ideal_sqrt)
      return el * el_r
    def to_sage(f, K = None, names='a', quotient=False):
      helements = tuple(f.elements[i].to_sage(K, names, quotient) for i in range(len(f.elements)))
      return (helements, f.powers)
    def trim(f, inplace=False):
      helements = []
      hpowers = []
      for i in range(len(f.elements)):
        if f.elements[i] != K.one() and f.powers[i] != 0:
          helements.append(f.elements[i])
          hpowers.append(f.powers[i])
      if inplace:
        f.elements = helements
        f.powers = hpowers
      if len(helements) == 0:
        return KC.one()
      return KC(helements, hpowers)
    def compactify(f):
      # Joining of duplicate elements.
      # TODO: make flag mark processed elements
      for i in range(len(f.elements)):
        if f.elements[i] == K.one():
          continue
        for j in range(i+1, len(f.elements)):
          if f.elements[i] == f.elements[j]:
            f.powers[i] += f.powers[j]
            f.elements[j] = K.one()
      f.trim(inplace=True)
      return f
    def pseudo_gcd(f, g):
      f.compactify()
      g.compactify()
      r = KC.one()
      for i in range(len(f.elements)):
        if f.powers[i] == 0:
          continue
        for j in range(len(g.elements)):
          if g.powers[j] == 0:
            continue
          if f.elements[i] == g.elements[j]:
          #if abs(f.elements[i].absnorm() / g.elements[j].absnorm()) == 1:
            r *= KC(f.elements[i]) ^ min(f.powers[i], g.powers[j])
            continue
      return r.trim()
    @staticmethod
    def random(bits = 8, exp_bound = 5, elems = 4, exp_neg = False):
      felements = tuple(K.random(bits) for i in range(elems))
      if exp_neg:
        fpowers = tuple(ZZ.random_element(-exp_bound,exp_bound) for j in range(elems))
      else:
        fpowers = tuple(ZZ.random_element(1,exp_bound) for j in range(elems))
      return KC(felements, fpowers)
    @staticmethod
    def base_field():
      return K
    @staticmethod
    def one():
      return KC(K.one())
    @staticmethod
    def zero():
      return KC(K.zero())

  return KC
