# Class for operations on ideal products of the form h prod_i P_i where h is an element of the field
# in power product representation and P_i are ideals from the factor base used in class group computations.

import field
import field_compact
import trees
import clgp
import fb
import sunits
import relations
import ideal_sqrt_cyc

from memoized import memoized

@memoized
def idealprods(d, d_parent = (), food = None, FB = None):
  n = len(d)
  N = 2^n
  KC = field_compact.field_compact(d)
  if food == None:
    food = trees.get_food()
  if FB == None:
    FB = clgp.get_FB(d, d_parent = d_parent, food = food)
  class IP:
    gens = d
    def __init__(f,*args):
      if len(args) == 1:
        g = args[0]
        if g.__class__ == KC:
          f.element = g
          f.powers = tuple([0]*len(FB))
          return
        if g.__class__ == KC:
          f.element = KC(g)
          f.powers = tuple([0]*len(FB))
          return
        if g.__class__ == fb.factorsb(d):
          f.element = KC.one()
          i = FB.index(g)
          f.powers = tuple([0] * i + [1] + [0] * (len(FB)-i-1))
          return
        # g is ideal product from subfield
        if (n >= 1 and g.__class__.gens == d[:-1]) or (n >= 2 and g.__class__.gens == d[:-2] + d[-1:]) or (n >= 2 and g.__class__.gens == d[:-2] + (d[-2]*d[-1],)):
          f.element = KC(g.element)
          f.powers = clgp.lift_e(g.powers, FB)
      
      assert len(args) == 2, "Wrong number of arguments!"
      
      g = args[0]
      
      if (n >= 1 and g.__class__.gens == d[:-1]) or (n >= 2 and g.__class__.gens == d[:-2] + d[-1:]) or (n >= 2 and g.__class__.gens == d[:-2] + (d[-2]*d[-1],)):
        f.element = KC(g)
        f.powers = clgp.lift_e(args[1], FB)
        return
      assert g.__class__.gens == d, f"Don't know how to embed elements and ideals from the field({g.__class__.gens}) to field({d})"

      f.element,f.powers = args # XXX: trusting caller to provide suitable values
      assert len(f.powers) == len(FB), "Wrong size of powers!"
      assert f.element.__class__ == KC, "Wrong type of element!"
      return
    def __repr__(f):
      return 'idealprod.idealprods(%s)(%s,%s)' % (d,f.element,f.powers)
    def __eq__(f,g):
      if g.__class__ != f.__class__: raise Exception('do not know how to compare %s to %s' % (f,g))
      if f.element.absnorm() / g.element.absnorm() not in [1,-1]:
        return False
      return f.powers == g.powers
    def __mul__(f,g):
      if f.__class__ != g.__class__:
        g = IP(g)
      helement = f.element * g.element
      hpowers = tuple(f.powers[i] + g.powers[i] for i in range(len(FB)))
      return IP(helement, hpowers)
    def __truediv__(f,g):
      if f.__class__ != g.__class__:
        g = IP(g)
      helement = f.element / g.element
      hpowers = tuple(f.powers[i] - g.powers[i] for i in range(len(FB)))
      return IP(helement, hpowers)
    def __pow__(f,e):
      if e == 0: return IP(KC.one(), [0]*len(FB))
      if e == 1: return f
      helement = f.element ^ e
      hpowers = tuple(vector(f.powers) * e)
      return IP(helement, hpowers)
    def sqrt(f):
      helement = f.element.sqrt(ideal_sqrt=True)
      hpowers = tuple(f.powers[i] / 2 for i in range(len(f.powers)))
      return IP(helement, hpowers)
    def absnorm(f):
      r = f.element.absnorm()
      for i in range(len(FB)):
        r *= FB[i].absolute_norm() ^ f.powers[i]
      return r
    def valuation(f, P):
      v = f.element.valuation(P)
      try:
        v += f.powers[FB.index(P)]
      except ValueError:
        pass
      return v
    def evaluate(f):
      # Evaluate product as Sage's ideal. Slow, for testing purposes only.
      K = KC.base_field().sage(food)
      r = prod([FB[i].to_sage(food)^f.powers[i] for i in range(len(f.powers))])
      return r * f.element.evaluate(ideal_sqrt=True).to_sage(K)
    def to_sage(f, K = None, names='a', quotient=False):
      helement = f.element.to_sage(K, names, quotient)
      return (helement, f.powers)

    def pseudo_gcd(f, g):
      # Computing sum of ideals. Result is correct only if f.element and g.element are represented using the same basis.
      helem = f.element.pseudo_gcd(g.element)
      hpowers = tuple(min(f.powers[i], g.powers[i]) for i in range(len(FB)))
      return IP(helem, hpowers)

    def generator(f, SU):
      '''
      Return generator of ideal product when this product is principal.
      
      Requires S-unit group SU.
      '''
      M = SU.relations_matrix()
      try:
        v = M.solve_left(vector(f.powers))
      except Exception as err:
        raise Exception(f"Ideal product is not principal! Error: {err}")
      s = prod([SU.gens()[i]^v[i] for i in range(len(v))], KC.one())
      return s * f.element
    
    def normalize_valuations(f, J, fb = None):
      '''
      Changes exponent vector b in I = alpha prod_i P_i^(b_i) in such way that v_(P_i)(J) = v_(P_i)(I) for each P_i in the factor base. 
      '''
      
      if fb == None:
        fb = FB
      
      assert fb != None

      v1 = vector([f.valuation(fb[j]) for j in range(len(fb))])
      v2 = vector([J.valuation(fb[j]) for j in range(len(fb))])
      v = v2 - v1
      if not v.is_zero():
          print("non-zero vector: {v}")
      hpowers = tuple(vector(f.powers) + v)

      return IP(f.element, hpowers)

    def reduce_mod_units(f, rollback=True):
      '''
      Given an ideal I = alpha * prod_i FB_i^(b_i), where alpha = prod alpha_j^(a_j), the method
      reduces elements alpha_1, alpha_2, ... modulo units.
      '''

      el = f.element.reduce_mod_units(rollback=rollback)
      I = IP(el, tuple(f.powers))
      return I

    def reduce_mod_sunits_base_exp(f, SU):
      A = SU.relations_matrix()
      FB0 = SU.factor_base(sage=False)
      FB = SU.factor_base(sage=True)

      hrem = KC.one()
      helement = KC.one()
      hpowers = vector(f.powers)
      for i in range(len(f.element.elements)):
        t = walltime()
        v = [f.element.powers[i] * f.element.elements[i].valuation(FB[j]) for j in range(len(FB))]
        print(f"[base_exp_vals_{i}|{walltime(t)} sec.]", flush=True)

        if fb.is_smooth(f.element.elements[i].absnorm(), FB0):
          hpowers += vector(v)
        else:
          t = walltime()
          w = clgp.approx_CVP(d, list(v))
          assert len(w) == len(SU.gens())
          print(f"[base_exp_acvp_{i}|{walltime(t)} sec.]", flush=True)
          
          t = walltime()
          s = prod([SU.gens()[i]^w[i] for i in range(len(w))], KC.one()).compute_roots(ideal_sqrt=True)
          print(f"[base_exp_s_prod_{i}|{walltime(t)} sec.]", flush=True)
          
          #el = KC(f.element.elements[i])^f.element.powers[i] / s
          #el = KC(el.reduce_mod_units(rollback=True).evaluate().reduce_mod_units(rollback=True))
          #el = KC(el.evaluate())
          
          t = walltime()
          el,r,su = (KC(f.element.elements[i])^f.element.powers[i]).rdiv(s, log_sunits=True, FB=FB)
          print(f"[base_exp_rdiv_{i}|{walltime(t)} sec.]", flush=True)

          helement *= el
          hrem *= r
          hpowers += vector(w) * A

      if hrem != KC.one() and (not vector(hrem.powers).is_zero()):
        t = walltime()
        r_v = vector([hrem.valuation(FB[j]) for j in range(len(FB))]) # TODO: slow, make separate class for S-units with powers
        hpowers -= vector(r_v)
        print(f"[base_exp_remainder(len={len(hrem.elements)})|{walltime(t)} sec.]", flush=True)

      return IP(helement, tuple(hpowers))

    def reduce_mod_sunits(f, SU):
      # A = SU.relations_matrix()
      FB = SU.factor_base(sage=True)

      hpowers = vector(f.powers)

      helement,su = f.element.reduce_mod_sunits(SU)

      if su != SU.one():
        v = vector([su.valuation(FB[j]) for j in range(len(FB))]) # TODO: slow, make separate class for S-units with powers
        hpowers += v
      
      return IP(helement, tuple(hpowers))
    
    def sqrt(f, SU):
      return ideal_sqrt_cyc.ideal_sqrt(d, f.element, f.powers, SU)
    
    def sqrt_rat(f, SU):
      return ideal_sqrt_cyc.ideal_sqrt_rat(d, f.element, f.powers, SU)

    def reduce_lin_dep(f, mod_sunits=False, SU=None, rollback=True):
      if mod_sunits:
        FB = SU.factor_base(sage=True)
        helement,su = f.element.reduce_lin_dep(mod_units=True, mod_sunits=mod_sunits, SU=SU, rollback=rollback)
      else:
        helement = f.element.reduce_lin_dep(mod_units=True, mod_sunits=mod_sunits, SU=SU, rollback=rollback)

      hpowers = vector(f.powers)
      
      if mod_sunits and su != SU.one():
        v = vector([su.valuation(FB[j]) for j in range(len(FB))]) # TODO: slow, make separate class for S-units with powers
        hpowers += v

      return IP(helement, tuple(hpowers))

    @staticmethod
    def random(bits = 8, exp_bound = 5, elems = 4, exp_neg = False):
      felement = KC.random(bits, exp_bound, elems, exp_neg)
      if exp_neg:
        fpowers = tuple(ZZ.random_element(-exp_bound,exp_bound) for j in range(len(FB)))
      else:
        fpowers = tuple(ZZ.random_element(1,exp_bound) for j in range(len(FB)))
      return IP(felement, fpowers)
    @staticmethod
    def base_field():
      return KC.base_field()
    @staticmethod
    def one():
      return IP(KC.one(), [0]*len(FB))
    @staticmethod
    def zero():
      return IP(KC.zero(), [0]*len(FB))

  return IP
