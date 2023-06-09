# Arithmetic using representation of field elements by power products.

import field
import ideal
import fb
import units
import char
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
      f.compactify(inplace=True)
      g.compactify(inplace=True)
      if len(f.elements) != len(g.elements): # We assume that compactify() joined duplicate entries.
        return False
      #if f.absnorm() != g.absnorm(): # FIXME: this can be expensive, make compact representation of norms and uncomment
      #  return False
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
        assert type(f.powers[i]) == Integer or (type(f.powers[i]) == Rational and f.powers[i].denominator()==1), f"Powers of type {f.powers[i].__class__} are not supported!" # call compute square roots first.
        e1,e0 = divmod(ZZ(f.powers[i]), 2) # f.powers[i] = e0 + 2*e1
        hpowers[i] = e1
        helements[i] = f.elements[i]
        if e0 == 1:
          el *= f.elements[i]
      h = KC(tuple(helements),tuple(hpowers))
      if ideal_sqrt:
        h = KC(ideal.idealsqrtshorten(len(d), 2^len(d), d, el)) * h
      else:
        h = KC(el.sqrt()) * h
      return h.trim(inplace=True, mod_units=ideal_sqrt)
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
    def altsymbols_log(f, low, high):
      if len(f.elements) == 0:
        return 0
      K = field.field(d)
      f_s = []
      f_d = []
      qs = []
      for i in range(len(f.elements)):
        a = f.elements[i]
        assert a != K.zero()
        f_d += [a.denom]
        aq = a.get_qs(low,high)
        f_s += [[aq[j][0] for j in range(len(aq))]]
        if qs == []:
          qs = [aq[j][1] for j in range(len(aq))]
      
      aq = []
      aqd = []
      vec = []
      for i in range(len(f.powers)):
        if (f.powers[i] != 0):
           vec += [f.powers[i]]
           aq += [f_s[i]]
           aqd += [f_d[i]]
      vec = vector(vec)
      denom = vec.denominator()

      s = char.altfinal2(aq, aqd, qs, denom*vec, denom)
      if 0 in s:
        raise Exception('zero symbol')
      return tuple(0 if si == 1 else 1 for si in s)

    def bit_size(f):
      r = 0
      for i in range(len(f.elements)):
        r += f.elements[i].bit_size() + ZZ(f.powers[i]).nbits()
      return r

    def absnorm(f, known_factors = []):
      #print(f"\nf.elements.absnorm() = {[f.elements[i].absnorm().factor() for i in range(len(f.elements))]}")
      #print(f"powers = {f.powers}")
      nrm = 1
      r = 1
      denom = 1
      l = 1
      for i in range(len(f.elements)):
        if f.powers[i].denominator() == 1:
          nrm *= f.elements[i].absnorm() ^ ZZ(f.powers[i])
        else:
          l = lcm(denom, ZZ(f.powers[i].denominator()))
          r = r^ZZ(l/denom) * f.elements[i].absnorm() ^ ZZ(f.powers[i].numerator() * l / f.powers[i].denominator())
          denom = l
      # reduce part with non-integer powers using known factors of norm
      for i in range(len(known_factors)):
        v = r.valuation(known_factors[i])
        if v == 0:
          continue
        assert Mod(v, denom) == 0
        nrm *= ZZ(known_factors[i]) ^ ZZ(v / denom)
        r = r / ZZ(known_factors[i]) ^ ZZ(v)

      res = QQ(nrm) * QQ(r^(1/denom))
      return res

    def norm_size(f):
      return abs(f.absnorm().numerator()) + abs(f.absnorm(denominator()))

    def lognorm(f, p=2, prec=None):
      if len(f.elements) == 0:
        return 0
      if prec == None:
        prec = units.logbits
      v = vector(f.elements[0].approxlog(prec=prec)) * f.powers[0]
      for i in range(1, len(f.elements)):
        v += f.powers[i] * vector(f.elements[i].approxlog(prec=prec))
      return v.norm(p=p)

    def valuation(f, P):
      v = 0
      for el,e in zip(f.elements, f.powers):
        if e == 0:
          continue
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

    def evaluate(f, ideal_sqrt=False, mod_units=False, rollback=True):
      # evaluates product. Very slow, it works only for small elements
      el = K.one()
      el_r = K.one()
      for i in range(len(f.elements)):
        if f.powers[i].denominator() == 1:
          el *= f.elements[i] ^ ZZ(f.powers[i])
        elif f.powers[i].denominator() == 2:
          el_r *= f.elements[i] ^ f.powers[i].numerator()
        else:
          raise NotImplementedError(f"Extraction of {f.powers[i].denominator()}-roots is not implemented!")
        if mod_units:
          el = el.reduce_mod_units(rollback=rollback)
          el_r = el_r.reduce_mod_units(rollback=rollback)
      if not ideal_sqrt:
        el_r = el_r.sqrt()
      else:
        el_r = ideal.idealsqrtshorten(len(d), 2^len(d), d, el_r)
      return el * el_r

    def compute_roots(f, ideal_sqrt=False): # expand rational exponents
      el = KC.one()
      el_r = KC.one()
      denom = 1
      for i in range(len(f.elements)):
        if f.powers[i].denominator() == 1:
          el *= KC(f.elements[i]) ^ ZZ(f.powers[i])
        elif Mod(f.powers[i].denominator(), 2) == 0:
          found = False
          if f.powers[i].denominator() == 2:
            if ideal_sqrt:
              try:
                rt = ideal.idealsqrtshorten(len(d), 2^len(d), d, f.elements[i])
                el *= KC(rt) ^ f.powers[i].numerator()
                found = True
              except:
                pass
            else:
              try:
                rt = f.elements[i].sqrt()
                el *= KC(rt) ^ f.powers[i].numerator()
                found = True
              except:
                pass
          if not found:
            l = lcm(denom, f.powers[i].denominator())
            el_r = el_r ^ (l / denom) * KC(f.elements[i]) ^ (f.powers[i].numerator() * l / f.powers[i].denominator())
            denom = l
        else:
          raise NotImplementedError(f"Extraction of {f.powers[i].denominator()}-roots is not implemented!")
      r = el
      if denom != 1:
        while Mod(denom, 2) == 0:
          el_r = el_r.sqrt(ideal_sqrt=ideal_sqrt)
          denom /= 2
      assert denom == 1, f"Extraction of {denom}-roots is not implemented!"
      r = r * el_r
      return r

    def to_sage(f, K = None, names='a', quotient=False):
      helements = tuple(f.elements[i].to_sage(K, names, quotient) for i in range(len(f.elements)))
      return (helements, f.powers)

    def trim(f, inplace=False, mod_units=False):
      helements = []
      hpowers = []
      for i in range(len(f.elements)):
        if f.elements[i] == K.one() or f.powers[i] == 0:
          continue
        if mod_units and f.elements[i].absnorm() in [1,-1]:
          ff = f.elements[i].reduce_mod_units(rollback=True)
          if ff not in [K.one(), -K.one()]: # TODO: check also full torsion for imaginary fields
            helements.append(ff)
            hpowers.append(f.powers[i])
          continue
        helements.append(f.elements[i])
        hpowers.append(f.powers[i])
      if len(helements) == 0:
        return KC.one()
      if inplace:
        f.elements = helements
        f.powers = hpowers
        return f
      return KC(tuple(helements), tuple(hpowers))

    def compactify(f, inplace=False, mod_units=False):
      '''
      Joins duplicate elements. If mod_units = True it also try to join elements that differ by a unit.
      '''

      helements = list(f.elements)
      hpowers = list(f.powers)
  
      # TODO: make flag to mark already processed elements
      for i in range(len(f.elements)):
        if helements[i] == K.one():
          continue
        for j in range(i+1, len(helements)):
          if helements[i] == helements[j]:
            hpowers[i] += hpowers[j]
            helements[j] = K.one()
          elif mod_units and (helements[j].absnorm() / helements[i].absnorm() in [1,-1]):
            fji = (helements[j] / helements[i]).reduce_mod_units(rollback=True)
            # We can't just exclude element fji since there are elements h of K s.t. 
            # N(h) = +-1, h is not integral, and h is not a unit in O_K.
            if fji.bit_size() < helements[j].bit_size():
              hpowers[i] += hpowers[j]  
              helements[j] = fji

      if inplace:
        f.elements = tuple(helements)
        f.powers = tuple(hpowers)
        h = f
      else:
        h = KC(tuple(helements), tuple(hpowers))

      return h.trim(inplace=True, mod_units=mod_units)

    def rdiv(f, g, bits=True, norms=False, log_units=False, log_sunits=False, l=2, prec=None, FB=None, mod_units=True):
      ''' 
      Given two elements f = prod f_i^a_i and g = prod g_i^b_i of the field divides elements f_i by g_j when 
      the size of f_i/g_j^sign(b_i) reduces with respect to given metric measuring a size of elements.
      
      The "size" is measured using the following metrics:
      1. ||ApproxLog(.)||_l, norm of logarithmic complex embeddings of the element (log_units = True).
      2. ||Log_S(.)||_l, norm of logarithmic S-embeddings of the element. Requires factor base (FB) to compute valuations (log_sunits = True, FB != None).
      3. bit_size(.), bit size of element's internal representation (bits = True).
      4. N(.), absolute norm of the element (norms=True).
      
      If mod_units = True the method reduces all intermediate results modulo units.
      
      Returns triple (f/s, g/s, s) where s s.t. size(f/s) <= size(f) and s divides (exactly) g.
      '''

      assert bits or norms or log_units or (log_sunits and FB!=None), "enable at least one metric"

      h = KC.one()
      s = KC.one()
      gpowers = list(g.powers)

      for i in range(len(f.elements)):
        el = f.elements[i]

        if vector(gpowers).is_zero():
          break

        j = 0
        while j < len(g.elements):
          el_old = el
          if gpowers[j] > 0:
            if log_units:
              l1 = vector(el.approxlog(prec=prec))
              l2 = vector(g.elements[j].approxlog(prec=prec))
              if RR((l1-l2).norm(l)) >= RR(l1.norm(l)):
                j += 1
                continue
            if norms:
              n1 = el.absnorm()
              n2 = g.absnorm()
              if n1/n2 >= n1:
                j += 1
                continue
            if log_sunits:
              v1 = vector([el.valuation(FB[k]) for k in range(len(FB))])
              v2 = vector([g.elements[j].valuation(FB[k]) for k in range(len(FB))])
              if (v1-v2).norm(l) >= v1.norm(l):
                j += 1
                continue
            el = el / g.elements[j]
            if mod_units:
              el = el.reduce_mod_units(rollback=True)
            if bits and (el.bit_size() >= el_old.bit_size()):
              el = el_old
              j += 1
              continue
            gpowers[j] -= f.powers[i]
            s *= KC(g.elements[j]) ^ f.powers[i] 
            j = 0
          elif gpowers[j] < 0:
            if log_units:
              l1 = vector(el.approxlog(prec=prec))
              l2 = vector(g.elements[j].approxlog(prec=prec))
              if RR((l1+l2).norm(l)) >= RR(l1.norm(l)):
                j += 1
                continue
            if norms:
              n1 = el.absnorm()
              n2 = g.absnorm()
              if n1*n2 >= n1:
                j += 1
                continue
            if log_sunits:
              v1 = vector([el.valuation(FB[k]) for k in range(len(FB))])
              v2 = vector([g.elements[j].valuation(FB[k]) for k in range(len(FB))])
              if RR((v1+v2).norm(l)) >= RR(v1.norm(l)):
                j += 1
                continue
            el = el * g.elements[j]
            if mod_units:
              el = el.reduce_mod_units(rollback=True)
            if bits and (el.bit_size() >= el_old.bit_size()):
              el = el_old
              j += 1
              continue
            s /= KC(g.elements[j]) ^ f.powers[i]
            gpowers[j] += f.powers[i]
            j = 0
          else:
            j += 1

        h *= KC(el)^ZZ(f.powers[i])

      h = h.trim(mod_units=mod_units)
      r = KC(g.elements, tuple(gpowers)).trim(mod_units=mod_units)
      return (h, r, s)

    def pseudo_gcd(f, g):
      f.compactify(inplace=True)
      g.compactify(inplace=True)
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

    def reduce_mod_units(f, rollback=True):
      # TODO: try to reduce all elements at once, not just per element
      helements = []
      hpowers = []
      for i in range(len(f.elements)):
#        if f.elements[i].absnorm() in [1,-1]: # doesn't work, since we need to check that element is integral
#          continue
        el = f.elements[i].reduce_mod_units(rollback=rollback)
        if rollback and (el.bit_size() >= f.elements[i].bit_size()):
          helements.append(f.elements[i])
          hpowers.append(f.powers[i])
        else:
          helements.append(el)
          hpowers.append(f.powers[i])
      return KC(tuple(helements), tuple(hpowers)).compactify(inplace=True, mod_units=True)

    def reduce_lin_dep(f, logbits=None, mod_units=False, mod_sunits=False, SU=None, rollback=False):

      if mod_sunits:
        assert SU != None

      S = f.elements
      if logbits==None:
        logbits = units.logbits*2
      
      t = walltime()
      
      M = matrix(QQ,[[QQ(i==j) for i in range(len(S))] + [QQ(v) for v in S[j].approxlog(prec=logbits)] for j in range(len(S))])
      
      M = M.LLL()

      H = [Mi[len(S):] for Mi in M]

      # obtaing transformation matrix
      T = [Mi[:len(S)] for Mi in M]
      T = matrix(ZZ, T)

      print(f"T (transform. matrix):\n{T}")

      MM = matrix(QQ,[[QQ(v) for v in S[j].approxlog(prec=logbits)] for j in range(len(S))])
      assert T * MM == matrix(H), f"M*MM = {T * MM}, H = {H}"

      v = vector(f.powers)

      T_inv = T^(-1)
      print(f"T^(-1) = {T_inv}")

      powers = tuple(v * T_inv)

      assert len(powers) == T.nrows()

      if mod_sunits:
        s = SU.one()
      
      r = KC.one()

      for i in range(T.nrows()):
        if powers[i] == 0:
          continue
        #if H[i].is_zero():
        #  continue
        if mod_sunits:
          s_el = SU.one()
        el = K.one()
        for j in range(T.ncols()):
          if mod_sunits:
            ex,s0 = S[j].pow_mod_sunits(T[i,j], SU)
            s_el *= s0
            el,s0 = el.mul_mod_sunits(ex, SU, mod_units=mod_units)
            s_el *= s0
          elif mod_units:
            el *= (S[j] ^ T[i,j]).reduce_mod_units(rollback=True)
          else:
            el *= S[j] ^ T[i,j]

        r *= KC(el) ^ powers[i]

        # We expect that all values in r.elements are independent, after LLL-reduction.
        # So, bit sizes can't decrease when we add new element el. Thus we perform early rollback.
        if rollback and r.bit_size() >= f.bit_size():
          print(f"[reduce_lin_dep(rollback)|{walltime(t)} sec.]")
          if mod_sunits:
            return (f, SU.one())
          return f
        
        if mod_sunits:
          s *= s_el ^ powers[i]
        if H[i].is_zero():
          #assert el == K.one(), f"{el.to_sage(K.sage())} != 1"
          print(f"H_i = 0, el = {el.to_sage(K.sage())}")

      print(f"[reduce_lin_dep(elements: {len(S)} => {len(r.elements)}, bits: {f.bit_size()} => {r.bit_size()}])|{walltime(t)} sec.]")

      if mod_sunits:
        return (r,s)
      return r

    def reduce_base_mod_sunits(f, SU, mod_units=True):
      '''
      Given element f = prod_i f_i^e_i of the field K reduces each component f_i modulo S-units.

      Returns pair (f/s, s) where s is an S-unit.
      '''
      h = KC.one()
      s = KC.one()
      for i in range(len(f.elements)):
        hi,si = f.elements[i].reduce_mod_sunits(SU, mod_units=mod_units)
        h *= KC(hi) ^ f.powers[i]
        s *= si ^ f.powers[i]
      return (h,s)

    def reduce_all_mod_sunits(f, SU, mod_units=True):
      '''
      Given element f = prod_i f_i^e_i of the field K reduces (mod S-units) all product at once.

      Returns pair (f/s, s) where s is an S-unit and size(f/s) <= size(s).
      Where size is the bit_size of elements.
      '''

      A = SU.relations_matrix()
      FB = SU.factor_base(sage=True)

      v = vector([f.valuation(FB[i]) for i in range(len(FB))])

      # f.absnorm() can be slow for non-reduced f.
      # if fb.is_smooth(f.absnorm(), FB0): 
      #   return (K.one(), SU.from_prime_prod(v))

      w = clgp.approx_CVP(d, list(v))
      assert len(w) == len(SU.gens())

      wA = vector(w) * A

      # reduction failed, resulting vector has larger norm
      if RR((vector(v)-wA).norm()) >= RR(vector(v).norm()):
        return (f, SU.one())
 
      s = prod([SU.gens()[i]^w[i] for i in range(len(w))], KC.one()).compute_roots(ideal_sqrt=mod_units)

      h,r,s = f.rdiv(s, bits=False, log_sunits=True, FB = FB)

      # TODO: reduce linear dependency using LLL, otherwise this is not very usefull procedure.

      return (h,s)

    def reduce_mod_sunits(f, SU, mod_units=True):
      '''
      Reduces element f = prod f_i^e_i modulo S-units in the following way:
      1) reduces base elements f_i individualy to remove S-unit components.
      2) reduces base elements f_j (j=1, ...) all together to remove dependence by S-unit between them.

      Returns pair (f/s, s) where s is an S-unit and size(f/s) <= size(s).
      Where size is the bit_size of elements.
      '''

      h,s = f.reduce_base_mod_sunits(SU, mod_units=mod_units)
      #h,s0 = h.reduce_all_mod_sunits(SU, mod_units=mod_units) # rdiv doesn't provide good reduction in this case
      #s *= s0

      # TODO: other reductions and verbose mode
      return (h,s)

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
