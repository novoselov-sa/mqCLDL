# This is a hack to run unittests from parent directory
import os
import sys
sys.path.append(os.path.dirname(__file__) + "/../")

import unittest
import field
import field_compact
import units
import random as rnd

class TestFieldCompactMethods(unittest.TestCase): 
  
  def setUp(self):
    self.fields = [
      (5, 13),
      (-11, -31),
      (-3, -7, -11),
      (-11, -31, -19),
      (-7, -11, -19, -23),
      (-11, -19, -31, -47),
      (-3, 5*17),
      (7*11, 23*43),
      (3*5, 7*13, 3*13)
    ]

  def test_eq(self):
    for d in self.fields:
      KC = field_compact.field_compact(d)
      K = field.field(d)

      a = KC.random()
      b = KC(list(a.elements).copy(), list(a.powers).copy())

      i = ZZ.random_element(0,len(a.elements))
      j = ZZ.random_element(0,len(a.elements))
      b.elements[i], b.elements[j] = b.elements[j], b.elements[i]
      b.powers[i], b.powers[j] = b.powers[j], b.powers[i]

      self.assertTrue(a == b)

      a = KC.random(elems=2)
      b = KC.random(elems=5)
      self.assertTrue(a*b == b*a)

  def test_mult(self):
    for d in self.fields:
      K1 = field_compact.field_compact(d)
      K2 = field.field(d)
      a = K2.random()
      b = K2.random()
      self.assertEqual(((K1(a)^2) * (K1(b)^3)).evaluate(), a^2 * b^3)

  def test_absnorm(self):
    for d in self.fields:
      K = field.field(d)
      KC = field_compact.field_compact(d)
      a = KC.random()
      self.assertEqual(a.absnorm(), a.evaluate().to_sage(K.sage()).absolute_norm())

  def test_valuation(self):
    for d in self.fields:
      K = field.field(d)
      KC = field_compact.field_compact(d)
      a = KC.random(bits = 2, exp_bound=2, elems=2)
      for p,e0 in a.absnorm().factor():
        for P,e1 in K.sage().factor(p):
          self.assertEqual(a.valuation(P), a.evaluate().to_sage(K.sage()).valuation(P))

  def test_factor(self):
    for d in self.fields:
      K = field.field(d)
      KC = field_compact.field_compact(d)

      a = KC.random(bits=2, exp_bound=2, elems=2)
      S = []
      for p in primes(20):
        for P,e in K.sage().factor(p):
          S.append(P)
      h,e = a.factor(S)
      aO = K.sage().ideal(a.evaluate().to_sage(K.sage()))
      hO = K.sage().ideal(h.evaluate().to_sage(K.sage()))
      I = prod([S[i] ^ e[i] for i in range(len(S))], hO)
      self.assertEqual(aO, I)

  def test_sqrt(self):
    for d in self.fields:
      K = field.field(d)
      KC = field_compact.field_compact(d)
      a = KC.random()
      s = a * a
      r = s.sqrt()
      self.assertIn(r.evaluate(), [a.evaluate(), -a.evaluate()])

      a = K.random()
      s = KC(a^2)
      self.assertIn(s.sqrt().evaluate(), [a, -a])

  def test_pow(self):
    for d in self.fields:
      K = field_compact.field_compact(d)
      K2 = field.field(d).sage()
      a = K.random()
      b = K.random()

      self.assertEqual((a^3 * b^3).evaluate(), ((a*b)^3).evaluate())
      self.assertEqual((a^5).evaluate().to_sage(), a.evaluate().to_sage()^5)

      c = (a^(-3)).to_sage(K2)
      c = prod([c[0][i]^c[1][i] for i in range(len(c[0]))])
      self.assertEqual(c, a.evaluate().to_sage(K2)^(-3))

  def test_symbols_log(self):
    for d in self.fields:
      KC = field_compact.field_compact(d)
      K = KC.base_field()
      r = 10

      a = KC.random()^2
      s = a.symbols_log(0,r)
      self.assertEqual(vector(s), vector([0]*r))

      a = KC(K.random()^2)
      s = a.symbols_log(0,r)
      self.assertEqual(vector(s), vector([0]*r))

      a = KC.random()
      b = KC.random()
      ab = a*b

      a_s = vector(GF(2), a.symbols_log(0,r))
      b_s = vector(GF(2), b.symbols_log(0,r))
      ab_s = vector(GF(2), ab.symbols_log(0,r))
      self.assertEqual(a_s + b_s, ab_s)

  def test_symbols2_log(self):
    for d in self.fields:
      KC = field_compact.field_compact(d)
      K = KC.base_field()
      r = 10

      a = KC.random()^2
      s = a.symbols2_log(0,r)
      self.assertEqual(vector(s), vector([0]*r))  

      a = KC(K.random()^2)
      s = a.symbols2_log(0,r)
      self.assertEqual(vector(s), vector([0]*r))

      a = KC.random()
      b = KC.random()
      ab = a*b

      a_s = vector(GF(2), a.symbols2_log(0,r))
      b_s = vector(GF(2), b.symbols2_log(0,r))
      ab_s = vector(GF(2), ab.symbols2_log(0,r))
      self.assertEqual(a_s + b_s, ab_s)
    
  def test_compactify(self):
    for d in self.fields:
      KC = field_compact.field_compact(d)
      K = KC.base_field()
      U = units.generators(d)

      a = K.random()

      # testing joining of equal elements
      b = KC((a,a), (1,1))
      b_c = b.compactify()
      self.assertTrue(len(b_c.elements) == 1 and len(b_c.powers) == 1)
      self.assertTrue(b_c.elements[0] == a)
      self.assertTrue(b_c.powers[0] == 2)

      # testing joining of elements that differ by a unit
      u = prod(rnd.choices(U, k=3)).element
      a = K.random()

      b = KC([a,a*u], (1,1))
      b_c = b.compactify(mod_units=True)
      if len(b_c.elements) == 2:
        # the case when reduction mod units failed
        self.assertTrue(b_c.powers[0] == 2 or b_c.powers[1] == 2)
        self.assertTrue(b_c.elements[0].absnorm() in [-1,1] or b_c.elements[1].absnorm() in [-1,1])
      else:
        self.assertTrue(b_c.powers[0] == 2)
        self.assertEqual(abs(b_c.elements[0].absnorm()), abs(a.absnorm()))
        
if __name__ == '__main__':
  unittest.main()
