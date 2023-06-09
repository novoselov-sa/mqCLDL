# This is a hack to run unittests from parent directory
import os
import sys
sys.path.append(os.path.dirname(__file__) + "/../")

import unittest
import field
import field_compact
import idealprod
import fb

class TestIdealProdMethods(unittest.TestCase): 
  
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
    self.primes = list(primes(10))
    self.abc = [chr(i) for i in range(ord('a'), ord('z')+1)]

  def test_mult(self):
    for d in self.fields:
      food = {}
      for i in range(len(d)):
        food[self.abc[i]] = d[i]
      K = field.field(d).sage(food)
      P = fb.factorsb(d, food=food)
      FB = []
      for p in self.primes:
        for I,e in K.factor(p):
          assert I.gens()[0] == p
          FB.append(P(p, str(I.gens()[1]), [1]*2^(len(d)+1)))

      IPS = idealprod.idealprods(d, FB=FB, food=food)
      a = IPS.random()
      b = IPS.random()
      self.assertEqual((a^2 * b^3), a^2 * b^3)

  def test_absnorm(self):
    for d in self.fields:
      food = {}
      for i in range(len(d)):
        food[self.abc[i]] = d[i]
      K = field.field(d).sage(food)
      P = fb.factorsb(d, food=food)
      FB = []
      for p in self.primes:
        for I,e in K.factor(p):
          assert I.gens()[0] == p
          FB.append(P(p, str(I.gens()[1]), [1]*2^(len(d)+1)))

      IPS = idealprod.idealprods(d, FB=FB, food=food)
      a = IPS.random(bits = 2, exp_bound=2, elems=1)
      self.assertEqual(abs(a.absnorm()), abs(a.evaluate().absolute_norm()))

  def test_valuation(self):
    for d in self.fields:
      if len(d) > 3: # tests are to slow for len(d) > 3
        continue
      food = {}
      for i in range(len(d)):
        food[self.abc[i]] = d[i]
      K = field.field(d).sage(food)
      P = fb.factorsb(d, food=food)
      FB = []
      for p in self.primes:
        for I,e in K.factor(p):
          assert I.gens()[0] == p
          FB.append(P(p, str(I.gens()[1]), [1]*2^(len(d)+1)))
      IPS = idealprod.idealprods(d, FB=FB, food=food)
      a = IPS.random(bits = 2, exp_bound=2, elems=1)
      for p,e0 in a.absnorm().factor():
        for I,e1 in K.factor(p):
          self.assertEqual(a.valuation(I), a.evaluate().valuation(I))

  def test_sqrt(self):
    for d in self.fields:
      food = {}
      for i in range(len(d)):
        food[self.abc[i]] = d[i]
      K = field.field(d).sage(food)
      P = fb.factorsb(d, food=food)
      FB = []
      for p in self.primes:
        for I,e in K.factor(p):
          assert I.gens()[0] == p
          FB.append(P(p, str(I.gens()[1]), [1]*2^(len(d)+1)))

      IPS = idealprod.idealprods(d, FB=FB, food=food)
      a = IPS.random(bits = 2, exp_bound=2, elems=2)
      s = a * a
      r = s.sqrt()
      self.assertEqual(r, a)

  def test_pow(self):
    for d in self.fields:
      food = {}
      for i in range(len(d)):
        food[self.abc[i]] = d[i]
      K = field.field(d).sage(food)
      P = fb.factorsb(d, food=food)
      FB = []
      for p in self.primes:
        for I,e in K.factor(p):
          assert I.gens()[0] == p
          FB.append(P(p, str(I.gens()[1]), [1]*2^(len(d)+1)))

      IPS = idealprod.idealprods(d, FB=FB, food=food) 
      a = IPS.random()
      b = IPS.random()

      self.assertEqual(a^3 * b^3, (a*b)^3)
      self.assertEqual(a^(-3) * b * a^2 * b^(-1), a^(-1))

if __name__ == '__main__':
  unittest.main()
