# This is a hack to run unittests from parent directory
import os
import sys
sys.path.append(os.path.dirname(__file__) + "/../")

import unittest
import field

class TestFieldMethods(unittest.TestCase):

  def test_absolute_polynomial(self):
    fields = [
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
    for d in fields:
      K = field.field(d)
      pol = K.absolute_polynomial().to_sage().change_ring(QQ)
      K2 = NumberField([x^2-d[i] for i in range(len(d))], names='a')
      self.assertEqual(K2.absolute_polynomial(), pol)

  def test_field_discriminant(self):
    self.assertEqual(field.field((2, 3, 5)).discriminant(), 3317760000)
    self.assertEqual(field.field((-2, 3, 5)).discriminant(), 3317760000)
    self.assertEqual(field.field((3, 5, 11)).discriminant(), 189747360000)
    self.assertEqual(field.field((5, 13, 17)).discriminant(), 1490902050625)
    self.assertEqual(field.field((-5, -13, 17)).discriminant(), 381670924960000)
    self.assertEqual(field.field((-3, -7, -11)).discriminant(), 2847396321)
    self.assertEqual(field.field((-3*5, -7, -11)).discriminant(), 1779622700625)
  
  def test_field_pow(self):
    fields = [
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
    for d in fields:
      K = field.field(d)
      self.assertEqual(K.zero()^0, K.one())
      v = ZZ.random_element(1000)
      a = ZZ.random_element(100)
      self.assertEqual(K.from_ZZ(v)^a, K.from_ZZ(v^a))
  
  def test_add_int(self):
    fields = [
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
    for d in fields:
      K = field.field(d)
      self.assertEqual(K.zero() + 1, K.one())
      a = ZZ.random_element(1000)
      b = ZZ.random_element(1000)
      self.assertEqual(K.from_ZZ(a) + b, K.from_ZZ(a+b))
      self.assertEqual(K.from_ZZ(a) - b, K.from_ZZ(a-b))

  def test_mul_int(self):
    fields = [
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
    for d in fields:
      K = field.field(d)
      self.assertEqual(K.zero()*2, K.zero())
      self.assertEqual(K.one()*0, K.zero())
      a = ZZ.random_element(1000)
      b = ZZ.random_element(1000)
      self.assertEqual(K.from_ZZ(a) * b, K.from_ZZ(a*b))
  
  def test_fixed_sqrt_mult(self):
    di = [2, -3, 5, 7*11, -13*17, 19*23*47]
    for p in primes(50):
      for a in di:
         for b in di:
            a_s = field.fixed_sqrt(a, p)
            b_s = field.fixed_sqrt(b, p)
            ab_s = field.fixed_sqrt(a*b, p)
            self.assertEqual(a_s*b_s, ab_s, f"choice of square roots should be compatible with multiplication, a = {a}, b = {b}")

  def test_fixed_sqrt_mult_reverse(self):
    di = reversed([2, -3, 5, 7*11, -13*17, 19*23*47])
    for p in primes(50):
      for a in di:
        for b in di:
            a_s = field.fixed_sqrt(a, p)
            b_s = field.fixed_sqrt(b, p)
            ab_s = field.fixed_sqrt(a*b, p)
            self.assertEqual(a_s*b_s, ab_s, f"choice of square roots should be compatible with multiplication, a = {a}, b = {b}")


if __name__ == '__main__':
  unittest.main()
