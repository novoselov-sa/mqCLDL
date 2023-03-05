# This is a hack to run unittests from parent directory
import os
import sys
sys.path.append(os.path.dirname(__file__) + "/../")

import unittest
import field
import prime_decomp

class TestPrimeDecompMethods(unittest.TestCase):

  def test_prime_decomp(self):
    d = (-3, -7, -11)
    K = field.field(d)
    proof.number_field(False)
    K2 = NumberField([x^2-d[i] for i in range(len(d))], names='a')

    idx = K.idx()

    for p in prime_range(2,100):
      #if abs(idx) % p == 0:
      #  continue
      D1 = prime_decomp.prime_decomp(d, p)
      D1 = [(K2.ideal(I[0][0], K2(str(I[0][1].to_sage(names='a')))), I[1]) for I in D1]
      D2 = K2.factor(p)

      for I in D1:
        self.assertIn(I, D2)

  def test_prime_decomp_fast_linear(self):
    d = (-7,-11,-19,-23)
    K = field.field(d)
    proof.number_field(False)
    K_sage = NumberField([x^2-d[i] for i in range(len(d))], names='a')

    idx = K.idx()

    for p in [163, 443, 463, 499, 653, 823, 883, 947]:
      if abs(idx) % p == 0:
        continue
      p_factors = prime_decomp.prime_decomp_fast(d, p)
      p_factors = [(K_sage.ideal(I[0][0], K_sage(str(I[0][1].to_sage(names='a')))), I[1]) for I in p_factors]
      p_factors_sage = K_sage.factor(p)

      for I in p_factors:
        self.assertIn(I, p_factors_sage)
      self.assertEqual(prod(I[0]^I[1] for I in p_factors), prod(I[0]^I[1] for I in p_factors_sage))


  def test_prime_decomp_fast_full(self):
    fields = [
      (5, 13),
      (-11, -31),
      (-11, -19, -31),
      (-11, -19, -23)
    ]
    for d in fields:
      K = field.field(d)
      disc = K.discriminant()
      idx = K.idx()
      
      proof.number_field(False)
      K2 = NumberField([x^2-d[i] for i in range(len(d))], names='a')

      for p in primes(250):
        #if abs(idx) % p == 0 or Mod(disc, p) == 0:
        if abs(idx) % p == 0:
          continue
        D1 = prime_decomp.prime_decomp_fast(d, p)
        D1 = [(K2.ideal(I[0][0], K2(str(I[0][1].to_sage(names='a')))), I[1]) for I in D1]
        D2 = K2.factor(p)

        for I in D1:
          self.assertIn(I, D2)
        self.assertEqual(prod(I[0]^I[1] for I in D1), prod(I[0]^I[1] for I in D2))

  def test_prime_decomp_fast_full_not_prime_di(self):
    fields = [
       (-3, 5*17),
       (7*11, 23*43),
       (3*5, 7*13, 3*13)
    ]
    for d in fields:
      K = field.field(d)
      proof.number_field(False)
      K_sage = NumberField([x^2-d[i] for i in range(len(d))], names='a')

      disc = K.discriminant()
      idx = K.idx()

      for p in primes(200):
        #if abs(idx) % p == 0 or Mod(disc, p) == 0:
        #  continue

        if abs(idx) % p == 0:
          continue
        D1 = prime_decomp.prime_decomp_fast(d, p)
        D1 = [(K_sage.ideal(F[0][0], K_sage(str(F[0][1].to_sage(names='a')))), F[1]) for F in D1]
        D2 = K_sage.factor(p)

        for F in D1:
          self.assertIn(F, D2)
        self.assertEqual(prod(F[0]^F[1] for F in D1), prod(F[0]^F[1] for F in D2))
  
  def test_prime_decomp_fast_gammas(self):
    #field.SQROOTS = {} - clearing SQROOT doesn't work without clearing memoized data 
    fields = [
       (-3, 5*17),
       (7*11, 23*43),
       (3, 5, 7),
       (-11, -19, -23),
       (-11, -31, -19),
       (3*5, 7*13, 3*13)
    ]
    for d in fields:
      K = field.field(d)

      n = len(d)
      idx = K.idx()
      disc = K.discriminant()

      for p in primes(250):
        #if abs(idx) % p == 0 or Mod(disc, p) == 0:
        #  continue

        if abs(idx) % p == 0:
          continue
        t1 = [K.gens()[i] for i in range(n) if Mod(d[i],p).is_square()]
        t2 = [K.gens()[i] for i in range(n) if not Mod(d[i],p).is_square()]
        theta1 = sum(t1, K.zero())
        theta2 = sum(t2, K.zero())
        theta = K.abs_gen()

        D1 = prime_decomp.prime_decomp_fast(d, p)
        for F in D1:
          I = F[0]
          p0, g, meta = I
          self.assertEqual(p0, p)
          gamma = meta[0]
          typ = meta[1]
          if typ == 'quadratic':
            theta1_a = theta1.apply_aut(gamma)
            theta2_a = theta2.apply_aut(gamma)
            theta1_a_p = theta1_a.evaluate_mod_ext(p)

            g0 = theta^2 -  theta * lift(theta1_a_p) * 2 + lift(theta1_a_p^2) - lift((theta2_a^2).evaluate_mod_ext(p))
            self.assertEqual(g, g0, f"Incorrect automorphism gamma = {gamma} in decomposition of p = {p}. ")

  def test_prime_decomp_slow_full(self):
    fields = [
      (5, 13),
      (-11, -31),
      (-11, -19, -31),
      (-11, -19, -23)
    ]
    for d in fields:
      K = field.field(d)
      disc = K.discriminant()
      idx = K.idx()
      
      proof.number_field(False)
      K2 = NumberField([x^2-d[i] for i in range(len(d))], names='a')

      for p in primes(250):
        #if abs(idx) % p == 0:
        #  print(f"d = {d}, p = {p}, p | idx = {idx} = {idx.factor()}")
        D1 = prime_decomp.prime_decomp_slow(d, p)
        D1 = [(K2.ideal(I[0][0], K2(str(I[0][1].to_sage(names='a')))), I[1]) for I in D1]
        D2 = K2.factor(p)

        for I in D1:
          self.assertIn(I, D2)
        self.assertEqual(prod(I[0]^I[1] for I in D1), prod(I[0]^I[1] for I in D2))

if __name__ == '__main__':
  unittest.main()
