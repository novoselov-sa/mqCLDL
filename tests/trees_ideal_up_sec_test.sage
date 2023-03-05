# This is a hack to run unittests from parent directory
import os
import sys
sys.path.append(os.path.dirname(__file__) + "/../")

import unittest
import field
import prime_decomp
import trees

class TestTreesNewSecFieldMethods(unittest.TestCase):
  def test_ideal_up_sec_field(self):
    #d_L = (-7, -11, -19, -23)
    #d_K = (-7, -11, (-19)*(-23))
    #primes = [163, 443, 463, 499, 653, 823, 883, 947, 1061, 1087, 1373, 1453, 1499, 1787, 1871]
    d_L = (-11, -19, -31)
    d_K = (-11, (-19)*(-31))
    primes = [47, 157, 163, 191, 311, 397, 419, 443, 467, 577]
    K = field.field(d_K)
    L = field.field(d_L)
    L_sage = NumberField([x^2 - di for di in d_L], names="a")
    K_sage = NumberField([x^2 - di for di in d_K])

    for p in primes:
      L.sq_basis_mod(p)
      Fp = GF(p)
      F = prime_decomp.prime_decomp_fast(d_K, p)
      for I,e in F:
        I_up = trees.ideal_up(d_K, d_L, I)

        for J,e in I_up:
          self.assertEqual(J[1], L.abs_gen() + lift(L.abs_gen().apply_aut(J[2][0]).evaluate_mod_ext(p)), "Wrong automorphism!")

        I_up = [[L_sage.ideal(J[0], L_sage(str(J[1].to_sage(names='a')))),e] for J,e in I_up]
        #print(f"I_up = {I_up}")

        p_factors_over_L = [J for J,e in L_sage.factor(p)]

        for J,e in I_up:
          self.assertIn(J, p_factors_over_L)
        
        # check that I = I1 * I2 in L
        J = prod(J1^e for J1,e in I_up)
        I_sage = L_sage.ideal(I[0], L_sage(str(L(I[1]).to_sage(names='a'))))

        #print "I = ", I_sage
        #print "I1 * I2 =", J
        self.assertEqual(I_sage, J)

  def test_ideal_up_sec_field_apply_aut(self):
    d_L = (-7, -11, -19, -23)
    d_K = (-7, -11, (-19)*(-23))

    K = field.field(d_K)
    L = field.field(d_L)
    L_sage = NumberField([x^2 - di for di in d_L], names="a")

    # primes p s.t. sqrt(d_i) in F_p
    for p in [
      163, 443, 463, 499, 653, 823, 883, 947, 1061, 1087, 1373, 1453, 1499, 1787, 1871
    ]:
      Fp = GF(p)
      L.sq_basis_mod(p)
      F = prime_decomp.prime_decomp_fast(d_K, p)
      for I,e in F:
        I_up = trees.ideal_up(d_K, d_L, I, aut=True)
        I_sigma = [I[0], I[1].apply_aut([0]*(len(d_K)-1) + [1]), I[2]]

        I_up = [[L_sage.ideal(J[0], L_sage(str(J[1].to_sage(names='a')))),e] for J,e in I_up]
        p_factors_over_L = [J for J,e in L_sage.factor(p)]

        for J,e in I_up:
          self.assertIn(J, p_factors_over_L)
        
        # check that sigma(I) = I1 * I2 in L
        J = prod(J1^e for J1,e in I_up)
        I_sage = L_sage.ideal(I_sigma[0], L_sage(str(L(I_sigma[1]).to_sage(names='a'))))

        #print "I = ", I_sage
        #print "I1 * I2 =", J
        self.assertEqual(I_sage, J)

  # Testing lift of prime ideals from
  #   Q[sqrt(d_1), ..., sqrt(d_{n-2}), sqrt(d_{n-1})*sqrt(d_{n}] 
  # to
  #   Q[sqrt(d_1), ..., sqrt(d_{n})].
  def test_ideal_up_non_linear_sec_field(self):
    d_L = (-11, -19, -31)
    d_K = (-11, (-19)*(-31))
    K = field.field(d_K)
    L = field.field(d_L)
    L_sage = NumberField([x^2 - di for di in d_L], names="a")

    for p in Primes()[2:100]:
      # p should not divide [O_L:Z[theta_L]], [O_K:Z[theta_K]], and the discriminants
      if Mod(K.idx(), p) == 0 or Mod(L.idx(), p) == 0:
        continue
      Fp = GF(p)
      #print(f"s_K = {s_K}")
      F = prime_decomp.prime_decomp_fast(d_K, p)

      for I,e in F:
        I_up = trees.ideal_up(d_K, d_L, I)

        for J,e in I_up:
          #print(f"J = {J}")
          p0,g0,m0 = J
          gamma,typ = m0
          theta = L.abs_gen()
          theta1 = sum([L.gens()[i] for i in range(len(d_L)) if Mod(d_L[i],p).is_square()], L.zero())
          theta2 = sum([L.gens()[i] for i in range(len(d_L)) if not Mod(d_L[i],p).is_square()], L.zero())
          theta1_a = theta1.apply_aut(gamma)
          theta2_a = theta2.apply_aut(gamma)
          theta1_a_p = theta1_a.evaluate_mod_ext(p)
          if typ == 'quadratic':
            theta2_a_sq = lift((theta2_a^2).evaluate_mod_ext(p))
            g = theta^2 - theta * 2 * lift(theta1_a_p) + lift(theta1_a_p^2) - theta2_a_sq
          else:
            g = theta + lift(theta1_a_p)
          self.assertEqual(g0, g, "Wrong automorphism in ideal representation!")

        I_up = [[L_sage.ideal(J[0], L_sage(str(J[1].to_sage(names='a')))), e] for J,e in I_up]
        p_L_factors_sage = [J for J,e in L_sage.factor(p)]

        for J,e in I_up:
          self.assertIn(J, p_L_factors_sage, "Wrong factor in ideal_up: it should be in the decomposition of p in O_L")

        I_sage = L_sage.ideal(I[0], L_sage(str(L(I[1]).to_sage(names='a'))))
        I_L_factors_sage = [F for F,e in I_sage.factor()]

        if len(I_up) == 1:
          self.assertEqual(len(I_L_factors_sage), 1)
          self.assertIn(I_up[0][0], I_L_factors_sage, "Wrong ideal in ideal_up: ideal is not above input ideal")
        else:
          self.assertTrue(I_up[0][0] in I_L_factors_sage or I_up[1][0] in I_L_factors_sage, "Wrong ideals in ideal_up: both ideals are not above input ideal")

        for J,e in I_up:
          self.assertIn(J, I_L_factors_sage, "Wrong factor in ideal_up: it should be in the decomposition of I in O_L")

        # check that I = I1 * I2 in L
        J = prod(J1^e for J1,e in I_up)
        self.assertEqual(I_sage, J, "Wrong factors in ideal_up: the product should be equal to the input ideal")
        #print(f"DEBUG: I = {I}, p = {p} passed")

  def test_ideal_up_non_linear_sec_field_apply_aut(self):
    d_L = (-5, -13, -17, -29)
    d_K = (-5, -13, (-17)*(-29))

    K = field.field(d_K)
    L = field.field(d_L)
    L_sage = NumberField([x^2 - di for di in d_L], names="a")

    for p in Primes()[2:100]:
      if Mod(K.idx(), p) == 0 or Mod(L.idx(), p) == 0:
        continue

      L.sq_basis_mod(p) # fix square roots in base field
      assert field.fixed_sqrt(d_L[-1], p)*field.fixed_sqrt(d_L[-2], p) == field.fixed_sqrt(d_L[-1]*d_L[-2], p)

      Fp = GF(p)
      F = prime_decomp.prime_decomp_fast(d_K, p)

      for I,e in F:
        I_up = trees.ideal_up(d_K, d_L, I, aut=True)
        I_sigma = [I[0], I[1].apply_aut([0]*(len(d_K)-1) + [1]), I[2]]

        I_up = [[L_sage.ideal(J[0], L_sage(str(J[1].to_sage(names='a')))),e] for J,e in I_up]
        p_factors_over_L = [J for J,e in L_sage.factor(p)]

        for J,e in I_up:
          self.assertIn(J, p_factors_over_L)

        # check that sigma(I) = I1 * I2 in L
        J = prod(J1^e for J1,e in I_up)
        I_sage = L_sage.ideal(I_sigma[0], L_sage(str(L(I_sigma[1]).to_sage(names='a'))))

        #print "I = ", I_sage
        #print "I1 * I2 =", J
        self.assertEqual(I_sage, J)


if __name__ == '__main__':
    unittest.main()
