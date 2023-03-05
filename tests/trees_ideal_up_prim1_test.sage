# This is a hack to run unittests from parent directory
import os
import sys
sys.path.append(os.path.dirname(__file__) + "/../")

import unittest
import field
import prime_decomp
import trees

class TestTreesNewMethods(unittest.TestCase):

  def test_ideal_up_prim_field_1(self):
    #d_L = (-7, -11, -19, -23)
    #d_K = (-7, -11, -19)
    #primes = [163, 443, 463, 499, 653, 823, 883, 947]
    d_L = (-11, -19, -31)
    d_K = (-11, -19)
    # primes s.t. sqrt(d_i) in Fp and p does not divide [O_L:Z[theta]] and [O_L:Z[theta_K]] and p spits completely in O_K.
    primes = [47, 157, 163, 191, 311, 397, 419, 443, 467, 577]

    K = field.field(d_K)
    L = field.field(d_L)
    L_sage = NumberField([x^2 - di for di in d_L], names="a")

    for p in primes:
      Fp = GF(p)
      F = prime_decomp.prime_decomp_fast(d_K, p)
      for I,e in F:
        I_up = trees.ideal_up(d_K, d_L, I)
        I_up = [[L_sage.ideal(J[0], L_sage(str(J[1].to_sage(names='a')))), e] for J,e in I_up]
        p_factors_over_L = [J for J,e in L_sage.factor(p)]

        for J,e in I_up:
          self.assertIn(J, p_factors_over_L)

        # check that I = I1 * I2 in L
        J = prod(J1^e for J1,e in I_up)
        I_sage = L_sage.ideal(I[0], L_sage(str(L(I[1]).to_sage(names='a'))))
        self.assertEqual(I_sage, J)
  
  # Testing conversion from:
  #   pO_K = prod (p, ( theta_s - (mu(theta_s') mod p))^2 - (mu(theta_s'')^2 mod p) )) for mu in Gal(K) / (mu(theta_s'') = -mu(theta_s''))
  # to
  #   pO_L = prod (p, ( theta - (gamma(theta') mod p))^2 - (gamma(theta'')^2 mod p) )) for gamma in Gal(L) / (gamma(theta'') = -gamma(theta'')).
  def test_ideal_up_non_linear_case_field_1(self):
    d_L = (-11, -19, -23, -31)
    d_K = (-11, -19, -23)
    K = field.field(d_K)
    L = field.field(d_L)
    L_sage = NumberField([x^2 - di for di in d_L], names="a")

    for p in Primes()[2:100]:
      # we select primes p s.t. p O_K splits into deg(K)/2 prime ideals
      s_K = [di for di in d_K if not is_square(Mod(di, p))]
      # linear->[linear] or linear -> quadratic
      #if len(s_K) == 0:
      #  continue

      # p should not divide [O_L:Z[theta_L]], [O_K:Z[theta_K]], and the discriminants
      if Mod(K.idx(), p) == 0 or Mod(L.idx(), p) == 0: #or Mod(K.discriminant(), p) == 0 or Mod(L.discriminant(), p) == 0:
        continue

      Fp = GF(p)
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
        
        #print(f"I_up = {I_up}")
        I_up = [[L_sage.ideal(J[0], L_sage(str(J[1].to_sage(names='a')))), e] for J,e in I_up]
        p_L_factors_sage = [J for J,e in L_sage.factor(p)]

        for J,e in I_up:
          self.assertIn(J, p_L_factors_sage, "Wrong factor in ideal_up: it should be in the decomposition of p in O_L")

        I_sage = L_sage.ideal(I[0], L_sage(str(L(I[1]).to_sage(names='a'))))
        I_L_factors_sage = [F for F,e in I_sage.factor()]

        #print(f"I_sage.factor() = {I_sage.factor()}")
        #print(f"I_up (sage) = {I_up}")

        if len(I_up) == 1:
          self.assertEqual(len(I_L_factors_sage), 1)
          self.assertIn(I_up[0][0], I_L_factors_sage, "Wrong ideals in ideal_up: ideals are not above input ideal")
        else:
          self.assertTrue(I_up[0][0] in I_L_factors_sage or I_up[1][0] in I_L_factors_sage, "Wrong ideals in ideal_up: both ideals are not above input ideal")

        for J,e in I_up:
          self.assertIn(J, I_L_factors_sage, "Wrong factor in ideal_up: it should be in the decomposition of I in O_L")

        # check that I = I1 * I2 in L
        J = prod(J1^e for J1,e in I_up)
        self.assertEqual(I_sage, J, "Wrong factors in ideal_up: the product should be equal to the input ideal")

if __name__ == '__main__':
    unittest.main()
