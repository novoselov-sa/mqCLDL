# This is a hack to run unittests from parent directory
import os
import sys
sys.path.append(os.path.dirname(__file__) + "/../")

import unittest
import field
import prime_decomp
import trees

class TestTreesNewSecFieldMethods(unittest.TestCase):
  def test_ideal_up_hard_sec_field(self):
    #fields =    [(-11, -19, -31),    (-7, -11, -19, -23),   (-3, -5, -7)] 
    #subfields = [(-11, (-19)*(-31)), (-7, -11, (-19)*(-23)), (-3, -7) ]
    fields =    [(-11, -19, -31)] 
    subfields = [(-11, (-19)*(-31))]

    for d_L, d_K in zip(fields, subfields):
      K = field.field(d_K)
      L = field.field(d_L)
      L_sage = NumberField([x^2 - di for di in d_L], names="a")
      K_sage = NumberField([x^2 - di for di in d_K])

      for p in primes(1000):
        if not (Mod(K.idx(),p) == 0 or Mod(L.idx(),p) == 0):
          continue
        print(f"d_L = {d_L}, d_K = {d_K}, p = {p}")
        L.sq_basis_mod(p)
        Fp = GF(p)
        F = prime_decomp.prime_decomp(d_K, p)
        
        for I,e in F:
          I_up = trees.ideal_up_hard(d_K, d_L, I)

          for J,e in I_up:
            p, alpha, meta = J
            if meta[1] != 'other':
              self.assertEqual(J[1], L.abs_gen() + lift(L.abs_gen().apply_aut(J[2][0]).evaluate_mod_ext(p)), "Wrong automorphism!")

          I_up = [[L_sage.ideal(J[0], L_sage(str(J[1].to_sage(names='a')))),e] for J,e in I_up]

          p_factors_over_L = [J for J,e in L_sage.factor(p)]

          for J,e in I_up:
            self.assertIn(J, p_factors_over_L)
          
          # check that I = I1 * I2 in L
          J = prod(J1^e for J1,e in I_up)
          I_sage = L_sage.ideal(I[0], L_sage(str(L(I[1]).to_sage(names='a'))))

          self.assertEqual(I_sage, J)

  def test_ideal_up_hard_sec_field_apply_aut(self):
    d_L = (-7, -11, -19, -23)
    d_K = (-7, -11, (-19)*(-23))

    #d_L = (-3, -5, -7, -11)
    #d_K = (-3, -5, 77)
  
    #d_L = (-3, -7, -11, 437)
    #d_K = (-3, -7, -4807)
    
    #d_L = (-3, -7, -11, 437)
    #d_K = (-3, -7, -11)

    #d_K = (-3, -385)
    #d_L = (-3, -5, 77)

    K = field.field(d_K)
    L = field.field(d_L)
    L_sage = NumberField([x^2 - di for di in d_L], names="a")

    print(f"K.idx() = {K.idx().factor()}")
    print(f"L.idx() = {L.idx().factor()}")

    for p in primes(1000):
      # we choose "hard" primes
      if not (Mod(K.idx(),p) == 0 or Mod(L.idx(),p) == 0):
        continue
      print(f"\np = {p}")
      Fp = GF(p)
      L.sq_basis_mod(p)
      F_K = prime_decomp.prime_decomp(d_K, p)
      for I,e in F_K:
        I_up = trees.ideal_up_hard(d_K, d_L, I, aut=True)

        #print(f"I_up (orig.) = {I_up}")

        # check correctness of the automorphisms
        for J,e in I_up:
          p0, g0, meta = J
          gamma,typ = meta
          if typ != 'other':
            g = trees.ideal_gen_from_aut(d_L, gamma, p, typ)
            self.assertEqual(g0, g, "Wrong automorphism in ideal representation!")

        I_up = [[L_sage.ideal(J[0], L_sage(str(J[1].to_sage(names='a')))),e] for J,e in I_up]

        # check that result of ideal_up belongs to factors of p
        p_factors_over_L = [J for J,e in L_sage.factor(p)]
        for J,e in I_up:
          self.assertIn(J, p_factors_over_L, "Wrong factor in ideal_up_hard: it should be in the decomposition of p in O_L")

        # check that sigma(I) = I1 * I2 in L, where sigma is the automorphism
        if I[2][1] == 'other':
          I_sigma = [I[0], I[1].apply_aut([0]*(len(d_K)-1) + [1]), I[2]]
        else:
          p0,g0,m0 = I
          gamma,typ = m0
          gamma_s = gamma[:-1] + [lift(Mod(gamma[-1] + 1, 2))]
          g_s = trees.ideal_gen_from_aut(d_K, gamma_s, p, typ)
          I_sigma = [p0, g_s, [gamma_s, typ]]
        
        I_sigma_sage = L_sage.ideal(p, L_sage(str(L(I_sigma[1]).to_sage(names='a'))))
        J = prod(J1^e for J1,e in I_up)

        #print(f"I = {I}\n")
        #print(f"e = {e}\n")
        #print(f"I_sigma = {I_sigma}\n")
        #print(f"I_sigma_sage = {I_sigma_sage}\n")
        #print(f"I_sigma_sage.factor() = {I_sigma_sage.factor()}\n")
        #print(f"I_up = {I_up}\n")
        #print(f"I1 * I2 = {J}\n")
        self.assertEqual(I_sigma_sage, J, "Wrong factors in ideal_up_hard: the product should be equal to the input ideal")


if __name__ == '__main__':
    unittest.main()
