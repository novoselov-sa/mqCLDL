# This is a hack to run unittests from parent directory
import os
import sys
sys.path.append(os.path.dirname(__file__) + "/../")

import unittest
import field
import prime_decomp
import trees

class TestTreesNewSecFieldMethods(unittest.TestCase):
  def test_ideal_up_hard_prim2_field(self):
    fields =    [(-3, -5, -7, -11, -13)] 
    subfields = [(-3, -5, -7, -13)]

    #fields =    [(-3, -5, 77)] 
    #subfields = [(-3, 77)]

    for d_L, d_K in zip(fields, subfields):
      K = field.field(d_K)
      L = field.field(d_L)
      L_sage = NumberField([x^2 - di for di in d_L], names="a")
      K_sage = NumberField([x^2 - di for di in d_K])

      for p in primes(1000):
        if not (Mod(K.idx(),p) == 0 or Mod(L.idx(),p) == 0):
          continue
        #print(f"d_L = {d_L}, d_K = {d_K}, p = {p}")
        print(f"p = {p}")

        L.sq_basis_mod(p)
        Fp = GF(p)
        F = prime_decomp.prime_decomp(d_K, p)

        for I,e in F:
          #print(f"I = {I}, e = {e}")
          I_up = trees.ideal_up_hard(d_K, d_L, I)

          #print(f"\nI_up = {I_up}")

          for J,e in I_up:
            p, g, meta = J
            if meta[1] != 'other':
              g0 = trees.ideal_gen_from_aut(d_K, meta[0], p, meta[1])
              self.assertEqual(g0, g, "Wrong automorphism!")

          I_up = [[L_sage.ideal(J[0], L_sage(str(J[1].to_sage(names='a')))),e] for J,e in I_up]

          p_factors_over_L = [J for J,e in L_sage.factor(p)]

          for J,e in I_up:
            self.assertIn(J, p_factors_over_L)

          # check that I = I1 * I2 in L
          J = prod(J1^e for J1,e in I_up)
          I_sage = L_sage.ideal(I[0], L_sage(str(L(I[1]).to_sage(names='a'))))

          self.assertEqual(I_sage, J)

if __name__ == '__main__':
    unittest.main()
