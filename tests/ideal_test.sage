# This is a hack to run unittests from parent directory
import os
import sys
sys.path.append(os.path.dirname(__file__) + "/../")

import unittest
import field
import ideal
import cryptosystem

class TestIdealMethods(unittest.TestCase):

  def test_hnf(self):
    fields = [
      (5, 13),
      (-11, -31),
      (-3, -7, -11),
      (-11, -31, -19),
      (-7, -11, -19, -23),
      (-11, -19, -31, -47),
      (7, 11, 19, 23),
      (-3, 5*17),
      (7*11, 23*43),
      (3*5, 7*13, 3*13)
    ]
    for d in fields:
      K = field.field(d)
      for i in range(10):
        I = ideal.random(d)
        I_sage = I.to_sage()
        # simple test by comparing the norms
        self.assertEqual(I_sage.absolute_norm(), I.q)
        # test by choosing random elements from one ideal and check that it belongs to another
        for j in range(10):
          #r = K_sage(str(I.random_element().to_sage(names="a")))
          r = I.random_element().to_sage(K.sage())
          self.assertIn(r, I_sage, f"random elements should belong to the result of conversion, but {r} doesn't belong to the ideal {I_sage}")

  def test_generator_real_field(self):
    fields = [
      (5, 13),
      (3, 5, 13),
      (7, 11, 19, 23),
      (7*11, 23*43),
      (3*5, 7*13, 3*13)
    ]
    for d in fields:
      K = field.field(d).sage()
      while True:
        params = cryptosystem.parameters(d)
        pubkey,seckey = cryptosystem.keygen(d,params)
        q,qs = pubkey
        try:
          I = ideal.ideals(d)(q,qs) 
          break
        except NotImplementedError: # sometimes keygen produces unsupported parameters
          pass
      g = I.generator()
      self.assertEqual(K.ideal(g.to_sage(K)), I.to_sage(), f"Wrong computation of generator for ideal I = {I} != <{g}>")

  def test_generator_imaginary_field(self):
    fields = [
      (-3, -7),
      (-11, -19, -23),
      (-3, -7, -11*-19),
      (-31, -43, -47, -59),
      (-5, -13, -17, -29),
      (-61, -67, -71)
    ]
    for d in fields:
      K = field.field(d).sage()
      #print(f"field: {d}")
      while True:
        params = cryptosystem.parameters(d)
        pubkey,seckey = cryptosystem.keygen(d,params)
        q,qs = pubkey
        try:
          I = ideal.ideals(d)(q,qs) 
          break
        except NotImplementedError: # sometimes keygen produces unsupported parameters
          pass
      g = I.generator()
      #print(f"I = {I} = < {g.to_sage()} >")
      self.assertEqual(K.ideal(g.to_sage(K)), I.to_sage(), f"Wrong computation of generator for ideal I = {I} != <{g}>")

if __name__ == '__main__':
  unittest.main()
