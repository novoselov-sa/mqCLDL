import nprofile
import ring
import field
import polynomial_ring
import time

for n in range(8,9):
  print n
  N = 2^n
  d = (-2,3,-5,7,-11,13,-17,19,-23,29)[:n]
  K = field.field(d)
  PR = polynomial_ring.polynomial_ring(d, "x")
  print "d = ", d
  tm = time.time()
  K2 = NumberField([x^2-di for di in d])
  if len(d) > 0:
    min_poly_sage = K2.absolute_polynomial()
    tm = time.time() - tm
    print "Min.poly. (sage):", min_poly_sage
    print "-> computed in ", tm

  tm = time.time()
  min_pol = K.absolute_polynomial()
  tm = time.time() - tm
  print "K.absolute_polynomial(): ", min_pol.to_sage()
  print "-> computed in ", tm
  #tm = time.time()
  #K2 = NumberField([x^2-di for di in d])
  if len(d) > 0:
  #  min_poly_sage = K2.absolute_polynomial()
  #  tm = time.time() - tm
  #  print "Min.poly. (sage):", min_poly_sage
  #  print "-> computed in ", tm
    assert min_poly_sage == min_pol.to_sage()

  for bits in range(20):
    for dg in range(5):
      pol1_c = [K(tuple(ZZ.random_element(1-2^bits,2^bits) for j in range(N)),ZZ.random_element(1,1+2^bits)) for i in range(dg+1)]
      pol2_c = [K(tuple(ZZ.random_element(1-2^bits,2^bits) for j in range(N)),ZZ.random_element(1,1+2^bits)) for i in range(dg+1)]
      pol1 = PR(pol1_c)
      pol2 = PR(pol2_c)
      #print "pol1: ", pol1.to_sage()
      #print "pol1.deg(): ", pol1.deg()
      #print "pol2: ", pol2.to_sage()
      #print "pol2.deg(): ", pol2.deg()

      pol3 = pol1*pol2
      #print "pol1*pol2: ", pol3.to_sage()
      #print "(pol1*pol2).deg(): ", pol3.deg()
      if pol3.deg() > 0:
        assert (pol1*pol2).deg() <= pol1.deg() + pol2.deg()
      #print "pol1*pol2(1): ", pol3.evaluate(K.one()).to_sage()
      assert pol1*pol2 == pol2*pol1
      #print "pol1+pol2 = ", (pol1+pol2).to_sage()
      assert pol1+pol2 == pol2+pol1
      assert pol1 - pol2 == pol1 + (-pol2)
      #print "(pol1 + pol2) - pol2 = ", ((pol1 + pol2) - pol2).to_sage()
      assert (pol1 + pol2) - pol2 == pol1

      assert pol1*PR.zero()== PR.zero()
      assert pol1-pol1 == PR.zero()
      #print "pol1*PR([K.one()]):", (pol1*PR.one()).to_sage()
      assert pol1*PR.one() == pol1
      assert -(-pol1) == pol1
      assert K.zero() + K.zero() == K.zero()
      assert K.zero() + K.one() == K.one()

      #f = K(tuple(ZZ.random_element(1-2^bits,2^bits) for j in range(N)),ZZ.random_element(1,1+2^bits))
      #print "f = ", f.to_sage()
      #print "f.minpoly() = ", f.minimal_polynomial().to_sage()
      #print "----------------"

print "done"
nprofile.output([polynomial_ring])