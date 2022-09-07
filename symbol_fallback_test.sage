import goodprime
import ideal

d = (2,3,5)

for pos in range(0,100):
  print d,pos
  q,s,sinv = goodprime.sqrtprime(d,pos)
  J = ideal.ideals(d)(q,s)
  assert abs(J.generator().absnorm()) == q
