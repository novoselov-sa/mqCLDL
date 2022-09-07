import field
import cryptosystem
import ideal

for n in range(7):
  print(n)
  N = 2^n

  for j in [1,n,n^2]:
    d = tuple(primes(j,j + 10000))[:n]
    for loop in range(10):
      params = cryptosystem.parameters(d)
      pubkey,seckey = cryptosystem.keygen(d,params)
      q,qs = pubkey
      I = ideal.ideals(d)(q,qs)
      attackg = I.generator()
      print('attack',attackg)
      K = field.field(d)
      g = K(seckey[0].c,1)
      assert g * attackg.divexact(g) == attackg
      assert attackg * g.divexact(attackg) == g
      print('secret',g)
      ratio = g.divexact(attackg)
      print('ratio',ratio)
      if n > 0:
        missing = ideal.exponents_rounded(d,ratio,1)
        print('missing',missing)
        check = ideal.shorten(d,ratio)
        one = K((1,)+(0,)*(N-1),1)
        assert check == one or check == -one
      print('canskipenumeration',g == attackg or g == -attackg)
