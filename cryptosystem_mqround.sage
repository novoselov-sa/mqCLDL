import sys
import field
import units
import ideal
import cryptosystem

loops = 10000

for n in range(2,9):
  N = 2^n

  B = prime_range(100)
  limit = len(B) - n + 1

  for i in range(limit):
    d = tuple(B[i:i + n])
    K = field.field(d)
    params = cryptosystem.parameters(d)

    U = units.mqgenerators_mod_torsion(d)

    numok = 0
    for loop in range(loops):
      pubkey,seckey = cryptosystem.keygen(d,params)
      g = seckey[0]
      logg = ideal.approxlog(d,K(g,1))

      ok = True
      for u in U:
        uu = sum(u.approxlog[j]^2 for j in range(N))
        ulogg = sum(u.approxlog[j]*logg[j] for j in range(N))
	if ulogg * 2 >= uu or ulogg * 2 <= -uu:
	  ok = False
	  break

      numok += ok

    print n,d[0],numok * 1.0 / loops
    sys.stdout.flush()
