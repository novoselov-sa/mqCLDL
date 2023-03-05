import cryptosystem
import ring

for n in range(0,8):
  print n
  N = 2^n

  d = (2,3,5,7,11,13,17,19,23,29)[:n]
  params = cryptosystem.parameters(d)

  for loop in range(10):
    pubkey,seckey = cryptosystem.keygen(d,params)
    expandedseckey = cryptosystem.expandseckey(d,seckey)

    output = '%s decrypt' % n

    for gap in [1.0,2.0,4.0,8.0,16.0,32.0,64.0,128.0,256.0,512.0,1024.0]:
      success = 0
      for loop2 in range(100):
        message = cryptosystem.randommessage(d,params,gap)
        ciphertext = cryptosystem.encrypt(d,message,pubkey)
        decrypt = cryptosystem.decrypt(d,ciphertext,expandedseckey)
        success += message == decrypt
      output += ' %.2f' % (success / 100.0)

    print output

  for loop in range(10):
    g = cryptosystem.keygen_g(d,params)
    q = abs(g.absnorm())
    for p in primes(2,params[0]):
      assert q % p != 0

  for badprimeslimit in [3,5,7]:
    params = badprimeslimit,prod(primes(2,badprimeslimit)),[0]*N
    for loop in range(10):
      g = cryptosystem.keygen_g(d,params)
      q = abs(g.absnorm())
      for p in primes(2,params[0]):
        assert q % p != 0

  for badprimeslimit in range(20):
    badprimes = prod(primes(2,badprimeslimit))
    if n <= 4:
      if badprimes^N <= 10000:
        target = set()
        for u in range(badprimes^N):
          u = tuple(((u//badprimes^j+badprimes//2)%badprimes)-badprimes//2 for j in range(N))
	  q = ring.ring(d)(u).absnorm()
	  if all(q % p != 0 for p in primes(2,badprimeslimit)):
            target.add(u)
        
        found = set()
        params = badprimeslimit,badprimes,[0]*N
        for loop in range(10*N*badprimes^N):
          g = cryptosystem.keygen_g(d,params)
          found.add(g.c)
          if len(found) == len(target): break
    
        assert found == target

  for p in primes(2,10):
    if n <= 4:
      if p^N <= 10000:
        d = tuple(ZZ.random_element(100) for i in range(n))

	target = set()
        for u in range(p^N):
          u = tuple((u//p^j)%p for j in range(N))
	  q = ring.ring_without_checking(d)(u).absnorm()
	  if q % p != 0:
	    target.add(tuple(GF(p)(uj) for uj in u))

	found = set()
        for loop in range(10*N*p^N):
	  g = cryptosystem.random_invertible_element(d,GF(p))
	  found.add(tuple(g))
	  if len(found) == len(target): break

	assert found == target
	print n,'invertible',p,d,len(found),'out of',p^N
