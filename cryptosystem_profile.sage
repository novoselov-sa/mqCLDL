import cryptosystem

n = 4
if len(sys.argv) > 1: n = ZZ(sys.argv[1])

scale = 0.5
if len(sys.argv) > 2: scale = RR(sys.argv[2])

d = ()
for p in primes(0,infinity):
  if len(d) >= n: break
  d += (p,)

params = cryptosystem.parameters(d,scale)
for loops in range(100):
  pubkey,seckey = cryptosystem.keygen(d,params)
  expandedseckey = cryptosystem.expandseckey(d,seckey)
  message = cryptosystem.randommessage(d,params,1024.0)
  ciphertext = cryptosystem.encrypt(d,message,pubkey)
  decrypt = cryptosystem.decrypt(d,ciphertext,expandedseckey)
  assert message == decrypt

import nprofile
import goodprime
import mult
import div
nprofile.output([goodprime,mult,div,cryptosystem])
