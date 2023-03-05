import sys
import field
import units
import ideal
import cryptosystem

nrange = range(2,9)
loops = 10
gloops = 1000

for n in nrange:
  N = 2^n

  B = prime_range(2*n*n)
  for loop in range(loops):
    d = tuple(sorted(Subsets(B,n).random_element()))
    K = field.field(d)
    assert len(d) == n
    print '% d',d
    sys.stdout.flush()
    params = cryptosystem.parameters(d)

    Q = units.mqgenerators_mod_torsion(d)

    D = [prod(d[t] for t in range(n) if i & (1 << t)) for i in range(N)]
    sqrtD = [sqrt(RR(Di)) for Di in D]
    logQ = matrix([q.approxlog for q in Q])
    Qdual = logQ.T * (logQ * logQ.T).inverse()
    
    Qcoeffs = [0]*len(Q)
    Qcoeffs2 = [0]*len(Q)
    
    for gloop in range(gloops):
      pubkey,seckey = cryptosystem.keygen(d,params)
      g = seckey[0]

      logg = ideal.approxlog(d,K(g,1))
      for i in range(len(Q)):
	u = Q[i]
        uu = sum(u.approxlog[j]^2 for j in range(N))
        ulogg = sum(u.approxlog[j]*logg[j] for j in range(N))
	Qcoeffs[i] += abs(ulogg * RR(1.0) / uu)

      glog = [RR(log(abs(sum(g.c[i]*sqrtD[i]*(-1)^(sum((i&j)>>t for t in range(n))) for i in range(N))))) for j in range(N)]
      glogaverage = sum(glog) / N
      glogprojected = [glog[i] - glogaverage for i in range(N)]
      eQ = vector(glogprojected) * Qdual
      for i in range(len(Q)):
        Qcoeffs2[i] += abs(eQ[i])

    for i in range(len(Q)):
      Qcoeffs[i] /= gloops
      Qcoeffs2[i] /= gloops
      ratio = Qcoeffs2[i] / Qcoeffs[i]
      assert 0.999 < ratio
      assert ratio < 1.001
    
    for i in range(len(Q)):
      prediction = sqrt(2/RR(pi))*1.11072/(sqrt(RR(N))*Q[i].approxlog[0])
      print n,D[i+1],Qcoeffs[i],prediction

    sys.stdout.flush()

