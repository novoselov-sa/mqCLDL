import goodprime
import goodprimecheat
import centermod
import mult

def evaluate(n,N,q,s,f):
  if n == 0: return f[0]
  n -= 1
  N /= 2
  g = [(f[j] + s[n] * f[j + N]) % q for j in range(N)]
  return evaluate(n,N,q,s,g)
  
def symbol2(d,f,denom,t):
  n = len(d)
  N = 2^n
  q,s,sinv = goodprimecheat.sqrtprime(d,t)
  v = evaluate(n,N,q,s,f)
  return v.jacobi(q) / denom.jacobi(q)

def symbols2(d,f,denom,low,high):
  n = len(d)
  N = 2^n
  if high < low + 1:
    raise Exception('%s must be below %s' % (low,high))
  if high == low + 1:
    return (symbol2(d,f,denom,low),)
  mid = low + 1
  while mid + (mid - low) < high:
    mid = mid + (mid - low)
  q0 = goodprimecheat.sqrtprime_product(d,low,mid)
  q1 = goodprimecheat.sqrtprime_product(d,mid,high)
  f0 = centermod.vector(f,q0)
  denom0 = denom % q0
  f1 = centermod.vector(f,q1)
  denom1 = denom % q1
  return symbols2(d,f0,denom0,low,mid) + symbols2(d,f1,denom1,mid,high)

def symbol(d,f,denom,t):
  n = len(d)
  N = 2^n
  q,s,sinv = goodprime.sqrtprime(d,t)
  v = evaluate(n,N,q,s,f)
  #print "symbol",t, v, q, denom
  return v.jacobi(q) / denom.jacobi(q)

def symbols(d,f,denom,low,high):
  n = len(d)
  N = 2^n
  if high < low + 1:
    raise Exception('%s must be below %s' % (low,high))
  if high == low + 1:
    return (symbol(d,f,denom,low),)
  mid = low + 1
  while mid + (mid - low) < high:
    mid = mid + (mid - low)
  q0 = goodprime.sqrtprime_product(d,low,mid)
  q1 = goodprime.sqrtprime_product(d,mid,high)
  f0 = centermod.vector(f,q0)
  denom0 = denom % q0
  f1 = centermod.vector(f,q1)
  denom1 = denom % q1
  return symbols(d,f0,denom0,low,mid) + symbols(d,f1,denom1,mid,high)

def altfinal(Aq, Aqd, Aqc, Aqcd, qs, vec, vecc, denom):
  chi = []
  for i in range(len(qs)):
    denom0 = denom
    q = qs[i]
    k = GF(q)
    newq = prod([ZZ(k((Aq[v][i]/Aqd[v]))).powermod(vec[v],q) for v in range(len(vec))]) % q
    newq *= prod([ZZ(k((Aqc[v][i]/Aqcd[v]))).powermod(vecc[v],q) for v in range(len(vecc))]) 
    newq = newq % q
    while denom0 > 1:
       denom0 /= 2
       #assert mod(newq,q).is_square()
       newq = mod(newq, q).sqrt()
    chi += [ZZ(newq).jacobi(q)]
  return tuple(chi)

def altfinal2(Aq, Aqd, qs, vec, denom):
  chi = []
  for i in range(len(qs)):
    denom0 = denom
    q = qs[i]
    k = GF(q)
    newq = prod([ZZ(k((Aq[v][i]/Aqd[v]))).powermod(vec[v],q) for v in range(len(vec))]) % q
    while denom0 > 1:
      denom0 /= 2
      #assert mod(newq,q).is_square()
      newq = mod(newq, q).sqrt()
    chi += [ZZ(newq).jacobi(q)]
  return tuple(chi)

def altsymbol(d, f, denom, t):
  n = len(d)
  N = 2^n
  q,s,sinv = goodprime.sqrtprime(d,t)
  v = evaluate(n,N,q,s,f)
  return v, q

def altsymbols(d, f, denom, low, high):
  n = len(d)
  N = 2^n
  if high < low + 1:
    raise Exception('%s must be below %s' % (low,high))
  if high == low + 1:
    return (altsymbol(d,f,denom,low),)
  mid = low + 1
  while mid + (mid - low) < high:
    mid = mid + (mid - low)
  q0 = goodprime.sqrtprime_product(d,low,mid)
  q1 = goodprime.sqrtprime_product(d,mid,high)
	#
	# centermod is defined in centermod.pyx
	# centermod.vector(v,q) takes the coos of vector v mod q in [-q/2, q/2) (not sure about the end points)
	#
  f0 = centermod.vector(f,q0)
  denom0 = denom % q0
  f1 = centermod.vector(f,q1)
  denom1 = denom % q1
  return altsymbols(d,f0,denom0,low,mid) + altsymbols(d,f1,denom1,mid,high)
