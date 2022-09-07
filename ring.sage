from memoized import memoized
import mult
import div
import norm
import sqrt
import char
import subsetprod
ssp = subsetprod
import powerprod

# assumes d is tuple of ZZ
# returns True if (c*any subset of d).is_square()
def has_any_twisted_square_product(c,d):
  if len(d) == 0: return c.is_square()
  if has_any_twisted_square_product(c,d[:-1]): return True
  if has_any_twisted_square_product(c*d[-1],d[:-1]): return True
  return False

# assumes d is tuple of ZZ
# returns True if (any nonempty subset of d).is_square()
def has_any_nontrivial_square_product(d):
  if len(d) == 0: return False
  if has_any_nontrivial_square_product(d[:-1]): return True
  if has_any_twisted_square_product(d[-1],d[:-1]): return True
  return False

@memoized
# assumes d is tuple of ZZ
# assumes not has_any_nontrivial_square_product(d)
def ring_without_checking(d):
  n = len(d)
  N = 2^n
  class R:
    def __init__(f,c):
      c = tuple(ZZ(cj) for cj in c)
      if len(c) != N: raise Exception('%s has wrong length' % str(c))
      f.c = c
    def __repr__(f):
      return 'ring.ring(%s)(%s)' % (d,f.c)
    def __eq__(f,g):
      if g.__class__ != f.__class__: raise Exception('do not know how to compare %s to %s' % (f,g))
      return f.c == g.c
    def __ne__(f,g):
      return not f == g
    def is_nonzero(f):
      return not all(fj == 0 for fj in f.c)
    def __add__(f,g):
      if g.__class__ != f.__class__: raise Exception('do not know how to compute %s + %s' % (f,g))
      return R(fj+gj for fj,gj in zip(f.c,g.c))
    def __neg__(f):
      return R(-fj for fj in f.c)
    def __sub__(f,g):
      if g.__class__ != f.__class__: raise Exception('do not know how to compute %s - %s' % (f,g))
      return R(fj-gj for fj,gj in zip(f.c,g.c))
    def __mul__(f,g):
      if parent(g) == ZZ:
        return R(g * fj for fj in f.c)
      if g.__class__ != f.__class__: raise Exception('do not know how to compute %s * %s' % (f,g))
      return R(mult.mult(n,N,d,f.c,g.c))
    def __getitem__(f, i):
      return f.c[i]
    def square(f):
      return R(mult.square(n,N,d,f.c))
    def sqrt(f):
      result,S = sqrt.squareroot(n,d,n,f.c)
      result = R(result)
      if result.square() != f: raise Exception('not a square')
      return result
    def divexact(f,g): # caller guarantees divisibility
      if g.__class__ != f.__class__: raise Exception('do not know how to compute %s / %s' % (f,g))
      return R(div.div(n,N,d,f.c,g.c))
    def __div__(f,g):
      result = f.divexact(g)
      assert f == g * result
      return result
    def normonestep(f,t):
      h = norm.normonestep(n,N,d,f.c,t)
      R1 = ring_without_checking(d[:t] + d[t+1:])
      return R1(h)
    def absnorm(f):
      return norm.absnorm(n,N,d,f.c)
    def symbol(f,t):
      return char.symbol(d,f.c,1,t)
    def symbols(f,low,high):
      return char.symbols(d,f.c,1,low,high)
    def quadnorms(f):
      if n == 0: return ()
      if n == 1: return (f,)
      m = 1
      while m + m < n: m = m + m
      f0 = f
      for j in range(m): f0 = f0.normonestep(0)
      f1 = f
      for j in range(n,m,-1): f1 = f1.normonestep(j-1)
      return f1.quadnorms() + f0.quadnorms()
  return R

def ring(d):
  d = tuple(ZZ(di) for di in d)
  if has_any_nontrivial_square_product(d):
    raise Exception('%s has subset with square product' % str(d))
  return ring_without_checking(d)

def subsetprod(d,g,e):
  R = ring(d)
  for gj in g:
    if gj.__class__ != R:
      raise Exception('%s not in ring(%s)',gj,d)
  n = len(d)
  N = 2^n
  h = ssp.subsetprod(n,N,d,[gj.c for gj in g],e)
  return tuple(R(hi) for hi in h)

def scaledpowerprod(d,g,denom,e,multiplier):
  R = ring(d)
  for gj in g:
    if gj.__class__ != R:
      raise Exception('%s not in ring(%s)',gj,d)
  n = len(d)
  N = 2^n
  h = powerprod.scaledpowerprod(n,N,d,[gj.c for gj in g],denom,e,multiplier)
  return tuple(R(hi) for hi in h)
