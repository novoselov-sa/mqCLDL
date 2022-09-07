# *****************************************************************************
# Univariate polynomials over the field Q[sqrt(d1), ..., sqrt(dn)].
# *****************************************************************************

import nprofile
profiling = ['minimal_polynomial', 'mult_sb']

import field
import numpy as np
from memoized import memoized

def trim_zeroes(arr):
  i = len(arr)-1
  while i > -1 and (not arr[i].is_nonzero()):
    i -= 1
  if i == -1:
    return []
  else:
    return arr[0:i+1]

@memoized
def polynomial_ring(d, var = "x"):
  n = len(d)
  N = 2^n
  var = var
  K = field.field(d)
  class PR:
    def __init__(pol, coeffs):
      pol.coeffs = trim_zeroes(coeffs)
    def deg(pol):
      return len(pol.coeffs)-1
    def __eq__(pol1,pol2):
      if pol1.__class__ != pol2.__class__: raise Exception('do not know how to compare %s to %s' % (pol1,pol2))
      return pol1.coeffs == pol2.coeffs
    def __add__(pol1, pol2):
      c = [pol1.coeffs[i] + pol2.coeffs[i] for i in range(min(len(pol1.coeffs),len(pol2.coeffs)))]
      if pol1.deg() > pol2.deg():
        c += pol1.coeffs[len(pol2.coeffs):]
      else:
        c += pol2.coeffs[len(pol1.coeffs):]
      return PR(c)
    def __neg__(pol):
      return PR([-pol.coeffs[i] for i in range(len(pol.coeffs))])
    def __sub__(pol1,pol2):
      return pol1 + (-pol2)
    def __mul__(pol1, pol2):
      return PR(mult_sb(d, pol1.coeffs, pol2.coeffs))
    def __repr__(pol):
      return 'polynomial_ring.polynomial_ring(%s,%s)(%s)' % (d, var, pol.coeffs)
    def __getitem__(pol, i):
      if i >= len(pol.coeffs) or i < 0:
        return K.zero()
      return pol.coeffs[i]
    def to_sage(pol, names="a", quotient=False):
      if pol.deg() == -1:
        return ZZ(0)
      c0 = pol[0].to_sage()
      R0 = c0.parent()
      R1 = PolynomialRing(R0, var)
      x = R1.gens()[0]
      res = R1(c0)
      for i in range(1,len(pol.coeffs)):
        res += pol.coeffs[i].to_sage(names=names, quotient=quotient) * x^i
      return res
    def evaluate(pol, val):
      if pol.deg() == -1:
        return K.zero()
      r = pol[-1]
      for i in range(pol.deg(),-1,-1):
        r = r*val + pol[i]
      return r

    @staticmethod
    def zero():
      return PR([K.zero()])
    @staticmethod
    def one():
      return PR([K.one()])
    @staticmethod
    def x():
      return PR([K.zero(), K.one()])
    @staticmethod
    def from_sage(pol):
      c = list(pol)
      return PR([ K(tuple([QQ(c[i]).numerator()] + [0]*(N-1)), QQ(c[i]).denominator()) for i in range(len(c))])
  return PR

@memoized
def absolute_polynomial(d):
  K = field.field(d)
  return K.abs_gen().minimal_polynomial()

@nprofile.profile
def mult_sb(d, pol1, pol2):
  K = field.field(d)
  if len(pol1) == 0 or len(pol2) == 0:
    return []
  m = len(pol1) + len(pol2) - 2
  c = []
  for k in range(m+1):
    a = K.zero()
    for i in range(max(0,k-len(pol2)+1), min(len(pol1)-1, k)+1):
      a += pol1[i] * pol2[k-i]
    c.append(a)
  return c

@nprofile.profile
# Minimal polynomial computation for the field element.
def minimal_polynomial(d, f):
  n = len(d)
  N = 2^n
  K = field.field(d)
  PR = polynomial_ring(d)
  res = PR([-f, K.one()])
  for i in range(1,N):
    mu = list(reversed(ZZ(i).bits()))
    mu = [0]*(n-len(mu)) + mu
    f2 = f.apply_aut(mu)
    if f2 != f:
      res *= PR([-f2, K.one()])
  return res