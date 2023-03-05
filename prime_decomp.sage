import nprofile
profiling = ['prime_decomp']

import field
import polynomial_ring
from memoized import memoized

@nprofile.profile
@memoized
def prime_decomp_old(d, p):
  K = field.field(d)
  assert Mod(K.idx(), p) != 0, f"Unimplemented case: p = {p} divides [O_K:Z[theta_K]] = {K.idx()}"
  theta = K.abs_gen()
  R = polynomial_ring.polynomial_ring(d)
  pol = K.absolute_polynomial()
  fp = pol.to_sage().change_ring(GF(p)).factor()
  res = []
  for F in fp:
    t = R.from_sage(F[0])
    e = F[1]
    I = (p, t.evaluate(theta))
    res.append([I, e])
  return res

@nprofile.profile
@memoized
def prime_decomp(d, p):
  K = field.field(d)
  if Mod(K.idx(), p) != 0:
    return prime_decomp_fast(d, p)
  else:
    return prime_decomp_slow(d, p)

def ZZ_to_aut(v, n):
  mu = list(reversed(ZZ(v).bits()))
  mu = [0]*(n-len(mu)) + mu # append zeroes, e.g. [1] => [0,0,0,1] for n=4
  return mu

# Prime decomposition using Sage implementation.
@nprofile.profile
@memoized
def prime_decomp_slow(d, p):
  K = field.field(d)
  theta = K.abs_gen()
  R = polynomial_ring.polynomial_ring(d)
  pol = K.absolute_polynomial().to_sage().change_ring(QQ)
  # Since factoring in relative field is slow, we use absolute field.
  K_sage = NumberField(pol, names="a")
  res = []
  for I,e in K_sage.factor(p):
    assert len(I.gens()) == 2
    assert I.gens()[0] == p
    alpha = I.gens()[1].polynomial()
    t = R.from_sage(alpha)
    I = (p, t.evaluate(theta), [[], "other"])
    res.append([I, e])
  return res

@nprofile.profile
@memoized
def prime_decomp_fast(d, p):
  n = len(d)
  N = 2^n
  K = field.field(d)
  theta = K.abs_gen()

  assert Mod(K.idx(), p) != 0, f"Unimplemented case: p = {p} divides [O_K:Z[theta_K]] = {K.idx()}"
  #assert Mod(K.discriminant(), p) != 0, f"Unimplemented case: p = {p} divides disc_K = {K.discriminant()}"

  #t1 = [for di in d if Mod(di,p).is_square()]
  #t2 = [for di in d if not Mod(di,p).is_square()]
  t1 = [K.gens()[i] for i in range(n) if Mod(d[i],p).is_square()]
  t2 = [K.gens()[i] for i in range(n) if not Mod(d[i],p).is_square()]
  theta1 = sum(t1, K.zero())
  theta2 = sum(t2, K.zero())

  res = []
  # linear case, p splits completely
  if theta2 == K.zero():
    for i in range(0,N):
      mu = ZZ_to_aut(i, n)
      alpha = K.abs_gen().apply_aut(mu)
      #I = [p, K.abs_gen() + lift(alpha.evaluate_mod(p)), [mu[-1], "linear"]]
      #I = [p, K.abs_gen() + lift(alpha.evaluate_mod(p)), [mu, "linear"]]
      I = [p, K.abs_gen() + lift(alpha.evaluate_mod_ext(p)), [mu, "linear"]]
      #res.append([I, 1])
      found = False
      for i in range(len(res)):
        J,e = res[i]
        assert p == J[0]
        if J[1] == I[1]:
          res[i][1] = e + 1
          found = True
          break
      if not found:
        res.append([I, 1])
  else:
    # quadratic case, p splits into deg(K)/2 prime ideals
    N1 = 2^len(t1)
    N2 = 2^(len(t2)-1)
    for i1 in range(N1):
      for i2 in range(N2):
        gamma1 = ZZ_to_aut(i1, len(t1))
        gamma2 = ZZ_to_aut(i2, len(t2))
        gamma = []
        for i in range(n):
          gi = K.gens()[i]
          try:
            j = t1.index(gi)
            gamma.append(gamma1[j])
          except ValueError:
            j = t2.index(gi)
            gamma.append(gamma2[j])
        # the following lines can be faster
        #theta1_a = sum([K.from_ZZ((-1)^gamma1[i])*t1[i] for i in range(len(t1))], K.zero()) #theta1.apply_aut(gamma)
        #theta2_a = sum([K.from_ZZ((-1)^gamma2[i])*t2[i] for i in range(len(t2))], K.zero()) #theta2.apply_aut(gamma)
        theta1_a = theta1.apply_aut(gamma)
        theta2_a = theta2.apply_aut(gamma)

        theta1_a_p = theta1_a.evaluate_mod_ext(p)
        theta2_a_sq = lift((theta2_a^2).evaluate_mod_ext(p))

        g = theta^2 - theta * 2 * lift(theta1_a_p) + lift(theta1_a_p^2) - theta2_a_sq

        #res.append([[p, g, [gamma, "quadratic"]], 1])
        I = [p, g, [gamma, "quadratic"]]
        found = False
        for i in range(len(res)):
          J,e = res[i]
          assert p == J[0]
          if J[1] == g:
            res[i][1] = e + 1
            found = True
            break
        if not found:
          res.append([I, 1])
  return res
