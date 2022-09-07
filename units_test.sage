proof.number_field(False)

from units import *
import ring

def evaluate(n,N,s,f):
  if n == 0: return f[0]
  n -= 1
  N /= 2
  return evaluate(n,N,s,[f[j] + s[n] * f[j + N] for j in range(N)])

for n in range(1,5):
  print(n)
  N = 2^n

  for loop in range(10):
    if loop < 5:
      while True:
        d = tuple(uniq(sorted([random_prime(30) for i in range(n)])))
        if len(d) == n:
          break
    else:
      while True:
        d = tuple(ZZ.random_element(-50,50) for i in range(n))
        if not ring.has_any_nontrivial_square_product(d):
          break

    print(d)

    w,t = torsion(d)
    check1 = t^w; check2 = units(d)(1)
    assert check1.element == check2.element
    assert check1.powers == check2.powers
    assert check1.approxlog == check2.approxlog
    for j in range(1,w):
      assert (t^j).element != check2.element

    U = generators_mod_torsion(d)

    inputs = len(U)
    for outputs in range(5):
      e = [[randrange(2) for j in range(inputs)] for i in range(outputs)]
      h = subsetprod(d,U,e)
      assert len(h) == outputs
      for hi,ei in zip(h,e):
        check = units(d)(1)
        for j in range(inputs):
          if ei[j]:
            check *= U[j]
        assert check.element == hi.element
        assert check.powers == hi.powers
        assert check.approxlog == hi.approxlog

    inputs = len(U)
    for outputs in range(5):
      e = [[randrange(-5,6) for j in range(inputs)] for i in range(outputs)]
      h = powerprod(d,U,e)
      assert len(h) == outputs
      for hi,ei in zip(h,e):
        check = units(d)(1)
        for j in range(inputs):
          check *= U[j]^ei[j]
        assert check.element == hi.element
        assert check.powers == hi.powers
        assert check.approxlog == hi.approxlog

    Q = mqgenerators_mod_torsion(d)
    assert len(Q) == len(U)

    regratio = approxregulator(d) * mqindex(d) / approxmqregulator(d)
    assert regratio >= 0.999 and regratio < 1.001

    M = matrix(QQ,[compresspowers(d,u.powers) for u in U])
    M = M.inverse()
    for i in range(len(M.rows())):
      check1 = powerprod(d,U,M.rows()[i:i+1])[0]^w
      check2 = Q[i]^w
      assert check1.element == check2.element
      assert check1.powers == check2.powers
      assert check1.approxlog == check2.approxlog

    x = var('x')
    K = NumberField([x^2-di for di in d],'a')
    U = tuple(evaluate(n,N,K.gens(),u.element.numer.c)/u.element.denom for u in U)
    assert w == len(K.roots_of_unity())
    F = K.unit_group().fundamental_units()
    for u in list(U) + list(F):
      assert u.norm() in [-1,1]
    phi = K.embeddings(ComplexField(1000))
    logU = matrix([[log(abs(phij(x))) for phij in phi] for x in U])
    logF = matrix([[log(abs(phij(x))) for phij in phi] for x in F])
    dual = logU.T * (logU * logU.T).inverse()
    M = logF * dual
    M = [[x.round() for x in row] for row in M.rows()]
    assert det(matrix(M)) in [-1,1]
    for Fi,Mi in zip(F,M):
      Ficheck = prod(Ui**ei for Ui,ei in zip(U,Mi))
      assert Ficheck^w == Fi^w

    regratio2 = approxregulator(d) / K.regulator()
    assert regratio2 >= 0.999 and regratio2 < 1.001
