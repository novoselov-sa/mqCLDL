import ideal
import relations
import fb
import trees

#food = {'a': 17, 'b': 29, 'c': 37, 'd': 41}
#food = {'a':7, 'b': 11, 'c': 19, 'd': 23}
#food = {'a':7, 'b': 79, 'c': 223, 'd': 229}
#food = {'a': 1297, 'b': 4759, 'c': 7057, 'd': 8761}
#food = {'a': 4759, 'b': 8761} # h_4759 = 13, h_8761 = 27
#food = {'a': 1297, 'b': 7057, 'c': 8761} # h_1297 = 11, h_4759 = 13, h_7057 = 21
#food = {"a": 5, "b": 13, "c": 17}
#food = {"a": 5, "b": 13, "c": 17, "d": 29}
#food = {"a": 5, "b": 13, "c": 17, "d": 29, "e": 37}
#food = {"a": 5, "b": 13, "c": 17, "d": 29, "e": 37, "g": 41}
#food = {"a": 5, "b": 13, "c": 17, "d": 29, "e": 37, "g": 41, "h": 53}
#food = {'a': -7, 'b': -11, 'c': -19}
#food = {'a': -7, 'b': -11, 'c': -19, 'd': -23}
#food = {'a':-11, 'b':-19, 'c':-23, 'd':-31} #Cl_K: C3 x C3 x C12 x C2520
#food = {'a': -11, 'b': -19, 'c': -23, 'd': -31, 'e': -43} #Cl_K: C2 x C2 x C4 x C4 x C12 x C48 x C96 x C3360 x C3360 x C443520
#food = {'a': 5, 'b': 79, 'c': 89, 'd': 97, 'e': 101}  #result: C2 x C2 x C2 x C2 x C2 x C2 x C2 x C2 x C4 x C4 x C4 x C24 x C48 x C48 x C96 x C96
#food = {'a': -7, 'c': -47, 'b': -31} # |Cl_K| = 510
#food = {'a': 5, 'b': 79, 'c': 89}
#food = {"a": -3, "b": -7, "c": -11, "d": -19}

food = trees.get_food() # load food from trees

d = tuple(sorted(food.values(), key=lambda di: abs(di)))
#d = [di for di in Primes()[2:100] if Mod(di, 4) == 1][0:8]

VERIFY_LIGHT = False # some light checks and checks based on class group computation for subfields with n = len(d) <= 5
VERIFY = False # checks based on class group computation for all subfields
VERIFY_HEAVY = False
proof.all(False)

pari.allocatemem(1024*1024*1024)
print(f"Pari stack size: {pari.stacksize()}")

def random_ideal(d, bound=10^6):
    while True:
        q = ZZ.random_element(3, bound)
        if Mod(q, 2) == 0:
            continue
        try:
            qs = [ZZ(sqrt(Mod(di, q))) for di in d]
        except:
            continue
        return ideal.ideals(d)(q, qs)

def is_smooth(b, FB):
    r = b
    for i in range(len(FB)):
        if Mod(r, FB[i].prime) == 0:
            r = r / (FB[i].prime ^ valuation(r, FB[i].prime))
            if r == 1:
                return True
    return False

def factor_over_FB(K, I, FB):
    N_I = pari.idealnorm(K.pari_nf(), I)
    e = []
    for i in range(len(FB)):
        if Mod(N_I, FB[i].prime) != 0:
            e.append(0)
        else:
            P = K.ideal(FB[i].prime, K(FB[i].elts))
            v = ZZ(pari.idealval(K.pari_nf(), I, P.pari_prime()))
            e.append(v)
    # verification
    #S = [K.ideal(i.prime, K(i.elts)) for i in FB]
    #print(K.ideal(I).factor())
    #for pr,e2 in K.ideal(I).factor():
    #    assert(pr in S)
    #    assert(e2 == e[S.index(pr)])
    return e

def strip_ones(a):
    for i in range(len(a)):
        if a[i] != 1:
            return a[i:]
    return []

def apply_aut_rec(K, K_gens, f, mu, i):
    n = len(K_gens)
    if i >= n:
        return f # f in QQ
    res = K.zero()
    g = K(K_gens[i])
    c = f.list()
    res = apply_aut_rec(K, K_gens, c[0], mu, i+1) + apply_aut_rec(K, K_gens, c[1], mu, i+1) * (-1)^mu[i] * g
    return res

def apply_aut(K, f, mu):
    K_gens = K.gens()
    return apply_aut_rec(K, K_gens, f, mu, 0)

def ideal_to_hnf(K, d, I):
    g = K.gens()
    for i in range(len(g)):
        assert g[i]^2 == d[i], f"{g[i]^2} != {d[i]}"
    # basis of the ring R = Z[sqrt(d_1), ..., sqrt(d_n)]
    F = [1]
    for gi in g:
        F += [gi*hi for hi in F]
    
    I_gens = [I.q]
    for i in range(len(d)):
        I_gens.append(g[i] - I.s[i])

    M = matrix()
    for i in range(len(I_gens)):
        for j in range(len(F)):
            r = matrix(pari.nfalgtobasis(K.pari_nf(), F[j]*I_gens[i]))
            if M.nrows() == 0:
                M = r
            else:
                M = M.stack(r)
    #print(f"M = {M}")
    return M.hermite_form(include_zero_rows=False)

def bbp_ideal_to_sage(d, I, K = None):
    if K == None:
        K = NumberField([x^2 - d[i] for i in range(len(d))], names=list(sorted(food.keys()))[:len(d)])
    I_hnf = ideal_to_hnf(K, d, I)
    I_pari = pari.idealhnf(K.pari_nf(), I_hnf.hermite_form(include_zero_rows=False).transpose())
    I_sage = K.ideal(I_pari)
    assert I_sage.absolute_norm() == I.q
    return I_sage

# Converting prod_i P_i^e_i => prod_i g_i^e'_i.
# Here P_1, P_2, ... are primes from the factor base and g_1, g_2, ... are class group generators.
# V is a matrix from SNF of relations matrix A, i.e. A = U B V.
def prod_primes_to_gens(e, V, B = None, reduce = True, strip_zeroes=False):
    e_g = [sum(int(V[i,j]*e[i]) for i in range(V.ncols())) for j in range(V.nrows())]
    if reduce and B != None:
        e_g = [int(lift(Mod(e_g[i], B[i,i]))) for i in range(len(e_g))]
    if strip_zeroes:
        for i in range(len(e_g)):
            if B[i,i] != 1:
                e_g = e_g[i:]
                break
    return e_g

# Converting prod_i prod_i g_i^e'_i => P_i^e_i.
# Here P_1, P_2, ... are primes from the factor base and g_1, g_2, ... are class group generators.
# V_inv = V^(-1) is a inverse of matrix V from SNF of relations matrix A, i.e. A = U B V.
def prod_gens_to_primes(e_g, V_inv):
    e = [sum(int(V_inv[i,j]*e_g[i]) for i in range(V_inv.ncols())) for j in range(V_inv.nrows())]
    return e

# We assume that the factor base for each subfield is precomputed using the script "trees_generation.sage"
# and it is stored in the "trees" folder.
def cl_dlog_quad(D, I, d_parent=[]):
    file = relations.convert_letters(d_parent, seed=1, food=food)
    FB = fb.load_primes(tuple([D]), file, food=food)
    M = load_relations([D])
    B,U,V=M.smith_form()
    V_inv = V^(-1)
    B_t = []
    for i in range(B.ncols()):
        #if B[i,i] != 1:
        B_t.append(B[i,i])
    #print(f"B_t = {strip_ones(B_t)}")

    print(f"Computing DLOG for the ideal I = {I} in group {strip_ones(B_t)} ...")

    K = QuadraticField(D, name="a")
    h = K.class_number(proof=False)
    #print(f"-> h_K = {h}")
    if h == 1:
        # This is a trivial case when all ideals lie in the same class.
        e = [0 for i in range(len(FB))]
        e_g = prod_primes_to_gens(e, V, B, strip_zeroes=True)
        #print(f"-> e = {e} (trivial case)")
        print(f"[done] Computing DLOG for the ideal I = {I} in group {strip_ones(B_t)} ... {e_g}\n")
        return e
    #print("-> FB:", [[K.ideal(FB[i].prime, K(FB[i].elts)), K.ideal(FB[i].prime, K(FB[i].elts)).is_principal(), K.ideal(FB[i].prime, K(FB[i].elts)).gens()] for i in range(len(FB))] )
    #print("-> FB (primes):", set([FB[i].prime for i in range(len(FB))]))
    
    # transformation of the input ideal I to the O_K-basis used by SageMath
    if Mod(D, 4) == 1:
        # check that Z-basis of O_K is (1, (sqrt(D)-1)/2))
        assert(K.pari_zk()[0] == 1 and K.pari_zk()[1] == K.pari_zk()[1].variable()/2 - 1/2)
        # transform the ideal of Z[sqrt(D)] to the ideal of O_K with respect to Pari/GP basis by putting sqrt(D) = 2*((sqrt(D)-1)/2) + 1 = (1, 2)
        I2_hnf = matrix([[I.q, 0], [2*I.q, I.q], [-I.s[0]+1, 2], [-I.s[0]+D, -2*I.s[0]]])
        # Sage uses row vectors in HNF, while Pari uses column vectors. So, we need to transpose the matrix.
        I2 = pari.idealhnf(K.pari_nf(), I2_hnf.hermite_form(include_zero_rows=False).transpose())
        N_I2 = pari.idealnorm(K.pari_nf(), I2)
        I2 = K.ideal(I2)
        assert(N_I2 == I2.absolute_norm())
        I2_alt = bbp_ideal_to_sage([D], I, K)
        assert I2_alt == I2, f"{I2_alt} != {I2}"
        # Sage and Pari use different bases for O_K and this could cause the problems.
        # The following code checks that everything is correctly transformed from Pari/GP to Sage.
        P_test = K.factor(ZZ.random_element(2^10, 2^16))[0][0]
        P_hnf = matrix(P_test.pari_hnf())
        P_hnf = pari.idealhnf(K.pari_nf(), P_hnf)
        assert K.ideal(P_hnf) == P_test, f"{K.ideal(P_hnf).pari_hnf()} != {P_test.pari_hnf()}"
    else:
        assert(K.ring_of_integers().ring_generators()[0] == K.gen())
        I2 = K.ideal([I.q, K.gen() - I.s[0]])
    #print("-> I converted, Z[sqrt(D)] => O_K: I2 =", I2.gens())
    #print("--> is_principal [pari]?", pari.bnfisprincipal(K.pari_bnf(), I2))
    #print("--> is_principal [sage]?", I2.is_principal(proof=False))

    if I2.is_principal(proof=False):
        e = [0 for i in range(len(FB))]
        e_g = prod_primes_to_gens(e, V, B, strip_zeroes=True)
        #print(f"-> e = {e} (trivial case 2)")
        #print(f"--> I2.ideal_class_log() = {I2.ideal_class_log()}")
        #assert set([ZZ(i) for i in I2.ideal_class_log()]) == set([0]), f"{set(I2.ideal_class_log())} != {set([0])}"
        print(f"[done] Computing DLOG for the ideal I = {I} in group {strip_ones(B_t)} ... {e_g}")
        return e

    # simple correctness check, norms should be the same
    N_I2 = I2.absolute_norm()
    #print("--> N(I2) =", N_I2)
    assert(N_I2 == I.q)
 
    # We look for a random ideal J s.t. LLL(I*J) factors completely over the factor base.
    while True:
        # Generation of random smooth ideal J.
        e1 = [ZZ.random_element(h) for i in range(len(FB))] # taking h_K seems to be wrong, since we have first write this ideal in terms of generators of class group
        #e1 = [ZZ.random_element(1000) for i in range(len(FB))]
        J = prod(K.ideal(FB[i].prime, K(FB[i].elts))^e1[i] for i in range(len(FB)))
        #print("-> J =", J.gens())
        #print("-> J = (...)")
        #print(f"--> e1 = {e1}")
        assert ZZ(pari.idealnorm(K.pari_nf(), J)) != 1

        I3 = I2 * J
        I3_red = pari.idealred(K.pari_nf(), I3.pari_hnf())
        N_I3_red = pari.idealnorm(K.pari_nf(), I3_red)
        #print("-> LLL(I2 * J) =", I3_red)
        #print(f"-> N(LLL(I2 * J)) = {N_I3_red} = {QQ(N_I3_red).factor()}")
        if N_I3_red == 1:
            #print(f"-> trivial reduction, e = e1 = {e1}")
            e = e1
            break
        smooth = is_smooth(N_I3_red, FB)
        if smooth:
            e2 = factor_over_FB(K, I3_red, FB)
            #print(f"-> smooth case: e = e1 - e2, e2 = {e2}")
            e = [e1[i] - e2[i] for i in range(len(FB))]
            break
    #print(f"-> e = {e}")
    # verification
    if VERIFY or VERIFY_LIGHT:
        Jt = prod(K.ideal(FB[i].prime, K(FB[i].elts))^e[i] for i in range(len(FB)))
        IJ = K.ideal(I2) * Jt
        CL_K = K.class_group(proof=False)
        assert CL_K(IJ) == CL_K.identity()
        # the following alternative is very slow even for small fields
        #assert IJ.is_principal(proof=False), f"Jt = {Jt.gens()}, I2 * Jt = {IJ.gens()}"
    # reduce vector modulo orders of class generators
    e_g = prod_primes_to_gens(e, V, B)
    #print(f"-> e_g = {e_g}")
    #print(f"--> I2.ideal_class_log() = {I2.ideal_class_log()}")
    
    
    #if set([ZZ(i) for i in I2.ideal_class_log()] + [0]) == set(e_g + [0]):
    #    print(f"-> {set(I2.ideal_class_log())} != {set(e_g)} [error possible]") 

    # check that we have correct matrices
    assert prod_gens_to_primes(prod_primes_to_gens(e, V, B, reduce=False), V_inv) == e
    e = prod_gens_to_primes(e_g, V_inv)
    #print(f"-> e (reduced) = {e}")
    print(f"[done] Computing DLOG for the ideal I = {I} in group {strip_ones(B_t)} ... {prod_primes_to_gens(e, V, B, strip_zeroes=True)}\n")
    return e

def prod_ideal(K, FB, e):
    assert len(FB) == len(e), f"Wrong length of factor base or input vector, {len(FB)} != {len(e)}"
    return prod(K.ideal(FB[i].prime, K(FB[i].elts))^int(e[i]) for i in range(len(FB)))

def load_relations(d):
    return relations.load_matrix(d)

# Solving the discrete logarithm problem for an element in cyclic group of 2-power order.
# Input: 
# 1) d, FB are the field and factor base
# 2) an element h = g_i^(ei*t)
# 3) the group <g_i^t>, where g_i is a generator of cyclic subgroup of the class group. 
# 4) r is s.t. #<g_i^t> = 2^r
# Output: l s.t. (g_i^t)^l = g_i^(ei*t)
# Note: in order to avoid evaluation of prime ideals products we work with exponents only.
def dlog_pow_2(ei, t, r, bi):
    #print(f"[begin] dlog_pow_2(ei = {ei}, t = {t}, r = {r}, bi = {bi})")
    l = 0
    gamma_e = 2^(r-1)*t # gamma = generator^(2^(r-1)) = (g_i^t)^(2^(r-1))
    for k in range(r):
        h_k_e = (-l*t + ei*t)*2^(r-1-k) # h_k = (generator^-l * h)^(2^(r-1-k)) = g_i^((-l*t + ei*t)*2^(r-1-k))
        # We have to find d_k in {0,1} s.t.
        # gamma^d_k = h_k, i.e. g_i^(d_k*t*2^(r-1)) = g_i^((-l*t + ei*t)*2^(r-1-k)) or
        # d_k*t*2^(r-1) = (-l*t + ei*t)*2^(r-1-k) mod b_i, where b_i = #<g_i>
        if Mod(h_k_e, bi) == 0:
            d_k = 0
        elif Mod(gamma_e, bi) == Mod(h_k_e, bi):
            d_k = 1
        else:
            raise Exception("DLOG_pow2: Something wrong in assumptions!")
        #print(f"-> dlog_pow_2: k = {k}, d_k = {d_k}, ei = {ei} = {Mod(ei, 2^r)} mod 2^r, {Mod(h_k_e, bi)}, {Mod(gamma_e, bi)}")
        l = l + 2^k * d_k
    #print(f"[done] dlog_pow_2({ei}, {t}, {r}, {bi}) => l = {l}")
    return l

def cyc_sqrt(e_i, b_i):
    if b_i == 0:
        return [0]
    # Generator of order 2 group is always non square
    if b_i == 2 and Mod(e_i, 2) == 1:
        print(f"[cyc_sqrt] sqrt doesn't exists for e_i = {e_i}, b_i = {b_i}")
        return None
    if b_i == 2 and Mod(e_i, 2) == 0:
        return [0, 1]
    if Mod(b_i, 2) == 1:
        return [ZZ(lift(Mod((b_i+1)/2 * e_i, b_i)))]
    r = valuation(b_i, 2)
    t = b_i / 2^r
    l = dlog_pow_2(e_i, t, r, b_i)
    # In the group of order 2^r we can have l = 1 | 2^r. If l != 1 and Mod(l,2) != 0 then look for an error in DLOG computation.
    if Mod(l, 2) != 0:
        print(f"[cyc_sqrt] sqrt doesn't exists for e_i = {e_i}, b_i = {b_i}")
        return None
    #assert Mod(l, 2) == 0, f"Result of DLOG_pow2 should be even (l = {l}), otherwise sqrt(J) does not exist!"
    sr_1 = e_i * (t+1) / 2 - t*l / 2
    sr_1 = ZZ(lift(Mod(sr_1, b_i)))

    sr_2 = sr_1 + b_i / 2
    sr_2 = ZZ(lift(Mod(sr_2, b_i)))
    return [sr_1, sr_2]

# Computation of square root (in class group) of smooth ideal class J = prod P[i]^e[i] for P[i] in the factor base.
# Requires precomputed relations (run testrelations.sage with the same 'food').
def smooth_ideal_sqrt(d, e, d_parent=()):
    #print(f"[begin] smooth_ideal_sqrt(d = {d})")
    #print(f"-> e = {e}")
    m = len(e)

    M = load_relations(d)
    #print("\n\nM:", M)
    B,U,V=M.smith_form()
    Ut = U^-1
    V_inv = V^-1
    #print("-> B:", strip_ones([B[i,i] for i in range(B.ncols())]))

    # Verifying correctnes of generators.
    # TODO: Find out why we should take the matrix V^-1 instead U^-1.
    # In [Buchmann-Düllmann, p.136, (4)] it is proposed to take U^-1 for transformation, but the code works only with V^(-1).
    # In [Biasse-Vredendaal, p.107, 3A] it is proposed to use V.
    if VERIFY_HEAVY or VERIFY or (VERIFY_LIGHT and len(d) <= 5):
        print("-> Verifying generators of class group ... ")
        if d_parent == []:
            file = relations.convert_letters(tuple(d), seed=1, food=food)
        else:
            file = relations.convert_letters(tuple(d_parent), seed=1, food=food)
        FB = fb.load_primes(tuple(d), file, food=food)
        print("--> Init field K ... ")
        K = NumberField([x^2 - d[i] for i in range(len(d))], names=list(sorted(food.keys()))[:len(d)])
        print("--> Computing Cl_K ... ")
        Cl_K = K.class_group(proof=False)
        print("---> ", Cl_K.elementary_divisors())
        B_t = []
        for i in range(B.ncols()):
            if B[i,i] != 1:
                B_t.append(B[i,i])
        if sorted(B_t) != sorted(Cl_K.elementary_divisors()):
            print(f"----> Group structure mismatch! We are working in subgroup {B_t} of class group instead of full group {sorted(Cl_K.elementary_divisors())}!")
            #open("cl_dlog_warns.log", "a").write(f"d = {d}, CL_K (sage) = {Cl_K.elementary_divisors()}, CL_K (BV) = {B_t}")
        if VERIFY_HEAVY:
            print("--> Computing explicit class group generators ... ")
            g = []
            for i in range(B.ncols()):
                gi_v = [V_inv[i,j] for j in range(V_inv.ncols())]
                gi = prod_ideal(K, FB, gi_v)
                #print(f"i = {i} / {B.ncols()}, gi_v = {gi_v}, order = {Cl_K(gi).order()}, B[i,i] = {B[i,i]}")
                assert Cl_K(gi).order() == B[i,i], f"Wrong generators, {Cl_K(gi).order()} != {B[i,i]}"
                g.append(gi)
                print(".", end="")
        print("--> [ok]")

    # Write the input ideal J in terms of class group generators
    #e_g = [sum(Ut[i,j]*e[i] for i in range(m)) for j in range(B.ncols())]
    #e_g = [sum(int(V[i,j]*e[i]) for i in range(m)) for j in range(B.ncols())]
    e_g = prod_primes_to_gens(e, V, B, reduce=True)
    #print(f"-> Input ideal exponents in terms of Cl_K gens, e_g = {e_g}")

    # check that ideal is written correctly in terms of generators of class group
    if VERIFY_HEAVY:
        print("-> Verifying input ideal representation in terms of class group generators ...")
        I1 = prod_ideal(K, FB, e)
        I2 = K.ideal(1)
        for i in range(len(e_g)):
            I2 = I2 * (g[i]^e_g[i])
        assert Cl_K(I1) == Cl_K(I2)
        print("--> [ok]")

    # check that transformation prod g_i^et_i => prod P_i^e_i is correct
    if VERIFY_HEAVY or VERIFY or (VERIFY_LIGHT and len(d) <= 5):
        print("-> Verifying transformation prod g_i^et_i => prod P_i^e_i ...")
        e2 = [sum(int(V_inv[i,j]*e_g[i]) for i in range(len(e_g))) for j in range(m)]
        I1 = prod_ideal(K, FB, e)
        I2 = prod_ideal(K, FB, e2)
        assert Cl_K(I1) == Cl_K(I2), "Wrong conversion prod g_i^et_i => prod P_i^e_i"
        IJ = prod_ideal(K, FB, [e[i] - e2[i] for i in range(len(e))])
        assert Cl_K(IJ) == Cl_K.identity()
        print("--> [ok]")

    # check that transform e_g[i] => (e_g[i] mod b_i) works, generators level
    if VERIFY_HEAVY:
        print("-> Verifying transform e_g[i] => (e_g[i] mod b_i) ...")
        for i in range(len(e_g)):
            e_g_red = int(lift(Mod(e_g[i], B[i,i])))
            #print(f"-> {e_g[i]} => {e_g_red}")
            assert Cl_K(g[i]^e_g[i]) == Cl_K(g[i]^e_g_red), f"Wrong assumption for i = {i}, {Cl_K(g[i]^(e_g[i]))} != {Cl_K(g[i]^e_g_red)}"
        print("--> [ok]")

    # check that transform e_g[i] => (e_g[i] mod b_i) works 2, products
    if VERIFY_HEAVY:
        print("-> Verifying transform e_g[i] => (e_g[i] mod b_i) 2 ...")
        I1 = prod_ideal(K, FB, e)
        I2 = K.ideal(1)
        for i in range(len(e_g)):
            e_g_red = int(lift(Mod(e_g[i], B[i,i])))
            I2 = I2 * (g[i] ^ e_g_red)
        assert Cl_K(I1) == Cl_K(I2)
        print("--> [ok]")

    #  check that transform e_g[i] => (e_g[i] mod b_i) works 3, full transform
    if VERIFY_HEAVY or VERIFY:
        print("-> Verifying transform e_g[i] => (e_g[i] mod b_i) 3 ...")
        e2_g = []
        for i in range(len(e_g)):
            e_g_red = int(lift(Mod(e_g[i], B[i,i])))
            e2_g.append(e_g_red)
        e2 = [sum(int(V_inv[i,j]*e2_g[i]) for i in range(len(e2_g))) for j in range(m)]
        print(f"-> e_g (cleaned) = {e2_g}")
        print(f"-> e (cleaned) = {e2}")
        I1 = prod_ideal(K, FB, e)
        I2 = prod_ideal(K, FB, e2)
        #print([[I1.is_principal(), Cl_K(I1), Cl_K(I1).order()], [I2.is_principal(), Cl_K(I2), Cl_K(I2).order()]])
        assert Cl_K(I1) == Cl_K(I2), "Wrong conversion prod g_i^et_i <=> prod P_i^e_i"
        IJ = prod_ideal(K, FB, [e[i] - e2[i] for i in range(len(e))])
        assert Cl_K(IJ) == Cl_K.identity()
        print("--> [ok]")

    # TODO: remove this or use for optimization
    # Assume that coefficients of e_g are reduced modulo b_i = #<g_i>.
    # Then the input ideal is principal only if e_g consist of zeroes.
    if set([ZZ(ei) for ei in e_g]) == {0}:
        principal = True
    else:
        principal = False
    #print(f"-> The input ideal is principal: {principal}")

    J_sqrts = []
    for i in range(len(e_g)):
        egi_sqrt = cyc_sqrt(e_g[i], B[i,i])
        if egi_sqrt == None: # Square root doesn't exist.
            return None
        if len(J_sqrts) == 0:
            J_sqrts = [[egi_sqrt[i]] for i in range(len(egi_sqrt))]
        else:
            J_sqrts = [J_sqrts[i] + [egi_sqrt[j]] for i in range(len(J_sqrts)) for j in range(len(egi_sqrt))]

    e_sqrts = []
    for J_sqrt in J_sqrts:
        e_sqrt = prod_gens_to_primes(J_sqrt, V_inv)
        #print(f"-> e_sqrt = {e_sqrt}")
        if VERIFY_HEAVY or VERIFY or (VERIFY_LIGHT and len(d) <= 5):
            print("-> Verifying smooth_ideal_sqrt result ...")
            I_v = prod_ideal(K, FB, e)
            J_sqrt_v = prod_ideal(K, FB, e_sqrt)
            IJ = prod_ideal(K, FB, [e[i] - 2*e_sqrt[i] for i in range(len(e))])
            assert Cl_K(I_v) == Cl_K(J_sqrt_v)^2, "Wrong computation of sqrt(I)"
            assert Cl_K(IJ) == Cl_K.identity(), "Wrong computation of sqrt(I)"
            #assert IJ.is_principal(proof=False), "Wrong computation of sqrt(I)" # slow even for small degree fields
            print("--> [ok]")
        e_sqrts.append(e_sqrt)
    #print(f"[done] smooth_ideal_sqrt(d = {d})\n")
    return e_sqrts

# Lift product of prime ideals prod(FB[i]^e[i] for i in 0..len(FB)) from subfield to the parent field using splitting data from in FB.
def lift_e(e, FB):
    currentp = 0
    e_lift = []
    assert len(FB)==len(e), 'len(FB)!=len(e)'
    for i in range(len(FB)):
        if currentp != FB[i].prime:
            e_lift += [e[i]*j for j in FB[i].powers]
            currentp = FB[i].prime
        else:
            c = len(FB[i].powers)
            e_lift[-c:] = [e_lift[-c:][h] + e[i]*FB[i].powers[h] for h in range(c)]
    return e_lift

# returns a vector e s.t. [ I * J ] = [ I * prod_i(P[i]^e[i] ] = [1], where P[i] are primes from factor base.  
def cl_dlog(d, I, d_parent = []):
    d = tuple(d)
    if len(d) == 1:
        return [cl_dlog_quad(d[0], I, d_parent=d_parent)] # should return smooth J, s.t. [I*J] = [1]

    M = load_relations(d)
    B,U,V = M.smith_form()
    B_t = [ZZ(B[i,i]) for i in range(B.ncols())]
    print(f"Computing DLOG for ideal I = {I} in the group {strip_ones(B_t)} ... ")

    #print(f"B_t = {strip_ones(B_t)}")
    if set(B_t) == {1}:
        res = [[0]*M.ncols()]
        e_g = prod_primes_to_gens(res[0], V, B, strip_zeroes=True)
        #print(f"h_K = 1, computation of dlog is trivial, result = {res}")
        #print(f"[done] cl_dlog(d = {d}, I = {I}, d_parent = {d_parent})\n")
        print(f"[done] Computing DLOG for ideal I = {I} in the group {strip_ones(B_t)} ... {e_g}\n")
        return res

    d_s  = d[:-1]
    d_t = d[:-2] + d[-1:]
    d_st = d[:-2] + (d[-2]*d[-1],)

    I_s = ideal.ideals(d_s)(I.q, I.s[:-1])
    I_t = ideal.ideals(d_t)(I.q, I.s[:-2] + I.s[-1:])
    I_st = ideal.ideals(d_st)(I.q, I.s[:-2] + (I.s[-2]*I.s[-1],))

    # loading precomputed ramification trees with factor base
    file = relations.convert_letters(d, seed=1, food=food)
    FB_s = fb.load_primes(tuple(d_s), file, food=food)
    FB_t = fb.load_primes(tuple(d_t), file, food=food)
    FB_st = fb.load_primes(tuple(d_st), file, food=food)

    dlog_s  = cl_dlog(d_s, I_s, d_parent=d)
    dlog_t  = cl_dlog(d_t, I_t, d_parent=d)
    dlog_st = cl_dlog(d_st, I_st, d_parent=d)

    #print(f"d_s = {d_s}, dlog_s = {dlog_s}")
    #print(f"d_t = {d_t}, dlog_t = {dlog_t}")
    #print(f"d_st = {d_st}, dlog_st = {dlog_st}")

    print(f"{len(dlog_s)*len(dlog_t)*len(dlog_st)} triples to enumerate")
    
    K_st = NumberField([x^2 - d_st[i] for i in range(len(d_st))], names=list(sorted(food.keys()))[:len(d_st)])
    FB_st_sage = [K_st.ideal(FB_st[i].prime, K_st(FB_st[i].elts)) for i in range(len(FB_st))]
    FB_st_conj = [K_st.ideal(FB_st[i].prime, apply_aut(K_st, K_st(FB_st[i].elts), [0]*(len(d_st)-1) + [1])) for i in range(len(FB_st))]

    # checks using Sage
    if d_parent == []:
        file = relations.convert_letters(d, seed=1, food=food)
    else:
        file = relations.convert_letters(d_parent, seed=1, food=food)
    FB = fb.load_primes(tuple(d), file, food=food)
    K = NumberField([x^2 - d[i] for i in range(len(d))], names=list(sorted(food.keys()))[:len(d)])
    #print(f"-> class group computation for checks using Sage, d = {d} ...")
    #Cl_K = K.class_group(proof=False)
    #print(f"--> {Cl_K}")
    I_sage = bbp_ideal_to_sage(d, I, K)
    #print(f"-> Computing dlog with Sage for checks ...")
    #dlog_sage = I_sage.ideal_class_log()
    #print(f"--> I.ideal_class_log() = {dlog_sage}")

    res = []

    for e_s in dlog_s:
        for e_t in dlog_t:
            for e_st in dlog_st:
                #print("-> e_s:", e_s)
                #print("-> e_t:", e_t)
                #print("-> e_st:", e_st, "\n")

                e_s_lift = lift_e(e_s, FB_s)
                e_t_lift = lift_e(e_t, FB_t)

                e_st_2 = []
                for i in range(len(e_st)):
                    j = FB_st_sage.index(FB_st_conj[i])
                    ei_conj = e_st[j]
                    #print(f"i = {i} => j = {j}| {FB_st[i].prime}, {FB_st[i].elts} => {FB_st[j].prime}, {FB_st[j].elts}")
                    #assert i != j or K_st(FB_st[i].elts).is_rational()
                    e_st_2.append(ei_conj)
                e_st = e_st_2
                #print(f"-> e_st.conj = {e_st}")

                # We assume that the automorphism sigma is already applied to the prime ideals in trees for the subfield field K_st.
                e_st_lift = lift_e(e_st, FB_st)

                #print(f"-> e_s_lift {d_s} => {d}: {e_s_lift}")
                #print(f"-> e_t_lift {d_t} => {d}: {e_t_lift}")
                #print(f"-> e_st_lift {d_st} => {d}: {e_st_lift}")
                
                assert(len(e_s_lift) == len(e_t_lift))
                assert(len(e_t_lift) == len(e_st_lift))

                e = [e_s_lift[i] + e_t_lift[i] - e_st_lift[i] for i in range(len(e_s_lift))]
                #print(f"-> e = {e} \n")

                e_sqrts = smooth_ideal_sqrt(d, e, d_parent=d_parent)
                if e_sqrts == None:
                    #print("--> sqrt doesn't exists, skipping this branch.")
                    continue

                #print("-> results")
                #print(f"--> d = {d}")
                #print(f"--> B = {strip_ones(B_t)}")
                check_sage = False
                for i in range(len(e_sqrts)):
                    e_sqrt = e_sqrts[i]
                    #print(f"--> result {i} [primes] = {e_sqrt}")
                    e_sqrt_g = prod_primes_to_gens(e_sqrt, V, B)
                    #print(f"--> result {i} [cl_gens] = {e_sqrt_g}")

                    #print("---> trying to find correct DLOG with LLL ...")
                    II = I_sage * prod_ideal(K, FB, e_sqrt)
                    II_red = pari.idealred(K.pari_nf(), II.pari_hnf())
                    N_II_red = pari.idealnorm(K.pari_nf(), II_red)
                    #print("----> LLL(I*sqrt(J)) =", II_red)
                    #print(f"----> N(LLL(I*sqrt(J))) = {N_II_red} = {QQ(N_II_red).factor()}")
                    if N_II_red == 1:
                        #print(f"----> correct dlog found [trivial case]!")
                        return [e_sqrt]
                    else:
                        smooth = is_smooth(N_II_red, FB)
                        if smooth:
                            II_red_e = factor_over_FB(K, II_red, FB)
                            e_sqrt_correct = [ZZ(e_sqrt[i] - II_red_e[i]) for i in range(len(e_sqrt))]
                            #print(f"----> found correct dlog [smooth case]: {e_sqrt_correct}!")
                            #print(f"-----> (in terms of gens): {prod_primes_to_gens(e_sqrt_correct, V, B)}")
                            e_g = prod_primes_to_gens(e_sqrt_correct, V, B, strip_zeroes=True)
                            print(f"[done] Computing DLOG for ideal I = {I} in the group {strip_ones(B_t)} ... {e_g}\n")
                            return [e_sqrt_correct]
                    
                    # checking using Sage, we need implementation of BB+ algortihm for HNF representation of ideals for efficient calculation
                    if II.is_principal(proof=False):
                        #print(f"----> correct dlog found with sage!")
                        e_g = prod_primes_to_gens(e_sqrt, V, B, strip_zeroes=True)
                        print(f"[done] Computing DLOG for ideal I = {I} in the group {strip_ones(B_t)} ... {e_g}\n")
                        return [e_sqrt]

                    #if set([ZZ(i) for i in dlog_sage] + [ZZ(0)]) == set(e_sqrt_g + [0]):
                    #    check_sage = True # correct square root belongs to list

                # print("-> Check results using sage ...")
                # assert check_sage, f"dlog check using Sage is failed! Correct variant doesn't belong to the results, field = {d}"
                # print("--> [ok]")

                # print("-> checking result using [I]==[J^-1] ...")
                # check_sage = False
                # for e_sqrt in e_sqrts:
                #     J_sq = prod_ideal(K, FB, e_sqrt)
                #     if Cl_K(I_sage) == Cl_K(J_sq^-1):
                #         check_sage = True
                #         break
                # assert check_sage
                # print("--> [ok]")
                for e_sqrt in e_sqrts:
                    e_sqrt_ZZ = [ZZ(e_sqrt[i]) for i in range(len(e_sqrt))]
                    #print(f"[DEBUG 42] e_sqrt_ZZ = {e_sqrt_ZZ}")
                    if not (e_sqrt_ZZ in res):
                        res.append(e_sqrt_ZZ)
    #print(f"[done] cl_dlog(d = {d}, I = {I}, d_parent = {d_parent})\n")
    res_g = [prod_primes_to_gens(r, V, B, strip_zeroes=True) for r in res] 
    print(f"[done] Computing DLOG for ideal I = {I} in the group {strip_ones(B_t)} ... {res_g}\n")
    return res

trees_food = trees.get_food()
if food != None:
    assert food == trees_food, f"Run trees generation and relation computation first! Trees are generated for {trees_food} != {food}."

bound = prod(di.abs() for di in d)
I = random_ideal(d, bound=bound)

# The ideal can be fixed as in the following examples.
#I = ideal.ideals([5, 229, 257])(669289,(303651, 6337, 202003))
#I = ideal.ideals([-7, -11, -19, -23])(2557,(143, 1040, 485, 1051))
#I = ideal.ideals([5, 13, 17, 29, 37, 41])(9032339,(2001020, 3845142, 269749, 1568238, 3382698, 560462))
#I = ideal.ideals((5, 13, 17, 29, 37))(744251,(275720, 225315, 139312, 8627, 29735))
#I = ideal.ideals((5, 13, 17, 29, 37))(271879,(49480, 115130, 12857, 3796, 90974))
#I = ideal.ideals((5, 13, 17, 29, 37))(879941,(341095, 173263, 46592, 397104, 411054))
#I = ideal.ideals((-7, -11, -19))(709,(213, 46, 119))

print(f"Target ideal: {I}")
K = NumberField([x^2 - d[i] for i in range(len(d))], names=list(sorted(food.keys()))[:len(d)])
I_hnf = ideal_to_hnf(K, d, I)
print(f"Target ideal (hnf): {I_hnf}")
I_sage = K.ideal(pari.idealhnf(K.pari_nf(), I_hnf.hermite_form(include_zero_rows=False).transpose()))
print(f"I_sage.norm = {I_sage.absolute_norm()}")
assert I_sage.absolute_norm() == I.q
#print(f"I_sage.ideal_class_log() = {I_sage.ideal_class_log()}")

dls = cl_dlog(d, I)

file = relations.convert_letters(d, seed=1, food=food)
FB = fb.load_primes(tuple(d), file, food=food)
M = load_relations(d)
B,U,V=M.smith_form()
V_inv = V^(-1)

print("Final results:")
for i in range(len(dls)):
    dl = dls[i]
    #print(f"i = {i}")
    #print(f"-> e (primes) = {dl}")
    print(f"-> i = {i}, {prod_primes_to_gens(dl, V, B, strip_zeroes=True)}")

# check correctness using Sage (slow):
print("\nChecks using Sage:")
for i in range(len(dls)):
    dl = dls[i]
    #print(f"i = {i}")
    e_g = prod_primes_to_gens(dl, V, B, strip_zeroes=True)
    print(f"i = {i}, checking e_g = {e_g} ...")
    J = prod_ideal(K, FB, dl)
    CL_K = K.class_group(proof=False)
    if CL_K(I_sage) != CL_K(J^(-1)):
        print("-> wrong result")
    else:
        print("-> ok")
    #I3 = prod([CL_K.gens()[i]^cl_sage[i] for i in range(len(CL_K.gens()))])
    #assert CL_K(I_sage) == CL_K(I3)
