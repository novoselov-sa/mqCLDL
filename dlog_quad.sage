import fb
import clgp
import verify
import field

# We assume that the factor base for each subfield is precomputed using the script "trees_generation.sage"
# and it is stored in the "trees" folder.
# For an ideal I retuns vector e = (e_1, ..., e_r) and g from K such that I*prod P_i^e_i = g O_K,
# where P_1, ..., P_r are prime ideals in factor base.
def dlog_quad(D, I, d_parent=(), food=None, reduce=True):
    FB = clgp.get_FB((D,), d_parent=d_parent, food=food)
    M,B,U,V,U_inv,V_inv = clgp.get_matrices((D,))
    
    comp_time = walltime()

    B_t = []
    for i in range(B.ncols()):
        B_t.append(B[i,i])

    print(f"Computing DLOG for the ideal I = {I} in group {clgp.strip_oz(B_t)} ...")
    clgp.check_clgp((D,), d_parent=d_parent)

    K = QuadraticField(D, name="a")

    h = K.class_number(proof=False)
    #print(f"-> h_K = {h}")
    if h == 1:
        # This is a trivial case when all ideals lie in the same class.
        e = [0 for i in range(len(FB))]
        e_g = clgp.primes_to_gens(e, V, B, strip_zeroes=True)
        #print(f"-> e = {e} (trivial case)")
        print(f"Finished computing DLOG for the ideal I = {I} in group {clgp.strip_oz(B_t)} ... {e_g}, {walltime(comp_time)} sec.")
        g = I.generator()
        #assert I.to_sage(K) == K.ideal(g.to_sage(K))
        print(f"-> generator: {g.to_sage()}")
        print("-> trivial class group\n")
        if verify.level() >= verify.LIGHT:
            print(f"-> Sage's result: {I.to_sage(K).ideal_class_log()}")
            print()
        assert I.to_sage(K) * fb.prod_ideal(K, FB, e) == K.ideal(g.to_sage(K))
        return (g, e)

    I2 = I.to_sage(K)

    if I2.is_principal(proof=False):
        e = [0 for i in range(len(FB))]
        e_g = clgp.primes_to_gens(e, V, B, strip_zeroes=True)
        #print(f"-> e = {e} (trivial case 2)")
        #print(f"--> I2.ideal_class_log() = {I2.ideal_class_log()}")
        print(f"Finished computing DLOG for the ideal I = {I} in group {clgp.strip_oz(B_t)} ... {e_g}, {walltime(comp_time)} sec.")
        g = I.generator()
        print(f"-> generator: {g.to_sage()}")
        #assert I.to_sage(K) == K.ideal(g.to_sage(K))
        print(f"-> case: principal input ideal\n")
        if verify.level() >= verify.LIGHT:
            print(f"-> Sage's result: {I.to_sage(K).ideal_class_log()}")
            print()
        assert I.to_sage(K) * fb.prod_ideal(K, FB, e) == K.ideal(g.to_sage(K))
        return (g, e)

    # simple correctness check, norms should be the same
    N_I2 = I2.absolute_norm()
    #print("--> N(I2) =", N_I2)
    assert(N_I2 == I.q)
    # TODO: add check by selecting random elements

    res_case = ""
 
    # We look for a random ideal J s.t. LLL(I*J) factors completely over the factor base.
    while True:
        # Generation of random smooth ideal J
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
            res_case = "trivial reduction, norm = 1"
            break
        smooth = fb.is_smooth(N_I3_red, FB)
        if smooth:
            e2 = fb.factor(K, I3_red, FB)
            #print(f"-> smooth case: e = e1 - e2, e2 = {e2}"
            e = [e1[i] - e2[i] for i in range(len(FB))]
            res_case = "smooth product"
            break
    print(f"-> e = {e}")

    if reduce:
        e = clgp.reduce_NP((D,), e)
        print(f"-> e (reduced) = {e}")

    # verification
    if verify.level() >= verify.LIGHT:
        Jt = prod(K.ideal(FB[i].prime, K(FB[i].elts))^e[i] for i in range(len(FB)))
        IJ = K.ideal(I2) * Jt
        CL_K = K.class_group(proof=False)
        assert CL_K(IJ) == CL_K.identity()
        # the following alternative is very slow even for small fields
        #assert IJ.is_principal(proof=False), f"Jt = {Jt.gens()}, I2 * Jt = {IJ.gens()}"

    # reduce vector modulo orders of class generators
    #e_g = clgp.primes_to_gens(e, V, B)
    e_g = clgp.primes_to_gens(e, V, B, strip_zeroes=True)
    #print(f"-> e_g = {e_g}")

    # check that we have correct matrices
    #assert clgp.gens_to_primes(clgp.primes_to_gens(e, V, B, reduce=False), V_inv) == e
    #e = clgp.gens_to_primes(e_g, V_inv)

    Jt = prod(K.ideal(FB[i].prime, K(FB[i].elts))^e[i] for i in range(len(FB)))
    IJ = K.ideal(I2) * Jt
    g = IJ.gens_reduced()
    if len(g) > 1: raise Exception('ideal is not principal')
    g = g[0]
    g0,g1 = list(g)
    s = lcm(QQ(g0).denom(),QQ(g1).denom())
    g = field.field((D,))((ZZ(g0*s),ZZ(g1*s)),s)

    print(f"Finished computing DLOG for the ideal I = {I} in group {clgp.strip_oz(B_t)} ... {clgp.primes_to_gens(e, V, B, strip_zeroes=True)}, {walltime(comp_time)} sec.")
    print(f"-> generator: {g.to_sage()}")
    print(f"-> case: {res_case}\n")
    if verify.level() >= verify.LIGHT:
        print(f"-> Sage's result: {I.to_sage(K).ideal_class_log()}")
        print()
    #clgp.check_dlog((D,), I2, e_g)
    #clgp.check_dlog_exact((D,), I, g, -vector(e), d_parent = d_parent, food=food)
    #assert I.to_sage(K) * fb.prod_ideal(K, FB, e) == K.ideal(g.to_sage(K))
    print()
    return (g,e)
