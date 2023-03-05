# Computing discrete logarithm using computation of square root in cyclic subgroups of CL_K decompose I using decomposition of I^2.
# Currently doesn't work incorrectly since we have incomplete relation matrices for small degree subfields.

import field
import ideal
import relations
import fb
import trees
import verify
import clgp
import dlog_quad
import smooth_ideal_sqrt_cyc
import field_compact
from memoized import memoized

from fpylll import GSO, IntegerMatrix, LLL

food = trees.get_food() # load food from trees

d = tuple(sorted(food.values(), key=lambda di: abs(di)))
#d = [di for di in Primes()[2:100] if Mod(di, 4) == 1][0:8]

verify.set(verify.LIGHT)

proof.all(False)

SLOW_SQRT_TEST = False

REDUCE_VECTORS = False # TODO: change principal factor

#pari.allocatemem(1024*1024*1024)
#print(f"Pari stack size: {pari.stacksize()}")

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

# Heuristic check for equality I == g * prod P_i^(e_i)], i = 1, ..., #FB.
def check_sqrt(K, FB, I, e, g):
    nrm = [K.ideal(FB[i].prime, K(FB[i].elts).absolute_norm()^e[i]) for i in range(len(e))]
    if I.q != g.absnorm() * nrm:
        return False

    g = g.to_sage(K)
    #return True
    for i in range(len(e)):
        # TODO: add check for e[i] < 0
        if e[i] > 0:
            P = K.ideal(FB[i].prime, K(FB[i].elts))
            el = I.random_element()
            if not (el/g in P):
            #if not (el in g*P):
                print(f"-> check for e[{i}] = {e[i]} [fail]")    
                return False
            else:
                print(f"-> check for e[{i}] = {e[i]} [ok]") 
        elif e[i] < 0:
            P = K.ideal(FB[i].prime, K(FB[i].elts))
            el = I.random_element()
            P_inv = clgp.prime_ideal_inverse(P)
            if not (el/g in P_inv):
            #if not (el in g*P_inv):
                print(f"-> check for e[{i}] = {e[i]} [fail]")    
                return False
            else:
                print(f"-> check for e[{i}] = {e[i]} [ok]")
    return True

# Heuristic check for equality I == g * prod P_i^(e_i)], i = 1, ..., #FB.
def ideal_eq_test(K, FB, I, e, g):
    assert len(FB) == len(e), f"Wrong exponent vector e, len(e) = {len(e)} != {len(FB)}"
    nrm = prod([K.ideal(FB[i].prime, K(FB[i].elts)).absolute_norm()^e[i] for i in range(len(e))])
    #print(f"[ideal_eq_test] I = {I}\n\te = {e}\n\tN(I) = {I.q.factor()}\n\tN(g J) = ({g.absnorm().factor()}) * ({nrm.factor()})")
    if abs(I.q) != abs(g.absnorm() * nrm): # compare upto units
        #print("\t[fail (norms)]")
        return False
    #else:
    #    print("\t[ok]")
    #    return True

    I_sage = I.to_sage(K)
    for i in range(len(e)):
        P = K.ideal(FB[i].prime, K(FB[i].elts))
        if g.valuation(P) + e[i] != I_sage.valuation(P):
            #print("\t[fail (valuations)]")
            return False
    return True

    #TODO: ideal elements selection test

    g = g.to_sage(K)
    #return True
    for i in range(len(e)):
        # TODO: add check for e[i] < 0
        if e[i] > 0:
            P = K.ideal(FB[i].prime, K(FB[i].elts))
            el = I.random_element()
            if not (el/g in P):
            #if not (el in g*P):
                print(f"-> check for e[{i}] = {e[i]} [fail]")    
                return False
            else:
                print(f"-> check for e[{i}] = {e[i]} [ok]") 
        elif e[i] < 0:
            P = K.ideal(FB[i].prime, K(FB[i].elts))
            el = I.random_element()
            P_inv = clgp.prime_ideal_inverse(P)
            if not (el/g in P_inv):
            #if not (el in g*P_inv):
                print(f"-> check for e[{i}] = {e[i]} [fail]")    
                return False
            else:
                print(f"-> check for e[{i}] = {e[i]} [ok]")
    return True

# Heuristic check for equality I^2 == g * prod P_i^(e_i)], i = 1, ..., #FB.
def ideal_sq_eq_test(K, FB, I, e, g):
    assert len(FB) == len(e), f"Wrong exponent vector e, len(e) = {len(e)} != {len(FB)}"
    nrm = prod([K.ideal(FB[i].prime, K(FB[i].elts)).absolute_norm()^e[i] for i in range(len(e))])
    #print(f"[ideal_sq_eq_test] I = {I}\n\tN(I^2) = {(I.q^2).factor()}\n\tN(g J) = ({g.absnorm().factor()}) * ({(nrm).factor()})")
    if abs(I.q)^2 != abs(g.absnorm() * nrm): # compare upto units
        #print("\t[fail (norms)]")
        return False
    
    I_sage = I.to_sage(K)
    for i in range(len(e)):
        P = K.ideal(FB[i].prime, K(FB[i].elts))
        if g.valuation(P) + e[i] != I_sage.valuation(P) * 2:
            #print("\t[fail (valuations)]")
            return False
    return True

    #TODO: ideal elements selection test

    g = g.to_sage(K)
    #return True
    for i in range(len(e)):
        # TODO: add check for e[i] < 0
        if e[i] > 0:
            P = K.ideal(FB[i].prime, K(FB[i].elts))
            el = I.random_element()
            if not (el/g in P):
            #if not (el in g*P):
                print(f"-> check for e[{i}] = {e[i]} [fail]")    
                return False
            else:
                print(f"-> check for e[{i}] = {e[i]} [ok]") 
        elif e[i] < 0:
            P = K.ideal(FB[i].prime, K(FB[i].elts))
            el = I.random_element()
            P_inv = clgp.prime_ideal_inverse(P)
            if not (el/g in P_inv):
            #if not (el in g*P_inv):
                print(f"-> check for e[{i}] = {e[i]} [fail]")    
                return False
            else:
                print(f"-> check for e[{i}] = {e[i]} [ok]")
    return True

# For an ideal I returns a pair (g,e) s.t. I J = I  * prod_i(P[i]^e[i]) = g O_K,
# where g is a field element in compact representation and P[i] are primes from the
# factor base (precomputed during class group computation). 
def cl_dlog(d, I, d_parent = ()):
    M,B,U,V,U_inv,V_inv = clgp.get_matrices(d)

    K = field.field(d).sage(food)
    KC = field_compact.field_compact(d)

    d = tuple(d)
    if len(d) == 1:
        # should return smooth ideal J and an element g, s.t. I*J = g O_K
        g, e = dlog_quad.dlog_quad(d[0], I, d_parent=d_parent, food=food)
        #return [(g, e)]
        return [(KC(g), [0]*len(e), e)]

    B_t = [ZZ(B[i,i]) for i in range(B.ncols())]
    print(f"Computing DLOG for ideal I = {I} in the group {clgp.strip_oz(B_t)} ... ")

    #print(f"B_t = {clgp.strip_oz(B_t)}")
    if set(B_t) == {1}:
        e = [0]*M.ncols()
        e_g = clgp.primes_to_gens(e, V, B, strip_zeroes=True)
        #print(f"h_K = 1, computation of dlog is trivial")
        g = I.generator()
        print(f"Finished computing DLOG for ideal I = {I} in the group {clgp.strip_oz(B_t)} ... {e_g}")
        print(f"-> I = < {g.to_sage(K)} >\n")
        #assert K.ideal(g.to_sage(K)) == I.to_sage(K), f"Wrong computation for trivial class group (d = {d})"
        return [(KC(g), [0]*len(e), e)]

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

    K_s = field.field(d_s).sage(food)
    K_t = field.field(d_t).sage(food)
    K_st = field.field(d_st).sage(food)

    # computing DLOGs in subfields
    print(f"Lift dlog into subfield {d_s}")
    dlog_s  = cl_dlog(d_s, I_s, d_parent=d)
    #assert I_s.to_sage(K_s) * fb.prod_ideal(K_s, FB_s, dlog_s[0][1]) == K_s.ideal(dlog_s[0][0].evaluate(ideal_sqrt=True).to_sage(K_s))

    print(f"Lift dlog into subfield {d_t}")
    dlog_t  = cl_dlog(d_t, I_t, d_parent=d)
    #assert I_t.to_sage(K_t) * fb.prod_ideal(K_t, FB_t, dlog_t[0][1]) == K_t.ideal(dlog_t[0][0].evaluate(ideal_sqrt=True).to_sage(K_t))

    print(f"Lift dlog into subfield {d_st}")
    dlog_st = cl_dlog(d_st, I_st, d_parent=d)
    #assert I_st.to_sage(K_st) * fb.prod_ideal(K_st, FB_st, dlog_st[0][1]) == K_st.ideal(dlog_st[0][0].evaluate(ideal_sqrt=True).to_sage(K_st))

    #print(f"d_s = {d_s}, dlog_s = {dlog_s}")
    #print(f"d_t = {d_t}, dlog_t = {dlog_t}")
    #print(f"d_st = {d_st}, dlog_st = {dlog_st}")

    print(f"{len(dlog_s)*len(dlog_t)*len(dlog_st)} triples lifted from subfields for the field {d}")

    if len(dlog_s) == 0 or len(dlog_t) == 0 or len(dlog_st) == 0:
        print(f"Failed to find dlog for the ideal I = {I}")
        return []

    K_st = NumberField([x^2 - d_st[i] for i in range(len(d_st))], names=list(sorted(food.keys()))[:len(d_st)])
    #FB_st_sage = [K_st.ideal(FB_st[i].prime, K_st(FB_st[i].elts)) for i in range(len(FB_st))]
    #FB_st_conj = [K_st.ideal(FB_st[i].prime, apply_aut(K_st, K_st(FB_st[i].elts), [0]*(len(d_st)-1) + [1])) for i in range(len(FB_st))]

    FB = clgp.get_FB(d, d_parent, food = food)

    I_sage = I.to_sage(K)

    res = []

    for h_s, h2_s, e_s in dlog_s:
        # Check correstness of DLOG computation in subfield
        assert ideal_eq_test(K_s, FB_s, I_s, -vector(e_s) + vector(h2_s), h_s)
        #assert I_s.to_sage(K_s) * fb.prod_ideal(K_s, FB_s, e_s) == K_s.ideal(h_s.evaluate(ideal_sqrt=True).to_sage(K_s)), f"incorrect dlog computation for the field {K_s}"

        for h_t, h2_t, e_t in dlog_t:
            # Check correstness of DLOG computation in subfield
            assert ideal_eq_test(K_t, FB_t, I_t, -vector(e_t) + vector(h2_t), h_t), f"incorrect dlog computation for the field {K_t}"
            #assert I_t.to_sage(K_t) * fb.prod_ideal(K_t, FB_t, e_t) == K_t.ideal(h_t.evaluate(ideal_sqrt=True).to_sage(K_t))

            for h_st, h2_st, e_st in dlog_st:
                # Check correstness of DLOG computation in subfield
                assert ideal_eq_test(K_st, FB_st, I_st, -vector(e_st) + vector(h2_st), h_st), f"incorrect dlog computation for the field {K_st}"
                #assert I_st.to_sage(K_st) * fb.prod_ideal(K_st, FB_st, e_st) == K_st.ideal(h_st.evaluate(ideal_sqrt=True).to_sage(K_st))

                #print("-> e_s:", e_s)
                #print("-> e_t:", e_t)
                #print("-> e_st:", e_st, "\n")

                e_s_lift = clgp.lift_e(e_s, FB_s)
                e_t_lift = clgp.lift_e(e_t, FB_t)

                g2_s = clgp.lift_e(h2_s, FB_s)
                g2_t = clgp.lift_e(h2_t, FB_t)

                # Since factor base contains sigma-conjugated primes, we need to convert exponents
                # TODO: find more effective way
                # e_st_2 = []
                # for i in range(len(e_st)):
                #     j = FB_st_sage.index(FB_st_conj[i])
                #     ei_conj = e_st[j]
                #     #print(f"i = {i} => j = {j}| {FB_st[i].prime}, {FB_st[i].elts} => {FB_st[j].prime}, {FB_st[j].elts}")
                #     #assert i != j or K_st(FB_st[i].elts).is_rational()
                #     e_st_2.append(ei_conj)
                # e_st = e_st_2
                #print(f"-> e_st.conj = {e_st}")

                # We assume that the automorphism sigma is applied to the prime ideals in trees for the subfield field K_st.
                e_st_lift = clgp.lift_e(e_st, FB_st)
                g2_st = clgp.lift_e(h2_st, FB_st)

                #print(f"-> e_s_lift {d_s} => {d}: {e_s_lift}")
                #print(f"-> e_t_lift {d_t} => {d}: {e_t_lift}")
                #print(f"-> e_st_lift {d_st} => {d}: {e_st_lift}")

                assert(len(e_s_lift) == len(e_t_lift))
                assert(len(e_t_lift) == len(e_st_lift))

                # applying norm equation
                e = [e_s_lift[i] + e_t_lift[i] - e_st_lift[i] for i in range(len(e_s_lift))]
                #print(f"-> e = {e} \n")

                # applying norm equation for second principal part
                g2 = [g2_s[i] + g2_t[i] - g2_st[i] for i in range(len(g2_s))]
                #print(f"g2 = {g2}")

                #print(f"h_s = {h_s.to_sage(K)}")
                #print(f"h_t = {h_t.to_sage(K)}")
                #print(f"h_st = {h_st.to_sage(K)}")
                #g_s = field.field(d)(h_s)
                #g_t = field.field(d)(h_t)
                #g_st = field.field(d)(h_st)
                g_s = KC(h_s)
                g_t = KC(h_t)
                g_st = KC(h_st)

                #print(f"g_s = {g_s}")
                #print(f"g_t = {g_t}")
                #print(f"g_st = {g_st}")

                #print(f"g_s = {g_s.to_sage(K)}, nrm. = {g_s.absnorm().factor()}")
                #print(f"g_t = {g_t.to_sage(K)}, nrm. = {g_t.absnorm().factor()}")
                #print(f"g_st.conj = {g_st.conj().to_sage(K)}, , nrm. = {g_st.conj().absnorm().factor()}")
                #print(f"g_st = {g_st.to_sage(K)}, , nrm. = {g_st.absnorm().factor()}")

                #print(f"g_st = {g_st.to_sage()}")
                #print(f"g (sage) = {g_s.to_sage(K)*g_t.to_sage(K) / g_st.conj().to_sage(K)}")
                #g = g_s * g_t / g_st
                g = g_s * g_t / g_st.conj()

                assert ideal_sq_eq_test(K, FB, I, -vector(e) + vector(g2), g) # N(I prod P_i^e_i) == N(g) prod P_i^(g2_i)

                #try:
                #    g_sqrt = g.sqrt(ideal_sqrt=True)
                #except:
                #    g_sqrt = g^(1/2)
                
                #g_sqrt = ideal.idealsqrtshorten(len(d),2^len(d),d,g)
                #print(f"sqrt(g O_K) = ({g_sqrt.to_sage()})")
                #print(f"N(g) = {g.absnorm().factor()}, N(sqrt(g)) = {g_sqrt.absnorm().factor()}")

                #g2_sqrt = vector(g2) / 2
                #print(f"sqrt(g2) = {g2_sqrt}")

                #CL = K.class_group(proof=False)
                #assert {CL(I.to_sage(K)^2) == CL(K.ideal(g.evaluate(ideal_sqrt=True).to_sage(K)) * fb.prod_ideal(K, FB, e))}
                #assert {CL(I.to_sage(K)^2) == CL(fb.prod_ideal(K, FB, e))}
                assert ideal_sq_eq_test(K, FB, I, -vector(e) + vector(g2), g)
                #assert I.to_sage(K)^2 * fb.prod_ideal(K, FB, e) == K.ideal(g.evaluate(ideal_sqrt=True).to_sage(K)), "Wrong generator"

                #print(f"I.is_principal() = {I.to_sage(K).is_principal()}")
                #print(f"(I^2).is_principal() = {(I.to_sage(K)^2).is_principal()}")

                e_g = clgp.primes_to_gens(e, V, B, strip_zeroes=True)
                print(f"Computing smooth square root for {e_g} in class group {clgp.strip_oz(B_t)} of the field {d}... ")
                e_sqrts = smooth_ideal_sqrt_cyc.smooth_ideal_sqrt(d, -vector(e), d_parent=d_parent, food=food)
                if e_sqrts == None:
                    #print("--> sqrt doesn't exists, skipping this branch")
                    print("[fail]")
                    continue
                print(f"Finished computing smooth square root for {e_g} in class group {clgp.strip_oz(B_t)}: {len(e_sqrts)} found")

                print(f"Selecting suitable root ...")
                # smooth_ideal_sqrt_cyc should return pairs (h, e_sqrt) s.t. J^-1 = h prod_i P_i^(2*e_sqrt[i])
                # FB_sage = [K.ideal(FB[i].prime, K(FB[i].elts)) for i in range(len(FB))]
                for h,e_sqrt in e_sqrts:
                    #print(f"h = {h}")
                    #print(f"e_sqrt = {e_sqrt}")
                    #print(f"-e = {-vector(e)}")
                    assert -vector(e) == vector(e_sqrt) * 2 + vector(h)
                    #nrm1 = prod([FB_sage[i].absolute_norm()^(2*e_sqrt[i]+h[i]) for i in range(len(e_sqrt))])
                    #nrm2 = prod([FB_sage[i].absolute_norm()^(-vector(e)[i]) for i in range(len(e))])
                    #assert nrm1 == nrm2, f"Field {d}. Wrong sqrt of smooth ideal: {e_sqrt}, N(h sqrt(J)^2) = {nrm1.factor()} != N(J) = {nrm2.factor()}"
                    #assert fb.prod_ideal(K, FB, h).is_principal()
                    #print(f"sqrt(g2) + sqrt(h) + sqrt(g): {vector(g2_sqrt) + vector(h)/2 + vector([g.valuation(FB_sage[i]) for i in range(len(FB_sage))])/2}")
                    #assert fb.prod_ideal(K, FB, vector(g2_sqrt) + vector(h)/2 + vector([g.valuation(FB_sage[i]) for i in range(len(FB_sage))])/2).is_principal()
                    #print(f"result: {(g_sqrt, vector(g2_sqrt) + vector(h)/2, -vector(e_sqrt))}")

                    # FIXME: The following generators should be computed using S-units calculations.
                    # Computing generator for prod_i P_i^h[i].
                    h_fb = fb.prod_ideal(K, FB, h)
                    #assert h_fb.is_principal()
                    hO_gen = h_fb.gens_reduced()
                    assert len(hO_gen) == 1
                    hO_gen = hO_gen[0]

                    # Computing generator for prod_i P_i^h[i].
                    g2_fb = fb.prod_ideal(K, FB, g2)
                    #assert g2_fb.is_principal()
                    g2O_gen = g2_fb.gens_reduced()
                    assert len(g2O_gen) == 1
                    g2O_gen = g2O_gen[0]

                    try:
                        hh = KC(field.field(d).from_sage(hO_gen))
                        g2gg = KC(field.field(d).from_sage(g2O_gen))

                        rt = (hh*g2gg*g).sqrt(ideal_sqrt=True)
                        #print(f"rt = {rt}")
                        e_g = clgp.primes_to_gens(-vector(e_sqrt), V, B, strip_zeroes=True)
                        print(f"Finished computing DLOG for ideal I = {I} in the group {clgp.strip_oz(B_t)} ... {e_g}")
                        if verify.level() > verify.LIGHT or (verify.level() == verify.LIGHT and len(d) <= 5):
                            print(f"-> Sage's result: {I_sage.ideal_class_log()}")
                            print()
                        #return [(g_sqrt, vector(g2_sqrt) + vector(h)/2, -vector(e_sqrt))]
                        return [(rt, [0]*len(FB), -vector(e_sqrt))]
                    except:
                        pass
                print("-> [fail]")
                #return []
                continue    
    return res

trees_food = trees.get_food()
if food != None:
    assert food == trees_food, f"Run trees generation and relation computation first! Trees are generated for {trees_food} != {food}."

bound = prod(di.abs() for di in d)
I = ideal.random(d, bound=bound)

# The ideal can be fixed as in the following examples.
#I = ideal.ideals([5, 229, 257])(669289,(303651, 6337, 202003))
#I = ideal.ideals([-7, -11, -19, -23])(2557,(143, 1040, 485, 1051))
#I = ideal.ideals((5, 13, 17, 29, 37, 41))(9032339,(2001020, 3845142, 269749, 1568238, 3382698, 560462))
#I = ideal.ideals((5, 13, 17, 29, 37))(744251,(275720, 225315, 139312, 8627, 29735))
#I = ideal.ideals((5, 13, 17, 29, 37))(271879,(49480, 115130, 12857, 3796, 90974))
#I = ideal.ideals((5, 13, 17, 29, 37))(879941,(341095, 173263, 46592, 397104, 411054))
#I = ideal.ideals((-7, -11, -19))(709,(213, 46, 119))
#I = ideal.ideals((-3, -7, -11, -19))(1213,(435, 197, 231, 206))
#I = ideal.ideals((5, 13, 17, 29))(5279,(1868, 218, 262, 1955))

print(f"Target ideal: {I}")
K = NumberField([x^2 - d[i] for i in range(len(d))], names=list(sorted(food.keys()))[:len(d)])
I_hnf = I.hnf(K)
print(f"Target ideal (hnf): {I_hnf}")
I_sage = I.to_sage(K)
print(f"I_sage.norm = {I_sage.absolute_norm()}")
assert I_sage.absolute_norm() == I.q

FB = clgp.get_FB(d, d_parent = (), food=food)

assert not fb.is_smooth(I.q, FB), "This code is for non-smooth ideals. For smooth ideal run dlog_smooth.sage"

dls = cl_dlog(d, I)

M,B,U,V,U_inv,V_inv = clgp.get_matrices(d)

print("Final results:")
for i in range(len(dls)):
    g,gg,dl = dls[i]
    print(f"-> i = {i}, {clgp.primes_to_gens(dl, V, B, strip_zeroes=True)}")
    # TODO: check that gg is principal (solve matrix equation with relation matrix)
    print(f"-> ideal eq. test: {ideal_eq_test(K, FB, I, -vector(dl) + vector(gg), g)}") # I = g * gg * prod(P_i^-dl[i])

if verify.level() > verify.LIGHT or (verify.level() == verify.LIGHT and len(d) <= 5):
    # check correctness using Sage (slow):
    print("\nChecks using Sage:")
    # This is not really helpful since Sage uses different generators and
    # the result may be different or have different order of elements.
    print(f"Sage's result = {I_sage.ideal_class_log()}")

    CL_K = K.class_group(proof=False)
    for i in range(len(dls)):
        g,gg,dl = dls[i]
        e_g = clgp.primes_to_gens(dl, V, B, strip_zeroes=True)
        print(f"i = {i}, checking {e_g} by direct computation of ideal product ...")
        J = fb.prod_ideal(K, FB, dl)
        if CL_K(I_sage) != CL_K(J^(-1)):
            print("-> [wrong result]")
        else:
            print("-> [ok]")
