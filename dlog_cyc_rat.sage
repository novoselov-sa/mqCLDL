# Computing discrete logarithm using computation of square root in cyclic subgroups of CL_K decompose I using decomposition of I^2.

import field
import ideal
import relations
import fb
import trees
import verify
import clgp
import dlog_quad
import ideal_sqrt_cyc
import field_compact
import idealprod
import units
import sunits
from memoized import memoized

from fpylll import GSO, IntegerMatrix, LLL

import cProfile
import pstats
from pstats import SortKey

profile = cProfile.Profile()

food = trees.get_food() # load food from trees


verify.set(verify.NONE)

proof.all(False)

# Folder to save results of DLOG computations for the field and its subfields. Set it to None if you don't want to store results
SAVE_ALL_RESULTS = "experiments/dlogs"

pari.allocatemem(600*1024^3)
print(f"Pari stack size: {pari.stacksize()}")

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
def ideal_eq_test(K, FB, I, e, g):
    assert len(FB) == len(e), f"Wrong exponent vector e, len(e) = {len(e)} != {len(FB)}"
    known_factors = [I.q] + [FB[i].prime for i in range(len(FB))]
    nrm = prod([K.ideal(FB[i].prime, K(FB[i].elts)).absolute_norm()^ZZ(e[i]) for i in range(len(e))])
    if abs(I.q) != abs(g.absnorm(known_factors=known_factors) * nrm): # compare upto units
        print("\t[fail (norms)]")
        return False
    #else:
    #    print("\t[ok]")
    #    return True

    I_sage = I.to_sage(K)
    for i in range(len(e)):
        P = K.ideal(FB[i].prime, K(FB[i].elts))
        if g.valuation(P) + e[i] != I_sage.valuation(P):
            print("\t[fail (valuations)]")
            return False
    
    #TODO: ideal elements selection test
    return True


# Heuristic check for equality I^2 == g * prod P_i^(e_i), i = 1, ..., #FB.
def ideal_sq_eq_test(K, FB, I, e, g):
    assert len(FB) == len(e), f"Wrong exponent vector e, len(e) = {len(e)} != {len(FB)}"
    nrm = prod([K.ideal(FB[i].prime, K(FB[i].elts)).absolute_norm()^ZZ(e[i]) for i in range(len(e))])
    known_factors = [I.q] + [FB[i].prime for i in range(len(FB))]
    #print(f"[ideal_sq_eq_test] I = {I}\n\tN(I^2) = {(I.q^2).factor()}\n\tN(g J) = ({g.absnorm().factor()}) * ({(nrm).factor()})")
    if abs(I.q)^2 != abs(g.absnorm(known_factors = known_factors) * nrm): # compare upto units
        #print("\t[fail (norms)]")
        return False
    
    I_sage = I.to_sage(K)
    for i in range(len(e)):
        P = K.ideal(FB[i].prime, K(FB[i].elts))
        if g.valuation(P) + e[i] != I_sage.valuation(P) * 2:
            #print("\t[fail (valuations)]")
            return False
    
    #TODO: ideal elements selection test
    return True

def save_result(d, G, I):
    if SAVE_ALL_RESULTS != None:
        fld = "_".join([str(d[i]) for i in range(len(d))])
        fn = f"{SAVE_ALL_RESULTS}/{fld}_dlog"
        with open(fn, 'w') as f: f.write(str(G))
        with open(f"{fn}_input", 'w') as f: f.write(str(I))\

def clean_results():
    trees.clear_dir("experiments/dlogs")

def cl_dlog(d, I, d_parent = ()):
    '''
    For an ideal I returns an ideal product g * prod_i(P[i]^e[i]) = I,
    where g is a field element in compact representation and P[i] are primes from the
    factor base (precomputed during class group computation)
    '''

    print(f"Init data for d = {d} ... ", end="", flush=True)
    t = walltime()
    print("[fields]", end=" ", flush=True)
    K = field.field(d).sage(food)
    KC = field_compact.field_compact(d)
    print("[clgp]", end=" ", flush=True)
    CL = clgp.clgp(d, d_parent=d_parent)

    print("[matrices]", end=" ", flush=True)
    M,B,U,V,U_inv,V_inv = clgp.get_matrices(d)

    #FB = clgp.get_FB(d, d_parent)
    FB = CL.factor_base()
    print("[fb]", end=" ", flush=True)

    IPS = idealprod.idealprods(d, food=food, FB=FB)
    print("[idealprod]", end=" ", flush=True)
    print(f"[ok] {walltime(t)} sec.", flush=True)

    d = tuple(d)
    if len(d) == 1:
        # should return smooth ideal J and an element g, s.t. I*J = g O_K
        g, e = dlog_quad.dlog_quad(d[0], I, d_parent=d_parent, food=food)
        G = IPS(KC(g), -vector(e))
        save_result(d, G, I)
        return [G]

    B_t = [ZZ(B[i,i]) for i in range(B.ncols())]
    print(f"Computing DLOG for ideal I = {I} in the group {clgp.strip_oz(B_t)} ... ", flush=True)
    comp_time = walltime()

    if set(B_t) == {1}:
        e = [0]*M.ncols()
        e_g = clgp.primes_to_gens(e, V, B, strip_zeroes=True)
        #print(f"h_K = 1, computation of dlog is trivial")
        g = I.generator()
        print(f"Finished computing DLOG for ideal I = {I} in the group {clgp.strip_oz(B_t)} ... {e_g}, {walltime(comp_time)} sec.")
        print(f"-> I = < {g.to_sage(K)} >\n")
        #assert K.ideal(g.to_sage(K)) == I.to_sage(K), f"Wrong computation for trivial class group (d = {d})"
        if verify.level() > verify.LIGHT or (verify.level() == verify.LIGHT and len(d) <= 5):
            print(f"-> Sage's result: {I.to_sage(K).ideal_class_log(proof=False)}")
            print()
        G = IPS(KC(g), e)
        save_result(d, G, I)
        return [G]

    d_s  = d[:-1]
    d_t = d[:-2] + d[-1:]
    d_st = d[:-2] + (d[-2]*d[-1],)

    I_s = ideal.ideals(d_s)(I.q, I.s[:-1])
    I_t = ideal.ideals(d_t)(I.q, I.s[:-2] + I.s[-1:])
    I_st = ideal.ideals(d_st)(I.q, I.s[:-2] + (I.s[-2]*I.s[-1],))

    # loading precomputed ramification trees with factor base
    CL_s = clgp.clgp(d_s, d_parent=d)
    CL_t = clgp.clgp(d_t, d_parent=d)
    CL_st = clgp.clgp(d_st, d_parent=d)

    FB_s = CL_s.factor_base()
    FB_t = CL_t.factor_base()
    FB_st = CL_st.factor_base()

    #K_s = field.field(d_s).sage(food)
    #K_t = field.field(d_t).sage(food)
    #K_st = field.field(d_st).sage(food)

    # computing DLOGs in subfields
    print(f"Lift dlog into subfield {d_s}")
    dlog_s  = cl_dlog(d_s, I_s, d_parent=d)

    print(f"Lift dlog into subfield {d_t}")
    dlog_t  = cl_dlog(d_t, I_t, d_parent=d)

    print(f"Lift dlog into subfield {d_st}")
    dlog_st = cl_dlog(d_st, I_st, d_parent=d)

    print(f"{len(dlog_s)*len(dlog_t)*len(dlog_st)} triples lifted from subfields for the field {d}")

    if len(dlog_s) == 0 or len(dlog_t) == 0 or len(dlog_st) == 0:
        print(f"Failed to find dlog for the ideal I = {I}")
        return []

    #K_st = NumberField([x^2 - d_st[i] for i in range(len(d_st))], names=list(sorted(food.keys()))[:len(d_st)])
    #FB_st_sage = [K_st.ideal(FB_st[i].prime, K_st(FB_st[i].elts)) for i in range(len(FB_st))]
    #FB_st_conj = [K_st.ideal(FB_st[i].prime, apply_aut(K_st, K_st(FB_st[i].elts), [0]*(len(d_st)-1) + [1])) for i in range(len(FB_st))]

    I_sage = I.to_sage(K)
    #FB_sage = [FB[i].to_sage(food) for i in range(len(FB))]

    res = []

    for H_s in dlog_s:
        # Check correstness of DLOG computation in subfield
        #assert ideal_eq_test(K_s, FB_s, I_s, H_s.powers, H_s.element), f"incorrect dlog computation for the field {d_s}"

        for H_t in dlog_t:
            # Check correstness of DLOG computation in subfield
            #assert ideal_eq_test(K_t, FB_t, I_t, H_t.powers, H_t.element), f"incorrect dlog computation for the field {d_t}"

            for H_st in dlog_st:
                # Check correstness of DLOG computation in subfield
                #assert ideal_eq_test(K_st, FB_st, I_st, H_st.powers, H_st.element), f"incorrect dlog computation for the field {d_st}"

                print(f"Applying norm equation... ", end="", flush=True)
                t = walltime()

                e_s_lift = clgp.lift_e(H_s.powers, FB_s)
                e_t_lift = clgp.lift_e(H_t.powers, FB_t)

                # We assume that the automorphism sigma is applied to the prime ideals in trees for the subfield field K_st.
                e_st_lift = clgp.lift_e(H_st.powers, FB_st)

                assert(len(e_s_lift) == len(e_t_lift))
                assert(len(e_t_lift) == len(e_st_lift))

                # applying norm equation
                e = [e_s_lift[i] + e_t_lift[i] - e_st_lift[i] for i in range(len(e_s_lift))]
                #e = [e_st_lift[i] - e_s_lift[i] - e_t_lift[i] for i in range(len(e_s_lift))]

                g_s = KC(H_s.element)
                g_t = KC(H_t.element)
                g_st = KC(H_st.element)

                g = g_s * g_t / g_st.conj()

                G = IPS(g, vector(e))

                print(f"{walltime(t)} sec.", flush=True)

                #print(f"I^2 = {G.to_sage(K)}")
                #print(f"-> G.element: {G.element.to_sage(K)}")
                #print(f"-> G.powers: {G.powers}")
                #print(f"-> G.element (norm factor): {[G.element.elements[i].absnorm().factor() ^ G.element.powers[i] for i in range(len(G.element.elements))]}")

                print(f"\n-> loading S-units for {d} ... ", end="", flush=True)
                t = walltime()
                SU = sunits.sunits(d, d_parent=d_parent)
                SU.preload()
                print(f"{walltime(t)} sec.", flush=True)

                clgp.check_clgp(d, d_parent, food=food)

                e_g = clgp.primes_to_gens(e, V, B, strip_zeroes=True, reduce=True)
                print(f"Computing square root for {e_g} in class group {clgp.strip_oz(B_t)} of the field {d}... ", flush=True)

                G = G.sqrt_rat(SU)
            
                #print(f"-> G.element: {G.element.to_sage(K)}")
                #print(f"-> G.powers: {G.powers}")
                
                profile.create_stats()
                pstats.Stats(profile).strip_dirs().sort_stats(SortKey.CUMULATIVE).print_stats(int(50))
                profile.enable()

                e_g = clgp.primes_to_gens(vector(G.powers), V, B, strip_zeroes=True, reduce=True)
                print(f"Finished computing DLOG for ideal I = {I} in the group {clgp.strip_oz(B_t)} ... {e_g}, {walltime(comp_time)} sec.\n", flush=True)
                
                save_result(d, G, I)
                
                if verify.level() > verify.LIGHT or (verify.level() == verify.LIGHT and len(d) <= 5):
                    tt = walltime()
                    print(f"-> Sage's result: {I_sage.ideal_class_log(proof=False)}, {walltime(tt)} sec.")
                    #print(f"-> is_principal(I)?: {I_sage.is_principal(proof=False)}")
                    #print(f"-> is_principal(I^2)?: {(I_sage^2).is_principal(proof=False)}")
                    print()
                
                return [G]
    return res

trees_food = trees.get_food()
if food != None:
    assert food == trees_food, f"Run trees generation and relation computation first! Trees are generated for {trees_food} != {food}."

d = tuple(sorted(food.values(), key=lambda di: abs(di)))
d_parent = ()
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
#I = ideal.ideals((-3, -7, -11, -19, -23))(46027,(14546, 4546, 9354, 18057, 6922))
#I = ideal.ideals((-3, -7, -11, -19, -23))(61051,(7926, 17874, 18728, 23706, 3468))
#I = ideal.ideals((5, 13, 17, 29, 37))(751481,(31543, 188087, 332403, 289351, 293067))

#d = (-3, -7, -11, -19, -23, -31)
#d_parent = (-3, -7, -11, -19, -23, -31, -43)
#I = ideal.ideals((-3, -7, -11, -19, -23, -31))(2034103,(397923, 250885, 121898, 170390, 923694, 761051))
#I = ideal.ideals((-3, -4807))(2034103,(397923, 912336))

print(f"Target ideal: {I}")
#K = NumberField([x^2 - d[i] for i in range(len(d))], names=list(sorted(food.keys()))[:len(d)])
#I_hnf = I.hnf(K)
#print(f"Target ideal (hnf): {I_hnf}")
#I_sage = I.to_sage(K)
#print(f"I_sage.norm = {I_sage.absolute_norm()}")
#assert I_sage.absolute_norm() == I.q

FB = clgp.get_FB(d, d_parent=d_parent)

assert not fb.is_smooth(I.q, FB), "This code is for non-smooth ideals. For smooth ideal run dlog_smooth.sage"

for i in range(len(FB)):
    assert gcd(FB[i].prime, I.q) == 1 # FIXME: change ideal_test_eq

clean_results()

profile.enable()
t = walltime()
dls = cl_dlog(d, I, d_parent = d_parent)
print(f"Computation time: {walltime(t)} sec.", flush=True)
profile.disable()

M,B,U,V,U_inv,V_inv = clgp.get_matrices(d)

print("Final results:")
for i in range(len(dls)):
    G = dls[i]
    print(f"-> i = {i}, {clgp.primes_to_gens(G.powers, V, B, strip_zeroes=True)}")
    K = field.field(d).sage(food)
    print(f"-> ideal eq. test: {ideal_eq_test(K, FB, I, G.powers, G.element)}") # I = g * gg * prod(P_i^-dl[i])

if verify.level() > verify.LIGHT or (verify.level() == verify.LIGHT and len(d) <= 5):
    # check correctness using Sage (slow):
    print("\nChecks using Sage:")
    I_sage = I.to_sage(K)
    # This is not really helpful since Sage uses different generators and
    # the result may be different or have different order of elements.
    t = walltime()
    print(f"Sage's result = {I_sage.ideal_class_log(proof=False)}, {walltime(t)} sec.")

    CL_K = K.class_group(proof=False)
    for i in range(len(dls)):
        G = dls[i]
        e_g = clgp.primes_to_gens(G.powers, V, B, strip_zeroes=True)
        print(f"i = {i}, checking {e_g} by direct computation of ideal product ...")
        J = fb.prod_ideal(K, FB, G.powers)
        if CL_K(I_sage) != CL_K(J^(-1)):
            print("-> [wrong result]")
        else:
            print("-> [ok]")
