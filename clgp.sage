import field
import relations
import fb
import verify
from memoized import memoized

from fpylll import GSO, IntegerMatrix, LLL

@memoized
def get_FB(d, d_parent = (), food = None):
    if tuple(d_parent) == ():
        file = relations.convert_letters(tuple(d), seed=1, food=food)
    else:
        file = relations.convert_letters(tuple(d_parent), seed=1, food=food)
    FB = fb.load_primes(tuple(d), file, food=food)
    return FB

# Converting prod_i P_i^e_i => prod_i g_i^e'_i.
# Here P_1, P_2, ... are primes from the factor base and g_1, g_2, ... are class group generators.
# V is a matrix from SNF of relations matrix A, i.e. A = U B V.
def primes_to_gens(e, V, B = None, reduce = False, strip_zeroes=False):
    e_g = [sum(ZZ(V[i,j]*e[i]) for i in range(V.ncols())) for j in range(V.nrows())]
    assert B.ncols() == len(e_g)
    if reduce and B != None:
        e_g = [ZZ(lift(Mod(e_g[i], B[i,i]))) for i in range(len(e_g))]
    if strip_zeroes:
        for i in range(len(e_g)):
            if B[i,i] != 1 and B[i,i] != 0:
                e_g = e_g[i:]
                break
        if (len(e_g) == 0) or (i == len(e_g)-1 and (e_g[i] == 0 or e_g[i] == 1)):
            e_g = []
    return e_g

# Converting prod_i prod_i g_i^e'_i => P_i^e_i.
# Here P_1, P_2, ... are primes from the factor base and g_1, g_2, ... are class group generators.
# V_inv = V^(-1) is a inverse of matrix V from SNF of relations matrix A, i.e. A = U B V.
def gens_to_primes(e_g, V_inv):
    e = [sum(ZZ(V_inv[i,j]*e_g[i]) for i in range(V_inv.ncols())) for j in range(V_inv.nrows())]
    return e

# Matrix of relations with SNF decomposition
@memoized
def get_matrices(d):
    M = relations.load_matrix(d)
    #print("\n\nM:", M)
    B,U,V=M.smith_form()
    U_inv = U^-1
    V_inv = V^-1
    return (M,B,U,V,U_inv,V_inv)

def inverse(d, e_g):
    M,B,U,V,U_inv,V_inv = get_matrices(d)
    #print(f"[clgp_inv]: B = {[B[i,i] for i in range(B.nrows()) if B[i,i] not in [0,1]]}")
    if len(e_g) == B.ncols():
        e_g = [ZZ(lift(Mod(-e_g[i], B[i,i]))) for i in range(len(e_g))]
    else:
        j = 0
        for i in range(B.ncols()):
            if B[i,i] != 1:
                j = i
                break 
        e_g = [ZZ(lift(Mod(-e_g[i], B[i+j,i+j]))) for i in range(len(e_g))]
    return e_g

# Computes principal ideals g_i^b_i, where g_i are class group generators and b_i = #<g_i>.
# Results are represented as ideals products of the form prod_i P_i^e_i.
@memoized
def get_gens_ids(d):
    M,B,U,V,U_inv,V_inv = get_matrices(d)
    g = [0]*B.nrows()
    for i in range(len(g)):
        if B[i,i] == 0 or B[i,i] == 1:
            g[0] = [0] * M.nrows()
            continue
        g[i] = [sum(M[j,k]*U[i,j] for j in range(M.nrows())) for k in range(M.ncols())]
    return g

# Reduces exponents (modulo orders of class group generators) of an element g of class group
# represented by vector e s.t. g = prod_i(g_i^e_i), where {g_i}_i are generators of class group. 
# Returns pair (h, e') s.t. prod_i(g_i^e_i) = h prod_i(g_i^e'_i), where h is vector that represents
# a principal ideal of the form prod_i(P_i^h_i). Here P_i are prime ideals in factor base.
def reduce_exp(d, e_g):
    M,B,U,V,U_inv,V_inv = get_matrices(d)
    gens = get_gens_ids(d)
    assert len(e_g) == B.nrows()
    res = [0] * len(e_g)
    h = vector(ZZ, M.ncols())
    for i in range(len(e_g)):
        if B[i,i] == 0 or B[i,i] == 1:
            res[i] = 0
            continue
        a,b = divmod(e_g[i], B[i,i])
        # e_g[i] = a*b_i + b
        res[i] = b
        h += a * vector(gens[i])
    return (list(h), res)

# Reduce exponent vector v using Nearest Plane algorithm
def reduce_NP(d, v):
    A = relations.load_matrix(d)
    A = matrix([r for r in A.rows() if not r.is_zero()])
    A0 = IntegerMatrix.from_matrix(A)
    #A0 = LLL.reduction(A0)
    M = GSO.Mat(A0, float_type="mpfr")
    _ = M.update_gso()
    w = M.babai(v)
    w = A0.multiply_left(w)
    return [ZZ(v[i]) - ZZ(w[i]) for i in range(len(v))]

# Computes generators of class group as explicit products of prime ideals.
# Extreemly slow, for testing purposes only. 
@memoized
def generators_explicit(d, d_parent = (), food = None):
    FB = get_FB(d, d_parent, food = food)
    K = field.field(d).sage(food)
    M,B,U,V,U_inv,V_inv = get_matrices(d)
    g = []
    for i in range(B.ncols()):
        gi_v = [V_inv[i,j] for j in range(V_inv.ncols())]
        gi = fb.prod_ideal(K, FB, gi_v)
        g.append(gi)
    return g

def strip_ones(a):
    for i in range(len(a)):
        if a[i] != 1:
            return a[i:]
    return []

def strip_oz(a):
    r = []
    for i in range(len(a)):
        if a[i] != 1 and a[i] !=0:
            r.append(a[i])
    return r

def ltrim(a, l=None, s=[0], j=0):
    r = []
    for i in range(len(a)):
        if not (a[i] in s):
            r = a[i:]
            break
    if l != None and j != None:
        if l > len(r):
            r = [j]*(l-len(r)) + r
    return r

# Lift product of prime ideals prod(FB[i]^e[i] for i in 0..len(FB)) from subfield to the parent field
# using splitting data from for prime ideals in the factor base (FB).
def lift_e(e, FB):
    currentp = 0
    e_lift = []
    assert len(FB)==len(e), f'{len(FB)}!={len(e)}'
    for i in range(len(FB)):
        if currentp != FB[i].prime:
            e_lift += [e[i]*j for j in FB[i].powers]
            currentp = FB[i].prime
        else:
            c = len(FB[i].powers)
            e_lift[-c:] = [e_lift[-c:][h] + e[i]*FB[i].powers[h] for h in range(c)]
    return e_lift

@memoized
def prime_ideal_inverse(P):
    return P^(-1)

# Checks equality [prod P[i]^e[i]] = [prod g[j]^e_g[j]] given vectors e and e_g.
# Here {P_i}_i are ideals generating class group, and {g_j}_j are s.t. CL_K = <g_1> x ... <g_r>.
# Uses Sage's methods, for testing purposes only.
def check_prods_p2g(d, e, e_g, d_parent = (), food=None):
    res = True
    if verify.level()==verify.HEAVY:
        K = field.field(d).sage(food)
        FB = get_FB(d, d_parent, food = food)
        Cl_K = K.class_group(proof=False)
        g = generators_explicit(d, d_parent, food = food)
        #print("Verifying input ideal representation in terms of class group generators ...")
        I1 = fb.prod_ideal(K, FB, e)
        I2 = K.ideal(1)
        for i in range(len(e_g)):
            I2 = I2 * (g[i]^e_g[i])
        if Cl_K(I1) != Cl_K(I2):
            res = False
    return res

# Checks that class group computed correctly using Sage's methods.
# For testing and debugging purposes.
def check_clgp(d, d_parent = (), food = None):
    res = True
    # TODO: Find out why we should take the matrix V^-1 instead U^-1.
    # In [Buchmann-Düllmann, p.136, (4)] it is proposed to take U^-1 for transformation, but the code works only with V^(-1).
    # In [Biasse-Vredendaal, p.107, 3A] it is proposed to use V.
    if (verify.level() > verify.LIGHT) or (verify.level() == verify.LIGHT and len(d) <= 5):
        print(f"Verifying class group of the field {d} ... ", end="")
        M,B,U,V,U_inv,V_inv = get_matrices(d)
        K = field.field(d).sage(food)
        #print("-> Computing Cl_K ... ")
        Cl_K = K.class_group(proof=False)
        #print("-> ", Cl_K.elementary_divisors())
        B_t = []
        for i in range(B.ncols()):
            if B[i,i] != 1:
                B_t.append(B[i,i])
        if sorted(B_t) != sorted(Cl_K.elementary_divisors()):
            print(f"[error]\n-> Group structure mismatch! We are working in the group {B_t} instead of the group {sorted(Cl_K.elementary_divisors())}!")
            res = False

        if verify.level() == verify.HEAVY:
            #print("-> Computing explicit class group generators ... ")
            g = generators_explicit(d, d_parent, food=food)
            for i in range(B.ncols()):
                if Cl_K(g[i]).order() != B[i,i]:
                    if res:
                        print("[error]")
                    print(f"-> Wrong generator, order check: {(Cl_K(g[i]).order())} != {B[i,i]}")
                    res = False
                    break
                #print(".", end="")
        if res:
            print("[ok]")
    return res

# Verify that [I] == [prod g[i]^e_g[i]]^(-1) using Sage's ideal_class_log method.
def check_dlog(d, I, e_g):
    res = True
    if verify.level() > verify.LIGHT or (verify.level() == verify.LIGHT and len(d) <= 5):
        e_g_inv = [ZZ(i) for i in inverse(d, e_g)]
        print(f"Verifying dlog {ltrim(e_g)} (eq. {ltrim(e_g_inv)}) ... ", end="")
        dl_sage = I.ideal_class_log()
        if ltrim(sorted(e_g_inv)) != ltrim(sorted(dl_sage)):
            res = False
            print(f"[error]\n-> {ltrim(e_g_inv)} != {dl_sage} (Sage)")
        else:
            print("[ok, ness.cond.]")
    return res

# Given vectors e and e_sqrt checks equality ([I] = [prod P_[i]^e[i]]) == [prod P_[i]^e_sqrt[i]]]^2 in the class group.
def check_smooth_sqrt(d, e, e_sqrt, d_parent = (), food=None):
    if verify.level() > verify.LIGHT or (verify.level() == verify.LIGHT and len(d) <= 5):
        K = field.field(d).sage(food)
        FB = get_FB(d, d_parent, food = food)
        Cl_K = K.class_group(proof=False)
        M,B,U,V,U_inv,V_inv = get_matrices(d)
        e_g = primes_to_gens(e, V, B, strip_zeroes=True)
        e_sqrt_g = primes_to_gens(e_sqrt, V, B, strip_zeroes=True)

        print(f"Verifying that {e_sqrt_g} is a square root of {e_g} in the class group ... ", end="")
        I_v = fb.prod_ideal(K, FB, e)
        J_sqrt_v = fb.prod_ideal(K, FB, e_sqrt)
        IJ = fb.prod_ideal(K, FB, [e[i] - 2*e_sqrt[i] for i in range(len(e))])
        assert Cl_K(I_v) == Cl_K(J_sqrt_v)^2, "Wrong computation of sqrt(I)"
        assert Cl_K(IJ) == Cl_K.identity(), "Wrong computation of sqrt(I)"
        #assert IJ.is_principal(proof=False), "Wrong computation of sqrt(I)" # slow even for small degree fields
        print("[ok]")

def check_smooth_sqrt_exact(d, h, e_sqrt, g, e, d_parent = (), food=None):
    K = field.field(d).sage(food)
    FB = get_FB(d, d_parent, food=food)

    #print(f"Verifying that I = {(h.to_sage(),e_sqrt)} is a square root of ideal II = {(g.to_sage(),e)} ... ", end="")
    print(f"Verifying square root computation by testing ideal equality I^2 = g J^(-1) == (h sqrt(J))^2 ... ", end="")
    I_square = K.ideal(g.to_sage(K)) * fb.prod_ideal(K, FB, e)
    I = K.ideal(h.to_sage(K)) * fb.prod_ideal(K, FB, e_sqrt)
    assert I_square == I^2

    I2 = K.ideal((h^2).to_sage(K)) * fb.prod_ideal(K, FB, vector(e_sqrt)*2)
    assert I_square == I2
    print("[ok]")

def check_smooth_sqrt_selection_exact(d, I, h, e_sqrt, d_parent = (), food=None):
    K = field.field(d).sage(food)
    FB = get_FB(d, d_parent, food=food)

    print(f"Verifying that during square root computation we selected correct square root ... ", end="")
    I2 = K.ideal(h.to_sage(K)) * fb.prod_ideal(K, FB, e_sqrt)
    assert I.to_sage(K) == I2
    print("[ok]")

def check_dlog_exact(d, I, h, e, d_parent = (), food=None):
    K = field.field(d).sage(food)
    FB = get_FB(d, d_parent, food=food)

    print(f"Verifying dlog by testing ideal equality I = h*J^(-1) ... ", end="")
    I2 = K.ideal(h.to_sage(K)) * fb.prod_ideal(K, FB, e)
    assert I.to_sage(K) == I2
    print("[ok]")