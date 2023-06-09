# Computation of square roots for S-smooth ideal I by reducing the problem to cyclic subgroups of the class group.
# Requires precomputed class group.

import clgp
import idealprod
import field
import field_compact
import saturation
import units

def dlog_pow_2(e, t, r, b):
    '''
    Solving the discrete logarithm problem for an element in cyclic group of 2-power order.
    
    Input: 
    1) d, FB are the field and factor base
    2) an element h = g^(e*t)
    3) the group <g^t>, where g is a generator of cyclic subgroup of the class group. 
    4) r is s.t. #<g^t> = 2^r
    
    Output: l s.t. (g^t)^l = g^(e*t)
    
    Note: in order to avoid evaluation of prime ideals products we work with exponents only.
    '''

    l = 0
    gamma_e = 2^(r-1)*t # gamma = generator^(2^(r-1)) = (g_i^t)^(2^(r-1))
    for k in range(r):
        h_k_e = (-l*t + e*t)*2^(r-1-k) # h_k = (generator^-l * h)^(2^(r-1-k)) = g_i^((-l*t + e*t)*2^(r-1-k))
        # We have to find d_k in {0,1} s.t.
        # gamma^d_k = h_k, i.e. g_i^(d_k*t*2^(r-1)) = g_i^((-l*t + e*t)*2^(r-1-k)) or
        # d_k*t*2^(r-1) = (-l*t + e*t)*2^(r-1-k) mod b_i, where b_i = #<g_i>
        if Mod(h_k_e, b) == 0:
            d_k = 0
        elif Mod(gamma_e, b) == Mod(h_k_e, b):
            d_k = 1
        else:
            raise Exception("DLOG_pow2: Something wrong in assumptions!")
        l = l + 2^k * d_k
    return l

def cyc_sqrts(e, b):
    '''
    Extraction of square root of element a in cyclic group <g> of order b.
    
    Input: an exponent e s.t. a = g^e.
    Output: list of square roots (exponents of group elements, not reduced) or None if such square roots doesn't exits.  
    '''

    if b == 0:
        return [0]
    # Generator of order 2 group is always non square.
    if b == 2 and Mod(e, 2) == 1:
        print(f"[cyc_sqrt] sqrt doesn't exists for e = {e}, b = {b}")
        return None
    if b == 2 and Mod(e, 2) == 0:
        return [0, 1]
    if Mod(b, 2) == 1:
        #return [ZZ(lift(Mod((b+1)/2 * e, b)))]
        return [ZZ((b+1)/2 * e)] # we expect that reduction will be done by caller
    r = valuation(b, 2)
    t = b / 2^r
    l = dlog_pow_2(e, t, r, b)
    # In the group of order 2^r we can have l = 1 | 2^r. If l != 1 and Mod(l,2) != 0 then look for an error in DLOG computation.
    if Mod(l, 2) != 0:
        print(f"[cyc_sqrt] sqrt doesn't exists for e = {e}, b = {b}")
        return None
    #assert Mod(l, 2) == 0, f"Result of DLOG_pow2 should be even (l = {l}), otherwise sqrt(J) does not exist!"
    sr_1 = e * (t+1) / 2 - t*l / 2 
    #sr_1 = ZZ(lift(Mod(sr_1, b))) # we expect that reduction will be done by caller

    sr_2 = sr_1 + b / 2
    #sr_2 = ZZ(lift(Mod(sr_2, b))) # we expect that reduction will be done by caller
    return [ZZ(sr_1), ZZ(sr_2)]

def smooth_ideal_sqrts(d, e, d_parent=(), food=None):
    '''
    Computation of square root (in class group) of smooth ideal class J = prod P[i]^e[i] where P[i] are prime ideals from the factor base.
    
    Returns list of pairs (h, e_sqrt) such that prod_i P[i]^h[i] is a principal ideal and
    J = (prod_i P[i]^h[i]) * (prod_i P[i]^e_sqrt[i])^2.
    
    Requires precomputed relations (run testrelations.sage with the same 'food').
    '''

    M,B,U,V,U_inv,V_inv = clgp.get_matrices(d)

    # Write the input ideal J in terms of class group generators
    e_g = clgp.primes_to_gens(e, V, B, reduce=True)
    #print(f"-> Input ideal exponents in terms of Cl_K gens, e_g = {e_g}")

    # check that ideal is written correctly in terms of generators of class group.
    assert clgp.check_prods_p2g(d, e, e_g, d_parent, food=food), "Wrong conversion prod P_i^e_i =>  prod g_i^et_i"

    J_sqrts = []
    for i in range(len(e_g)):
        egi_sqrt = cyc_sqrts(e_g[i], B[i,i])
        if egi_sqrt == None: # square root doesn't exist.
            return None
        if len(J_sqrts) == 0:
            J_sqrts = [[egi_sqrt[i]] for i in range(len(egi_sqrt))]
        else:
            J_sqrts = [J_sqrts[i] + [egi_sqrt[j]] for i in range(len(J_sqrts)) for j in range(len(egi_sqrt))]

    e_sqrts = []
    for J_sqrt in J_sqrts:
        h, J_sqrt = clgp.reduce_exp(d, J_sqrt) # returns h = prod_i P_i and J_sqrt = prod_j g_j s.t. [(h J_sqrt)^2] = [J]

        e_sqrt = clgp.gens_to_primes(J_sqrt, V_inv)
        assert clgp.check_prods_p2g(d, e_sqrt, J_sqrt, d_parent, food=food), "Wrong conversion from prod g[i]^e_g[i] to P[i]^e[i]"

        #clgp.check_smooth_sqrt(d, e, e_sqrt, d_parent, food=food)

        # principal multiplier
        h = vector(e) - 2*vector(e_sqrt)

        e_sqrts.append((h, e_sqrt))
    return e_sqrts

def cyc_sqrt(e, b):
    '''
    Extraction of square root of element a in cyclic group <g> of order b.
    
    Input: an exponent e for an element a = g^e.
    Output: a pair (c,d) s.t. (g^(c * x + d))^2 = a for any x in {0,1}.
    '''

    if b == 0:
        return (0, 0)
    # Generator of order 2 group is always non square.
    if b == 2 and Mod(e, 2) == 1:
        print(f"[cyc_sqrt] sqrt doesn't exists for e = {e}, b = {b}")
        return None
    if b == 2 and Mod(e, 2) == 0:
        return (1, 1)
    if Mod(b, 2) == 1:
        return (0, ZZ((b+1)/2 * e))
    r = valuation(b, 2)
    t = b / 2^r
    l = dlog_pow_2(e, t, r, b)
    # In the group of order 2^r we can have l = 1 | 2^r. If l != 1 and Mod(l,2) != 0 then look for an error in DLOG computation.
    if Mod(l, 2) != 0:
        print(f"[cyc_sqrt] sqrt doesn't exists for e = {e}, b = {b}")
        return None
    sr = e * (t+1) / 2 - t*l / 2 

    return (t * b / 2, ZZ(sr))

def smooth_ideal_sqrt(d, e):
    '''
    Computation of square root (in class group CL_K = <g_1> ... <g_t>) of smooth ideal class J = prod P_i^e_i where P_i are prime ideals from the factor base.
    
    Returns the list of pairs of ideals (A_1,B_1), ..., (A_t,B_t) s.t. J / prod_j (A_j^x_j * B_j)^2 is a principal ideal for any x_j in {0,1},
    ord([A_j^x_j]) = 1 or 2, and A_j^x_j *B_j is square root in the cyclic subgroup <g_j> of class group.
    '''

    CL = clgp.clgp(d)
    KC = field_compact.field_compact(d)
    K = field.field(d)

    M,B,U,V,U_inv,V_inv = CL.relations_matrix_SNF()
    IPS = idealprod.idealprods(d, FB=CL.factor_base())

    # Write the input ideal J in terms of class group generators
    e_g = clgp.primes_to_gens(e, V, B, reduce=True)
    #e_g = clgp.primes_to_gens(e, V, B)

    res = []
    for i in range(len(e_g)):
        if B[i,i] == 1 or B[i,i] == 0:
            res.append((IPS.one(), IPS.one()))
            continue

        sqrt_g = cyc_sqrt(e_g[i], B[i,i])
        if sqrt_g == None: # Square root doesn't exist.
            return None
        a,b = sqrt_g
        a_p = clgp.gens_to_primes([0]*i + [a] + [0]*(B.nrows()-i-1), V_inv)
        b_p = clgp.gens_to_primes([0]*i + [b] + [0]*(B.nrows()-i-1), V_inv)
        A_i = IPS(KC.one(), a_p)
        B_i = IPS(KC.one(), b_p)
        res.append((A_i, B_i))
    return res

def ideal_sqrt(d, h, e, SU):
    '''
    Computation of square root of an ideal I = h * prod P_i^e_i where P_i are prime ideals from the factor base.
    
    Returns J = g prod_i P_i^f_i such that I = J^2. Raises exception if such a square root doesn't exist.
    '''

    print("[ideal_sqrt]", end="", flush=True)

    CL = clgp.clgp(d)
    KC = field_compact.field_compact(d)

    #M,B,U,V,U_inv,V_inv = CL.relations_matrix_SNF()
    IPS = idealprod.idealprods(d, FB=CL.factor_base())
    
    print(f"[smooth_ideal_sqrt|", end="", flush=True)
    t = walltime()
    sqrts = smooth_ideal_sqrt(d, e)
    print(f"{walltime(t)} sec.]", end=" ")

    print(f"[gens|", end="", flush=True)
    t = walltime()
    beta_prod = IPS(KC.one(), tuple(e))
    alphas = []
    for i in range(len(sqrts)):
        a,b = sqrts[i]
        beta_prod /= b^2
        gen = (a^2).generator(SU).compute_roots(ideal_sqrt=True)
        alphas.append(gen)
    print(f"{walltime(t)} sec.]", end=" ")

    print(f"[beta|", end="", flush=True)
    t = walltime()
    beta = beta_prod.generator(SU).compute_roots(ideal_sqrt=True)
    print(f"{walltime(t)} sec.]", end=" ")

    print(f"[units|", end="", flush=True)
    t = walltime()
    U = units.generators(d)
    print(f"{walltime(t)} sec.]", end=" ")

    print(f"[find_sq|", end="", flush=True)
    t = walltime()
    x = saturation.find_square_exp(h*beta, alphas + list(U))
    # TODO: we can also use S-units, but it seems like this is unnessary
    print(f"{walltime(t)} sec.]", end=" ")
    
    print(f"[prods|", end="", flush=True)
    t = walltime()
    g = h * beta / prod([alphas[i]^ZZ(x[i]) for i in range(len(alphas))], KC.one())
    G_f = prod([sqrts[i][0]^ZZ(x[i]) * sqrts[i][1] for i in range(len(sqrts))], IPS.one())
    print(f"{walltime(t)} sec.]", end=" ")

    print(f"[sqrt|", end="", flush=True)
    t = walltime()
    rt = (g^(1/2)).compute_roots(ideal_sqrt=True)
    print(f"{walltime(t)} sec.]", end=" ")

    G = IPS(rt, G_f.powers)
    return G

def ideal_sqrt_rat(d, h, e, SU):
    '''
    Computation of square root of an ideal I = h * prod P_i^e_i where P_i are prime ideals from the factor base.
    
    Returns J = g prod_i P_i^f_i such that I = J^2. Raises exception if such a square root doesn't exist.
    
    Assumes that h can have rational exponents in the compact representation.
    '''

    print("[ideal_sqrt_rat]", end="", flush=True)

    CL = clgp.clgp(d)
    KC = field_compact.field_compact(d)

    FB = CL.factor_base()

    IPS = idealprod.idealprods(d, FB=FB)
    
    print(f"[smooth_ideal_sqrt|", end="", flush=True)
    t = walltime()
    sqrts = smooth_ideal_sqrt(d, e)
    print(f"{walltime(t)} sec.]", end=" ")

    print(f"[gens|", end="", flush=True)
    t = walltime()
    beta_prod = IPS(KC.one(), tuple(e))
    alphas = []
    for i in range(len(sqrts)):
        a,b = sqrts[i]
        beta_prod /= b^2
        gen = (a^2).generator(SU)
        alphas.append(gen)
    print(f"{walltime(t)} sec.]", end=" ")

    print(f"[beta|", end="", flush=True)
    t = walltime()
    beta = beta_prod.generator(SU)
    print(f"{walltime(t)} sec.]", end=" ")

    print(f"[units|", end="", flush=True)
    t = walltime()
    U = units.generators(d)
    U = [KC(U[i].element) for i in range(len(U))]
    print(f"{walltime(t)} sec.]", end=" ")

    print(f"[find_sq_rat|", end="", flush=True)
    t = walltime()
    
    T = alphas + list(U)

    x = saturation.find_square_rat_exp(h*beta, T, r = len(FB))
    # TODO: we can also use S-units, but it seems like this is unnessary
    print(f"{walltime(t)} sec.]", end=" ")
    
    print(f"[prods|", end="", flush=True)
    t = walltime()
    #g = h * beta / prod([alphas[i]^ZZ(x[i]) for i in range(len(alphas))], KC.one())
    #g /= prod([U[i]^ZZ(x[i+len(alphas)]) for i in range(len(U))], KC.one())
    g = h * beta / prod([T[i]^ZZ(x[i]) for i in range(len(T))], KC.one())
    G_f = prod([sqrts[i][0]^ZZ(x[i]) * sqrts[i][1] for i in range(len(sqrts))], IPS.one())
    print(f"{walltime(t)} sec.]", end=" ")

    print(f"[sqrt|", end="", flush=True)
    t = walltime()
    rt = g^(1/2)
    print(f"{walltime(t)} sec.]", end=" ")

    G = IPS(rt, G_f.powers)
    return G
