# Computation of square roots for S-smooth ideal I by reducing the problem to cyclic subgroups of the class group.
# Requires precomputed class group.

import clgp

# Solving the discrete logarithm problem for an element in cyclic group of 2-power order.
# Input: 
# 1) d, FB are the field and factor base
# 2) an element h = g_i^(ei*t)
# 3) the group <g_i^t>, where g_i is a generator of cyclic subgroup of the class group. 
# 4) r is s.t. #<g_i^t> = 2^r
# Output: l s.t. (g_i^t)^l = g_i^(ei*t)
# Note: in order to avoid evaluation of prime ideals products we work with exponents only.
def dlog_pow_2(ei, t, r, bi):
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
        l = l + 2^k * d_k
    return l

# Extraction of square root of element a in cyclic group <g_i> of order b_i.
# Input: an exponent e_i s.t. a = g_i^e_i.
# Output: list of square roots (exponents of group elements, not reduced) or None if such square roots doesn't exits.  
def cyc_sqrt(e_i, b_i):
    if b_i == 0:
        return [0]
    # Generator of order 2 group is always non square.
    if b_i == 2 and Mod(e_i, 2) == 1:
        print(f"[cyc_sqrt] sqrt doesn't exists for e_i = {e_i}, b_i = {b_i}")
        return None
    if b_i == 2 and Mod(e_i, 2) == 0:
        return [0, 1]
    if Mod(b_i, 2) == 1:
        #return [ZZ(lift(Mod((b_i+1)/2 * e_i, b_i)))]
        return [ZZ((b_i+1)/2 * e_i)] # we expect that reduction will be done by caller
    r = valuation(b_i, 2)
    t = b_i / 2^r
    l = dlog_pow_2(e_i, t, r, b_i)
    # In the group of order 2^r we can have l = 1 | 2^r. If l != 1 and Mod(l,2) != 0 then look for an error in DLOG computation.
    if Mod(l, 2) != 0:
        print(f"[cyc_sqrt] sqrt doesn't exists for e_i = {e_i}, b_i = {b_i}")
        return None
    #assert Mod(l, 2) == 0, f"Result of DLOG_pow2 should be even (l = {l}), otherwise sqrt(J) does not exist!"
    sr_1 = e_i * (t+1) / 2 - t*l / 2 
    #sr_1 = ZZ(lift(Mod(sr_1, b_i))) # we expect that reduction will be done by caller

    sr_2 = sr_1 + b_i / 2
    #sr_2 = ZZ(lift(Mod(sr_2, b_i))) # we expect that reduction will be done by caller
    return [ZZ(sr_1), ZZ(sr_2)]

# Computation of square root (in class group) of smooth ideal class J = prod P[i]^e[i] for P[i] in the factor base.
# Returns list of pairs (h, e_sqrt) such that prod_i P[i]^h[i] is a principal ideal and
# J = (prod_i P[i]^h[i]) * (prod_i P[i]^e_sqrt[i])^2.
#
# Requires precomputed relations (run testrelations.sage with the same 'food').
def smooth_ideal_sqrt(d, e, d_parent=(), food=None):
    #print(f"-> e = {e}")

    M,B,U,V,U_inv,V_inv = clgp.get_matrices(d)

    clgp.check_clgp(d, d_parent, food=food)

    # Write the input ideal J in terms of class group generators
    e_g = clgp.primes_to_gens(e, V, B, reduce=True)
    #print(f"-> Input ideal exponents in terms of Cl_K gens, e_g = {e_g}")

    # check that ideal is written correctly in terms of generators of class group.
    #assert clgp.check_prods_p2g(d, e, e_g, d_parent, food=food), "Wrong conversion prod P_i^e_i =>  prod g_i^et_i"

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
        h, J_sqrt = clgp.reduce_exp(d, J_sqrt) # returns h = prod_i P_i and J_sqrt = prod_j g_j s.t. [(h J_sqrt)^2] = [J]

        #print(f"h = {h}\nJ_sqrt = {J_sqrt}\ne = {e}\ne_g = {e_g}")
        e_sqrt = clgp.gens_to_primes(J_sqrt, V_inv)
        #assert clgp.check_prods_p2g(d, e_sqrt, J_sqrt, d_parent, food=food), "Wrong conversion from prod g[i]^e_g[i] to P[i]^e[i]"
        #print(f"-> e_sqrt = {e_sqrt}")
        #clgp.check_smooth_sqrt(d, e, e_sqrt, d_parent, food=food)

        # principal multiplier
        h = vector(e) - 2*vector(e_sqrt)

        e_sqrts.append((h, e_sqrt))
    return e_sqrts
