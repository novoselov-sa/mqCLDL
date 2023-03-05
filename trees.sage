import os
import nprofile
import field
import prime_decomp

profiling = ['apply_aut','ideal_below','build_tree','build_tree_par','L_to_K_sec', 'primes_above', 'valuation']

NCPUS = 7

USE_RAMIFIED_PRIMES = True # include primes p s.t. p | disc(K)
USE_EXCEPTIONAL_PRIMES = True # include primes p s.t. p | [O_K : Z[theta_K]]

VERIFY = False # for debugging purposes

def field_id(d, parent_id=None, names=None):
    if d == None or len(d) == 0:
        return None
    if names == None:
        names = list(build_alphabet(name="lower"))
    res = "_".join(names[0:len(d)])
    if parent_id != None:
        res += " < " + parent_id
    return res

def compute_bound(d):
    for di in d[2:]:
        if di % 4 != 1: print('one of the generators:', di, 'is not 1 mod 4!')
    r = 0
    if (d[0] % 4 == 1 & d[1] % 4 == 1): r = 0
    if ((d[0] % 4 == 2 | d[0] % 4 == 3) & d[1] % 4 == 1): r = 2
    if ( d[0] % 4 == 2 & d[1] % 4 == 3): r = 3
    disc = (2**r *prod(d))**(2**(len(d) - 1))
    B = floor(12*log(disc).n()^2)
    return B

def determine_bound_B(d): # function from file maketree.sage
    mod4 = True
    r = 3
    for di in d[1:]:
            if di%4 != 1: mod4 = False
    if (d[0]%4 == 1) & mod4: r = 0
    if (d[0]%4 == 2) & mod4: r = 2
    if (d[0]%4 != 2) & ((d[0]%4 == 3) | ~mod4): r = 2
    disc = (2^r*prod(d))^(2^len(d) - 1)
    return floor(12*log(disc).n()^2)

# In the list a of the form [[a_1, e_1], ..., [a_k, e_k]]
# this function finds i s.t. elt = a_i[j] and it sets e_i = e_i + 1.
# If such i doesn't exist the function append element [elt, 1]
# to the list.
def append_or_inc(a, elt, j = None):
    for i in range(len(a)):
        if j != None and a[i][0][j] == elt[j]:
            a[i][1] += 1
            return
        if j == None and a[i][0] == elt:
            a[i][1] += 1
            return
    a.append([elt, 1])

# Computing second generator of ideal I = (p, a) corresponding to the automorphism mu in prime decomposition of p.
def ideal_gen_from_aut(d_K, mu, p, typ):
    assert typ == 'linear' or typ == 'quadratic'
    K = field.field(d_K)
    theta = K.abs_gen()
    theta1 = sum([K.gens()[i] for i in range(len(d_K)) if Mod(d_K[i],p).is_square()], K.zero())
    theta2 = sum([K.gens()[i] for i in range(len(d_K)) if not Mod(d_K[i],p).is_square()], K.zero())
    theta1_a = theta1.apply_aut(mu)
    theta2_a = theta2.apply_aut(mu)
    theta1_a_p = theta1_a.evaluate_mod_ext(p)
    if typ == 'quadratic':
        theta2_a_sq = lift((theta2_a^2).evaluate_mod_ext(p))
        g = theta^2 - theta * 2 * lift(theta1_a_p) + lift(theta1_a_p^2) - theta2_a_sq
    else:
        g = theta + lift(theta1_a_p)
    return g

# Checks that automorphism in ideal is correct.
def check_ideal_aut(d_K, I):
    K = field.field(d_K)
    p, g, meta = I
    mu, typ = meta
    if typ == 'other':
        return mu == []
    g0 = ideal_gen_from_aut(d_K, mu, p, typ)
    return g0 == g

# Computing ideals above prime ideal I using reversed version of [Cohen1993Adv, Alg. 2.5.3].
def ideal_up_hard(d_K, d_L, I, aut=False):
    n = len(d_L)
    assert len(d_K) == n-1
    K = field.field(d_K)
    L = field.field(d_L)
    p, alpha, meta = I
    mu = meta[0]
    typ = meta[1]
    assert typ == 'other' or Mod(L.idx(), p) == 0 # check that we are in hard case.
    F_L = prime_decomp.prime_decomp(d_L, p)
    K_sage = NumberField([x^2 - di for di in d_K], names="a")
    alpha_sage = K_sage(str(alpha.to_sage(names="a")))
    I_sage = K_sage.ideal(p, alpha_sage)
    res = []
    for J,e in F_L:
        p0, beta, meta0 = J
        assert p0 == p
        b = beta.relnorm(d_K)
        b_sage = K_sage(str(b.to_sage(names="a")))
        pr = K_sage.ideal(p, b_sage)
        if I_sage == pr or (Mod(L.discriminant(), p) == 0 and b_sage.valuation(I_sage) > 0):
            if not aut:
                res.append([J, 1])
            else:
                gamma, typ0 = meta0
                if typ0 == 'other':
                    J0 = [p, beta.apply_aut([0]*(len(d_L)-1) + [1]), meta0]   
                else:
                    gamma_new = gamma[:-1] + [lift(Mod(gamma[-1] + 1, 2))]
                    beta_new = ideal_gen_from_aut(d_L, gamma_new, p, typ0)
                    J0 = [p, beta_new, [gamma_new, typ0]]
                res.append([J0, 1])
            if len(res) == 2:
                break
    #print(f"ideal_up_hard: p = {p}, len(res) = {len(res)}")
    if len(res) == 1:
        F_K = prime_decomp.prime_decomp(d_K, p)
        # Since [L:K]=2 and this extension is Galois, we have the case I O_L = P^2 for
        # some prime ideal P only when the number of factors in p O_K and p O_L coincides,
        # the exponents in p O_L are equal to 2, and exponents in p O_L are equal to 1.
        if len(F_K) == len(F_L) and F_L[0][1] == 2 and F_K[0][1] != 2:
            res[0][1] = 2
    assert res != [], f"error in ideal_up_hard, d_K = {d_K}, d_L = {d_L}, p = {p}"
    return res

# Computes ideals above I (from K) in extension L of K.
# Returns a list of (I_1, e_1), ... (I_k, e_k) s.t. I * O_L = I_1^e_1 * ... * I_k^e_k.
# If aut==True, it also applies automorphism sqrt(d_L[-1]) -> -sqrt(d_L[-1]) for usage of ideal in norm equation.
def ideal_up(d_K, d_L, I, aut=False):
    n = len(d_L)
    assert len(d_K) == n-1
    K = field.field(d_K)
    L = field.field(d_L)
    p, alpha, meta = I
    mu = meta[0]
    typ = meta[1]

    if typ == 'other' or Mod(L.idx(), p) == 0:
        if USE_EXCEPTIONAL_PRIMES:
            return ideal_up_hard(d_K, d_L, I, aut)
        else:
            raise NotImplementedError(f"fast ideal_up is not implemented for exceptional primes") 
    
    s = Mod(alpha.const(), p)
    if d_L[-1]*d_L[-2] == d_K[-1]:
        if typ == 'linear' and is_square(Mod(d_L[-2],p)) and is_square(Mod(d_L[-1],p)):
            d1_r = field.fixed_sqrt(d_L[-1], p)
            d2_r = field.fixed_sqrt(d_L[-2], p)
            d1d2_r = field.fixed_sqrt(d_L[-1]*d_L[-2], p)
            assert d1_r * d2_r == d1d2_r, f"sqrt(d_L[-1])*sqrt(d_L[-2]) != sqrt(d_L[-1]*d_L[-2]), check correctness of field.fixed_sqrt()"

            if mu[-1] == 0 and not aut:
                W = [[0,1], [1,0]]
            elif mu[-1] == 1 and not aut:
                W = [[1,1], [0,0]]
            elif mu[-1] == 0 and aut:
                W = [[0,0], [1,1]]
            elif mu[-1] == 1 and aut:
                W = [[1,0], [0,1]]
            
            theta = L.abs_gen()
            res = []
            for w in W:
                gamma = mu[:-1] + w
                theta_a = theta.apply_aut(gamma)
                #res.append([p, theta + lift(theta_a.evaluate_mod_ext(p)), [gamma, "linear"]])
                I = [p, theta + lift(theta_a.evaluate_mod_ext(p)), [gamma, "linear"]]
                append_or_inc(res, I, 1)
            return res
        elif typ == 'linear': # K_st [linear] -> K [quadratic]"
            theta = L.abs_gen()
            t1_L = [L.gens()[i] for i in range(n) if Mod(d_L[i],p).is_square()]
            t2_L = [L.gens()[i] for i in range(n) if not Mod(d_L[i],p).is_square()]

            theta1_L = sum(t1_L, L.zero())
            theta2_L = sum(t2_L, L.zero())

            d1_r = field.fixed_sqrt(d_L[-1], p)
            d2_r = field.fixed_sqrt(d_L[-2], p)
            d1d2_r = field.fixed_sqrt(d_L[-1]*d_L[-2], p)
            assert d1_r * d2_r == d1d2_r, f"sqrt(d_L[-1])*sqrt(d_L[-2]) != sqrt(d_L[-1]*d_L[-2]), check correctness of field.fixed_sqrt()"

            mu2 = [lift(Mod(mu[i] + 1, 2)) for i in range(len(mu))]
            if mu2[-1] == 0 and not aut:
                w = [1,1]
            elif mu2[-1] == 1 and not aut:
                w = [1,0]
            elif mu2[-1] == 0 and aut:
                w = [1,0]
            elif mu2[-1] == 1 and aut:
                w = [1,1]
            
            gamma = mu2[:-1] + w
            res = []
            theta1_a = theta1_L.apply_aut(gamma)
            theta1_a_p = theta1_a.evaluate_mod_ext(p)
            theta2_a = theta2_L.apply_aut(gamma)
            theta2_a_sq = lift((theta2_a^2).evaluate_mod_ext(p))
            g = theta^2 - theta * 2 * lift(theta1_a_p) + lift(theta1_a_p^2) - theta2_a_sq
            #res.append([p, g, [gamma, "quadratic"]])
            I = [p, g, [gamma, "quadratic"]]
            append_or_inc(res, I, 1) #  ideals can coincide, e.g. when p | disc(K)
        elif typ == "quadratic": # K_st [quadratic] -> K [quadratic]"
            theta = L.abs_gen()
            t1_L = [L.gens()[i] for i in range(n) if Mod(d_L[i],p).is_square()]
            t2_L = [L.gens()[i] for i in range(n) if not Mod(d_L[i],p).is_square()]

            theta1_L = sum(t1_L, L.zero())
            theta2_L = sum(t2_L, L.zero())

            d1_r = field.fixed_sqrt(d_L[-1], p)
            d2_r = field.fixed_sqrt(d_L[-2], p)
            d1d2_r = field.fixed_sqrt(d_L[-1]*d_L[-2], p)
            assert d1_r * d2_r == d1d2_r, f"sqrt(d_L[-1])*sqrt(d_L[-2]) != sqrt(d_L[-1]*d_L[-2]), check correctness of field.fixed_sqrt()"

            if mu[-1] == 0 and not aut:
                W = [[1,1], [0,0]]
            elif mu[-1] == 1 and not aut:
                W = [[1,0], [0,1]]
            elif mu[-1] == 0 and aut:
                W = [[1,0], [0,1]]
            elif mu[-1] == 1 and aut:
                W = [[1,1], [0,0]]
            
            res = []
            for w in W:
                theta1_a = theta1_L.apply_aut(mu[:-1] + w)
                theta1_a_p = theta1_a.evaluate_mod_ext(p)
                theta2_a = theta2_L.apply_aut(mu[:-1] + w)
                theta2_a_sq = lift((theta2_a^2).evaluate_mod_ext(p))
                g = theta^2 - theta * 2 * lift(theta1_a_p) + lift(theta1_a_p^2) - theta2_a_sq
                gamma = mu[:-1] + w
                #res.append([p, g, [gamma, "quadratic"]])
                I = [p, g, [gamma, "quadratic"]]
                append_or_inc(res, I, 1) #  ideals can coincide, e.g. when p | disc(K)
        else:
            raise NotImplementedError(f"ideal_up is not implemented for type = {typ}")
    elif d_L[:-1] == d_K:
        if typ == 'quadratic':
            res = []
            theta = L.abs_gen()
            t1_K = [L.gens()[i] for i in range(n-1) if Mod(d_L[i],p).is_square()]
            t2_K = [L.gens()[i] for i in range(n-1) if not Mod(d_L[i],p).is_square()]
            theta1_K = sum(t1_K, L.zero())
            theta2_K = sum(t2_K, L.zero())

            for bit in [0,1]:
                if not is_square(Mod(d_L[-1],p)):
                    theta1_a = theta1_K.apply_aut(mu + [0])
                    theta2_a = theta2_K.apply_aut(mu + [0]) + L.gens()[-1] * (-1)^(bit)
                    theta1_a_p = theta1_a.evaluate_mod_ext(p)
                else:
                    theta1_a = theta1_K.apply_aut(mu + [0]) + L.gens()[-1] * (-1)^(bit)
                    theta2_a = theta2_K.apply_aut(mu + [0])
                    theta1_a_p = theta1_a.evaluate_mod_ext(p)

                theta2_a_sq = lift((theta2_a^2).evaluate_mod_ext(p))
                g = theta^2 - theta * 2 * lift(theta1_a_p) + lift(theta1_a_p^2) - theta2_a_sq
                gamma = mu + [bit]
                #res.append([p, g, [gamma, "quadratic"]])
                I = [p, g, [gamma, "quadratic"]]
                append_or_inc(res, I, 1) #  ideals can coincide, e.g. when p | disc(K)
        elif typ == 'linear' and not is_square(Mod(d_L[-1],p)):
            #raise NotImplementedError(f"ideal_up is not implemented for non-linear factors: case K_s [linear] -> K [quadratic]")
            theta = L.abs_gen()
            t1_L = [L.gens()[i] for i in range(n) if Mod(d_L[i],p).is_square()]
            t2_L = [L.gens()[i] for i in range(n) if not Mod(d_L[i],p).is_square()]

            theta1_L = sum(t1_L, L.zero())
            theta2_L = sum(t2_L, L.zero())

            res = []
            gamma = [lift(Mod(mu[i] + 1, 2)) for i in range(len(mu))] + [0] # + [1] also should works in this case
            theta1_a = theta1_L.apply_aut(gamma)
            theta1_a_p = theta1_a.evaluate_mod_ext(p)
            theta2_a = theta2_L.apply_aut(gamma)
            #print(f"theta1_a = {theta1_a.to_sage()}, theta2_a = {theta2_a.to_sage()}")
            theta2_a_sq = lift((theta2_a^2).evaluate_mod_ext(p))
            g = theta^2 - theta * 2 * lift(theta1_a_p) + lift(theta1_a_p^2) - theta2_a_sq

            #res.append([p, g, [gamma, "quadratic"]])
            I = [p, g, [gamma, "quadratic"]]
            append_or_inc(res, I, 1) #  ideals can coincide, e.g. when p | disc(K)
        elif typ == 'linear':
            #res = [
            #    #[p, L.abs_gen() + L.from_ZZ(lift(s + Mod(d_L[-1], p).sqrt())), [0, "linear"]],
            #    #[p, L.abs_gen() + L.from_ZZ(lift(s - Mod(d_L[-1], p).sqrt())), [1, "linear"]]
            #    [p, L.abs_gen() + lift(Mod(s + field.fixed_sqrt(d_L[-1], p), p)), [mu + [0], "linear"]],
            #    [p, L.abs_gen() + lift(Mod(s - field.fixed_sqrt(d_L[-1], p), p)), [mu + [1], "linear"]]
            #]
            res = []
            for bit in [0,1]:
                I = [p, L.abs_gen() + lift(Mod(s + (-1)^bit * field.fixed_sqrt(d_L[-1], p), p)), [mu + [bit], "linear"]]
                append_or_inc(res, I, 1)
        else:
            raise NotImplementedError(f"ideal_up is not implemented for type = {typ}")
    elif d_L[:-2] + (d_L[-1],) == d_K:
        if typ == 'quadratic':
            res = []
            theta = L.abs_gen()
            t1_K = [L.gens()[i] for i in range(n-2) if Mod(d_L[i],p).is_square()]
            t2_K = [L.gens()[i] for i in range(n-2) if not Mod(d_L[i],p).is_square()]
            if Mod(d_L[-1],p).is_square():
                t1_K.append(L.gens()[-1])
            else:
                t2_K.append(L.gens()[-1])
            theta1_K = sum(t1_K, L.zero())
            theta2_K = sum(t2_K, L.zero())

            for bit in [0,1]:
                if not is_square(Mod(d_L[-2],p)):
                    theta1_a = theta1_K.apply_aut(mu[:-1] + [0] + [mu[-1]])
                    theta2_a = theta2_K.apply_aut(mu[:-1] + [0] + [mu[-1]]) + L.gens()[-2] * (-1)^(bit)
                    theta1_a_p = theta1_a.evaluate_mod_ext(p)
                else:
                    theta1_a = theta1_K.apply_aut(mu[:-1] + [0] + [mu[-1]]) + L.gens()[-2] * (-1)^(bit)
                    theta2_a = theta2_K.apply_aut(mu[:-1] + [0] + [mu[-1]])
                    theta1_a_p = theta1_a.evaluate_mod_ext(p)

                theta2_a_sq = lift((theta2_a^2).evaluate_mod_ext(p))
                g = theta^2 - theta * 2 * lift(theta1_a_p) + lift(theta1_a_p^2) - theta2_a_sq
                gamma = mu[:-1] + [bit] + [mu[-1]]
                #res.append([p, g, [gamma, "quadratic"]])
                I = [p, g, [gamma, "quadratic"]]
                append_or_inc(res, I, 1)
        elif not is_square(Mod(d_L[-2],p)) and typ == 'linear':
            #raise NotImplementedError(f"ideal_up is not implemented for non-linear factors: case K_t [linear] -> K [quadratic]")
            res = []
            theta = L.abs_gen()
            t1_K = [L.gens()[i] for i in range(n-2) if Mod(d_L[i],p).is_square()]
            t2_K = [L.gens()[i] for i in range(n-2) if not Mod(d_L[i],p).is_square()]
            if Mod(d_L[-1],p).is_square():
                t1_K.append(L.gens()[-1])
            else:
                t2_K.append(L.gens()[-1])
            theta1_K = sum(t1_K, L.zero())
            theta2_K = sum(t2_K, L.zero())

            gamma = [lift(Mod(mu[i] + 1, 2)) for i in range(len(mu)-1)] + [0] + [lift(Mod(mu[-1] + 1, 2))]

            bit = 0
            #theta1_a = theta1_K.apply_aut(mu[:-1] + [0] + [mu[-1]])
            theta1_a = theta1_K.apply_aut(gamma)
            #theta2_a = theta2_K.apply_aut(mu[:-1] + [0] + [mu[-1]]) + L.gens()[-2] * (-1)^(bit)
            theta2_a = theta2_K.apply_aut(gamma) + L.gens()[-2] * (-1)^(bit)
            theta1_a_p = theta1_a.evaluate_mod_ext(p)
            theta2_a_sq = lift((theta2_a^2).evaluate_mod_ext(p))
            g = theta^2 - theta * 2 * lift(theta1_a_p) + lift(theta1_a_p^2) - theta2_a_sq
            #res.append([p, g, [gamma, "quadratic"]])
            I = [p, g, [gamma, "quadratic"]]
            append_or_inc(res, I, 1)
        elif typ == 'linear':
            res = []
            for bit in [0,1]:
                I = [p, L.abs_gen() + lift(Mod(s + (-1)^bit * field.fixed_sqrt(d_L[-2], p), p)), [mu[:-1] + [bit] + [mu[-1]], "linear"]]
                append_or_inc(res, I, 1)
        else:
            raise NotImplementedError(f"ideal_up is not implemented for type = {typ}")
    else:
        raise NotImplementedError(f"ideal_up is not implemented for the case {d_K} -> {d_L}")
    return res

# Computes lcm of [O_K:Z[theta_K]] and indices of all subfields of K.
def glob_idx(d_K):
    K = field.field(d_K)
    idx_K = K.idx()
    if len(d_K) == 1:
        return idx_K
    d_s = d_K[:-1]
    d_t = d_K[:-2] + (d_K[-1],)
    d_st = d_K[:-2] + (d_K[-1] * d_K[-2],)
    return lcm(lcm(lcm(idx_K, glob_idx(d_s)), glob_idx(d_t)), glob_idx(d_st))

def find_prime_ideal(K, p, alpha, elts, names):
    alpha_sage = K(str(alpha.to_sage(names=names)))
    I_sage = K.ideal(p, alpha_sage) # this can be slow
    for i in range(len(elts)):
        beta = elts[i]
        beta_sage = K(str(beta.to_sage(names=names)))
        P_sage = K.ideal(p, beta_sage)
        if P_sage == I_sage:
            return i
    return None

@nprofile.profile
def build_tree(d_K, tree, norm_bound=None, primes_bound=None, aut=None, aut2=None, names=None, parent_id=None, primes=None, append_primes=[]):
    K = field.field(d_K)
    K_id = field_id(d_K, parent_id=parent_id, names=names)
    pol = K.absolute_polynomial().to_sage().change_ring(QQ)
    pol_disc = pol.discriminant()
    disc_K = K.discriminant()
    #if g_disc == None: #discriminant of the root
    #    g_disc = disc_K
    idx = pol_disc / disc_K
    #g_idx = lcm(g_idx, idx)

    print("Branch: {}\n-> disc(K) = {}\n-> index = {}".format(K_id, disc_K, idx))

    if norm_bound == None:
        norm_bound = compute_bound(d_K) # GRH
    if primes_bound == None:
        primes_bound = norm_bound

    # List of primes with simple decomposition using [Coh'93, Th. 4.8.13].
    # If USE_EXCEPTIONAL_PRIMES = True then we include also hard primes. 
    # Generated once for the root of the tree.
    if primes == None:
        g_idx = glob_idx(d_K)
        print(f"-> g_idx = {g_idx}")
        primes = []
        p = 2
        while p <= primes_bound:
            p = next_prime(p)
            if Mod(g_idx, p) == 0 and not USE_EXCEPTIONAL_PRIMES:
                print(f"p = {p}, p | [O_K : Z[theta_K]], skipped")
                continue
            if Mod(disc_K, p) == 0 and not USE_RAMIFIED_PRIMES:
                print(f"p = {p}, p | disc(K), skipped")
                continue
            sq = True
            for di in d_K:
                if not is_square(Mod(di, p)):
                    sq = False
                    break
            if sq:
                primes.append(p)
            else:
                #print(f"p = {p}, does not split completely")
                primes.append(p)
    
    for p in append_primes:
        if not p in primes:
            if Mod(g_idx, p) == 0 and not USE_EXCEPTIONAL_PRIMES:
                print(f"[additional primes] p = {p}, p | [O_K : Z[theta_K]], skipped")
            elif Mod(disc_K, p) == 0 and not USE_RAMIFIED_PRIMES:
                print(f"[additional primes] p = {p}, p | disc(K), skipped")
            else:
                primes.append(p)

    tree[K_id] = {}

    if len(d_K) == 1:
        K2 = NumberField([x^2-d_K[0]], names=names[0])
        tree[K_id] = {'primes': {} , 'subfields': []}
        #p = 2
        #while p <= primes_bound:
        for p in primes:
            #if Mod(g_disc, p) == 0 or Mod(g_idx, p) == 0:
            #    p = next_prime(p)
            #    continue
            if not USE_EXCEPTIONAL_PRIMES:
                assert Mod(idx, p) != 0, "There are bad primes in the list: {}!".format(p)
            if not USE_RAMIFIED_PRIMES:
                assert Mod(disc_K, p) != 0, "There are bad primes in the list: {}!".format(p)
            #F = prime_decomp.prime_decomp_fast(d_K, p)
            F = prime_decomp.prime_decomp(d_K, p)
            for i,el in enumerate(F):
                I,e = el
                alpha = I[1]
                mu = I[2][0]
                typ = I[2][1]
                if aut2 != None:
                    alpha = alpha.apply_aut(aut2)
                    if typ != 'other':
                        assert len(aut2) == len(mu)
                        mu = [(aut2[i] + mu[i]) % 2 for i in range(len(mu))]
                # S: apply automorphism s (seems that it is not necessary):
                #if aut2 != None:
                #   #assert Mod(alpha.const(), p) == (-1)^mu * Mod(d_K[0], p).sqrt()
                #    assert Mod(alpha.const(), p) == (-1)^mu[-1] * field.fixed_sqrt(d_K[0], p)
                #   #mu = (mu + 1) % 2
                #    mu[-1] = (mu[-1] + 1) % 2
                #    #alpha = K.abs_gen() + K.from_ZZ(lift((-1)^mu[-1] * Mod(d_K[0], p).sqrt())) #K(tuple(0,lift((-1)^mu * Mod(d_K[0], p).sqrt())),  ZZ(1))
                #    alpha = K.abs_gen() + lift(Mod((-1)^mu[-1] * field.fixed_sqrt(d_K[0], p), p))
                #    assert K2.ideal(p, K2(str(I[1].apply_aut(aut2).to_sage(names=names)))) == K2.ideal(p, K2(str(alpha.to_sage(names=names))))
                
                # FIXME: can't exclude some primes since they may be used in other subfields. Resort to prime bound only
                #alpha_K2 = K2(str(alpha.to_sage(names=names)))
                #J = K2.ideal(p, alpha_K2) # this can be slow
                #if J.absolute_norm() <= norm_bound:
                if True: #J.absolute_norm() <= norm_bound or alpha_K2.is_rational():
                    if not (p in tree[K_id]['primes']):
                        tree[K_id]['primes'][p] = {'elements': [], 'powers': [], 'aut': [], 'type': []}
                    tree[K_id]['primes'][p]['elements'].append(alpha)
                    tree[K_id]['primes'][p]['aut'].append(mu)
                    tree[K_id]['primes'][p]['type'].append(typ)
                    #print("(ideal) I = ({},{}), N(I)={}".format(p, alpha.to_sage(names=names), J.absolute_norm()))
                else:
                    print("(excluded ideal) I = ({},{}), N(I)={}".format(p, alpha.to_sage(names=names), J.absolute_norm()))
            p = next_prime(p)
        return tree

    K_names = names
    #s = [1] * (len(d_K) - 2)  + [ 1, -1]
    #t = [1] * (len(d_K) - 2)  + [-1,  1]
    #st = [1] * (len(d_K) - 2) + [-1, -1]
    s = [0] * (len(d_K) - 2)  + [0, 1]
    t = [0] * (len(d_K) - 2)  + [1, 0]
    st = [0] * (len(d_K) - 2) + [1, 1]
    #
    d_s = d_K[:-1]
    K_s_names = K_names[:-1]
    K_s_id = field_id(d_s,  parent_id=K_id, names=K_s_names)
    #
    d_t = d_K[:-2] + (d_K[-1],)
    K_t_names = K_names[:-2] + [K_names[-1]]
    K_t_id = field_id(d_t,  parent_id=K_id, names=K_t_names)
    #
    d_st = d_K[:-2] + (d_K[-1] * d_K[-2],)
    K_st_names = K_names[:-2] + [K_names[-2] + K_names[-1]]
    K_st_id = field_id(d_st,  parent_id=K_id, names=K_st_names)
    #
    tree[K_id]['subfields'] = [K_s_id, K_t_id, K_st_id]
    #tree[K_id]['subfields'] = [K_s_id, K_t_id]

    build_tree(d_t, tree, primes_bound=primes_bound, norm_bound=norm_bound,
                          parent_id=K_id,
                          names=K_t_names,
                          aut=t,
                          primes=primes)
    build_tree(d_s, tree, primes_bound=primes_bound, norm_bound=norm_bound,
                          parent_id=K_id,
                          names=K_s_names,
                          aut=s,
                          primes=primes)
    aut2 = [0] * (len(d_K) - 2)  + [1] # aut2 = res_{K/K_st}(s)
    build_tree(d_st, tree, primes_bound=primes_bound, norm_bound=norm_bound,
                           parent_id=K_id,
                           names=K_st_names,
                           aut=st, aut2=aut2,
                           primes=primes)

    K_primes = {}
    K_s_primes = tree[K_s_id]['primes']
    print("K_s primes ideal_up ...")
    for p, data in K_s_primes.items():
        if not USE_RAMIFIED_PRIMES:
            assert Mod(field.field(d_s).discriminant(), p) != 0, "There are bad primes in the list: {}!".format(p)
        if not USE_EXCEPTIONAL_PRIMES:
            assert Mod(field.field(d_s).idx(), p) != 0, "There are bad primes in the list: {}!".format(p)
        if K_primes.get(p) == None:
            K_primes[p] = {'elements': [], 'powers': [], 'aut': [], 'type': []}
        for i,alpha_s in enumerate(data['elements']):
            mu = data['aut'][i]
            typ = data['type'][i]
            I = [p, alpha_s, [mu, typ]]
            I_up = ideal_up(d_s, d_K, I)
            elts = K_primes[p]['elements']
            J,e = I_up[0]
            g = len(I_up)
            vec = [0] * (len(data['elements']) * g)
            for J,e in I_up:
                elts.append(J[1])
                K_primes[p]['aut'].append(J[2][0])
                K_primes[p]['type'].append(J[2][1])
                vec[len(elts)-1] = e #1
            K_s_primes[p]['powers'].append(vec)

    K_t_primes = tree[K_t_id]['primes']
    print("K_t primes ideal_up ...")
    for p, data in K_t_primes.items():
        if not USE_RAMIFIED_PRIMES:
            assert Mod(field.field(d_t).discriminant(), p) != 0, "There are bad primes in the list: {}!".format(p)
        if not USE_EXCEPTIONAL_PRIMES:
            assert Mod(field.field(d_t).idx(), p) != 0, "There are bad primes in the list: {}!".format(p)
        if K_primes.get(p) == None:
            K_primes[p] = {'elements': [], 'powers': [], 'aut': [], 'type': []}
        for i,alpha_t in enumerate(data['elements']):
            mu = data['aut'][i]
            typ = data['type'][i]
            I = [p, alpha_t, [mu, typ]]
            I_up = ideal_up(d_t, d_K, I)
            elts = K_primes[p]['elements']
            auts = K_primes[p]['aut']
            typs = K_primes[p]['type']
            #print(f"I_up = {I_up}")
            J,e = I_up[0]
            g = len(I_up)
            vec = [0] * (len(data['elements']) * g)

            for J,e in I_up:
                found = False
                if J[1] in elts:
                    vec[elts.index(J[1])] = e #1
                    found = True
                elif USE_EXCEPTIONAL_PRIMES and ( Mod(field.field(d_t).idx(), p) == 0 or Mod(field.field(d_K).idx(), p) == 0 ):
                    K_sage = NumberField([x^2 - di for di in d_K], names=K_names)
                    i = find_prime_ideal(K_sage, p, J[1], elts, K_names)
                    if i != None:
                        vec[i] = e
                        found = True
                if not found:
                    elts.append(J[1])
                    auts.append(J[2][0])
                    typs.append(J[2][1])
                    vec[len(elts)-1] = e
            K_t_primes[p]['powers'].append(vec)

    print("K_st primes ideal_up ...")
    K_st_primes = tree[K_st_id]['primes']
    for p, data in K_st_primes.items():
        if not USE_RAMIFIED_PRIMES:
            assert Mod(field.field(d_st).discriminant(), p) != 0, f"There are bad primes in the list: {p}!"
        if not USE_EXCEPTIONAL_PRIMES:
            assert Mod(field.field(d_st).idx(), p) != 0, f"There are bad primes in the list: {p}!"
        if K_primes.get(p) == None:
            K_primes[p] = {'elements': [], 'powers': [], 'aut': [], 'type': []}
        for i,alpha_st in enumerate(data['elements']):
            mu = data['aut'][i]
            typ = data['type'][i]
            #print "alpha_st =", alpha_st.to_sage(names=K_st_names)
            I_st = [p, alpha_st, [mu, typ]]
            #print "I_st[1] =", I_st[1].to_sage(names=K_st_names)
            I_up = ideal_up(d_st, d_K, I_st, aut=True)

            if VERIFY:
                #******************************* checks ****************************
                ids2 = ideal_up(d_st, d_K, I_st)
                K2 = NumberField([x^2-di for di in d_K], names=K_names)
                I_st_K2 = K2.ideal(p, K2(str(K(alpha_st).to_sage(names=K_names))))
                I_st_K2_aut = K2.ideal(p, K2(str(K(alpha_st.apply_aut(aut2)).to_sage(names=K_names))))
                if len(ids2) == 2:
                    J1,J2 = ids2
                    J1_K2 = K2.ideal(p, K2(str(J1[0][1].to_sage(names=K_names)))) ^ J1[1]
                    J2_K2 = K2.ideal(p, K2(str(J2[0][1].to_sage(names=K_names)))) ^ J2[1]
                elif len(ids2) == 1:
                    J1 = ids2[0]
                    J1_K2 = K2.ideal(p, K2(str(J1[0][1].to_sage(names=K_names)))) ^ J1[1]
                    J2_K2 = K2.ideal(1)
                else:
                    raise Exception("Wrong length of ids2")
                #print "J2_K2 * J1_K2: ", (J2_K2 * J1_K2).gens_two()
                #print "I_st_K2: ", (I_st_K2).gens_two()

                assert J2_K2 * J1_K2 == I_st_K2
                #*******************************************************************

            elts = K_primes[p]['elements']
            auts = K_primes[p]['aut']
            typs =  K_primes[p]['type']

            #print(f"I_up = {I_up}")

            J,e = I_up[0]
            g = len(I_up)
            vec = [0] * (len(data['elements']) * g)
            #print(f"-> g = {g}, len(data['elements']) = {len(data['elements'])}")
            #print(f"-> len(elts) = {len(elts)}")
            for I,e in I_up:
                #print(f"I[1] = {I[1].to_sage(names=K_names)}, e = {e}")
                assert check_ideal_aut(d_K, I)
                alpha = I[1]
                found = False
                if alpha in elts:
                    vec[elts.index(alpha)] = e #1
                    found = True
                elif USE_EXCEPTIONAL_PRIMES and ( Mod(field.field(d_st).idx(), p) == 0 or Mod(field.field(d_K).idx(), p) == 0 ):
                    K_sage = NumberField([x^2 - di for di in d_K], names=K_names)
                    i = find_prime_ideal(K_sage, p, alpha, elts, K_names)
                    if i != None:
                        vec[i] = e
                        found = True
                if not found:
                    elts.append(alpha)
                    auts.append(I[2][0])
                    typs.append(I[2][1])
                    vec[len(elts)-1] = e #1
            K_st_primes[p]['powers'].append(vec)

    if parent_id == None:
        for p,data in K_primes.items():
            for i in range(len(data['elements'])):
                bits = (2^i).digits(2)
                vec = bits + [0]*(2^(len(d_K))-len(bits))
                data['powers'].append(vec)

    tree[K_id]['primes'] = K_primes

    return tree

def field_in_BV_format(tree, K_id, names=None):
    res = "Field : " + K_id.split(" < ")[0] + "\n"
    primes = sorted(tree[K_id]['primes'].items(), key=lambda itm: itm[0])
    for i_p in range(len(primes)):
        p,vals = primes[i_p]
        for i_e in range(len(vals['elements'])):
            elem = vals['elements'][i_e].to_sage(names=names, quotient=False)
            if elem == 0:
                e = str(p)
            else:
                e = elem
            pws = "[ " + ", ".join([str(i) for i in vals['powers'][i_e]]) + " ]"
            res += "{}\n{}\n{}\n".format(p, e, pws)
    return res

def conv_tree_to_BV_format_step(tree, L_id, names=None):
    res = ""
    for K_id in tree[L_id]['subfields']:
        res += field_in_BV_format(tree, K_id, names=names)
    res += field_in_BV_format(tree, L_id, names=names)
    return res

def conv_tree_to_BV_format(tree, L_id, names=None):
    res = {}
    if len(tree[L_id]['subfields']) != 0:
        res[L_id.split(" < ")[0] + ".christine"] = conv_tree_to_BV_format_step(tree, L_id, names=names)
        for K_id in tree[L_id]['subfields']:
            res.update(conv_tree_to_BV_format(tree, K_id, names=names))
    return res

def clear_dir(folder, exclude=[".gitkeep"]):
    import os, shutil
    for filename in os.listdir(folder):
        if filename in exclude:
            continue
        file_path = os.path.join(folder, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
        except Exception as e:
            print('Failed to delete %s. Reason: %s' % (file_path, e))

def clean_trees(folder="trees"):
    clear_dir(folder)

def save_tree(tree, field_id, names=None, folder = "trees/"):
    files = conv_tree_to_BV_format(tree, field_id, names=names)
    for filename, content in files.items():
        #print(filename)
        #print(content)
        file = open(folder + filename, "w")
        file.write(content)
        file.close()

def save_food(food, folder = "trees/"):
    open(f"{folder}/food.txt", "w").write(str(food))

def get_food(folder = "trees/"):
    if os.path.exists(f"{folder}/food.txt"):
        food = eval(open("trees/food.txt").read())
        res = {}
        for k,v in food.items():
            res[k] = ZZ(v)
        food = res
    else:
        food = None
    return food
