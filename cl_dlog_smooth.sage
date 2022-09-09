import ideal
import relations
import fb
import trees

#food = {'a': 17, 'b': 29, 'c': 37, 'd': 41}
#food = {'a':7, 'b': 11, 'c': 19, 'd': 23}
#food = {'a':7, 'b': 79, 'c': 223, 'd': 229}
#food = {'a': 1297, 'b': 4759, 'c': 7057, 'd': 8761}
#food = {'a': 4759, 'b': 8761} # h_4759 = 13, h_8761 = 27
#food = {'a': 1297, 'b': 8761, 'c': 7057} # h_1297 = 11, h_4759 = 13, h_7057 = 21
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

#food = {"a": -3, "b": -7, "c": -11, "d": -19}
#food = {"a": -3, "b": -7, "c": -11, "d": -19, "e": -23}
food = trees.get_food()

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
    print(f"[factor_over_FB] N_I = {N_I}")
    e = []
    for i in range(len(FB)):
        if Mod(N_I, FB[i].prime) != 0:
            e.append(0)
        else:
            P = K.ideal(FB[i].prime, K(FB[i].elts))
            v = ZZ(pari.idealval(K.pari_nf(), I, P.pari_prime()))
            e.append(v)
            print(f"P = {P}, v = {v}")
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
    #K_gens = K.gens()
	#print(K, K_gens)
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

def prod_ideal(K, FB, e):
    assert len(FB) == len(e), f"Wrong length of factor base or input vector, {len(FB)} != {len(e)}"
    return prod(K.ideal(FB[i].prime, K(FB[i].elts))^int(e[i]) for i in range(len(FB)))

def load_relations(d):
    fn = "relations/" + "_".join([str(di) for di in sorted(d)]) + "_relations"
    with open(fn, "r") as f:
        data = f.read()
        return matrix(eval(data))

trees_food = trees.get_food()
if food != None:
    assert food == trees_food, f"Run trees generation and relation computation first! Trees are generated for {trees_food} != {food}."

bound = prod(di.abs() for di in d)
I = random_ideal(d, bound=bound)
#I = ideal.ideals([5, 229, 257])(669289,(303651, 6337, 202003))
#I = ideal.ideals([-7, -11, -19, -23])(2557,(143, 1040, 485, 1051))
#I = ideal.ideals([5, 13, 17, 29, 37, 41])(9032339,(2001020, 3845142, 269749, 1568238, 3382698, 560462))
#I = ideal.ideals((5, 13, 17, 29, 37))(744251,(275720, 225315, 139312, 8627, 29735))
#I = ideal.ideals((5, 13, 17, 29, 37))(271879,(49480, 115130, 12857, 3796, 90974))
#I = ideal.ideals((5, 13, 17, 29, 37))(879941,(341095, 173263, 46592, 397104, 411054))
#I = ideal.ideals((-7, -31, -47))(9907,(536, 2192, 4326))

#I = ideal.ideals((-3, -7, -11, -19))(1831,(486, 121, 582, 789))
#I = ideal.ideals((-3, -7, -11, -19))(3733,(1836, 632, 1451, 1333))

print(f"Target ideal: {I}")
K = NumberField([x^2 - d[i] for i in range(len(d))], names=list(sorted(food.keys()))[:len(d)])
I_hnf = ideal_to_hnf(K, d, I)
print(f"Target ideal (hnf): {I_hnf}")
I_sage = K.ideal(pari.idealhnf(K.pari_nf(), I_hnf.hermite_form(include_zero_rows=False).transpose()))
print(f"I_sage.norm = {I_sage.absolute_norm()}")
assert I_sage.absolute_norm() == I.q

file = relations.convert_letters(d, seed=1, food=food)
FB = fb.load_primes(tuple(d), file, food=food)
M = load_relations(d)
B,U,V=M.smith_form()
V_inv = V^(-1)

primes = [FB[i].prime for i in range(len(FB))]
flag = True
for p,e in factor(I.q):
    if not p in primes:
        print(f"p = {p} | I.q is not in factor base, run class group generation including this prime")
        flag = False
assert flag

print(f"CL_K: {[B[i,i] for i in range(B.nrows())]}")

w = walltime()
#print(f"FB = {FB}")
e = factor_over_FB(K, I_sage, FB)
print(f"e = {e}")

e_g = prod_primes_to_gens(e, V, B)
print(f"cl_dlog: {e_g}, time: {walltime(w)}", flush=True)

print("Checking using Sage ...") # Unreachable for n >= 6
if len(d) >= 6:
    print("-> DISABLED")
    exit(0)

w = walltime()
cl_sage = I_sage.ideal_class_log(proof=False)
cl_sage = [ZZ(i) for i in cl_sage]
print(f"-> [sage] cl_dlog = {cl_sage}, time: {walltime(w)}")
print(f"-> [sage] I.factor() = {I_sage.factor()}")

e_rev = prod_gens_to_primes(e_g, V_inv)
print(f"e_rev = {e_rev}")
I2 = prod_ideal(K, FB, e_rev)
CL_K = K.class_group(proof=False)
assert CL_K(I_sage) == CL_K(I2)

#print([CL_K.gens()[i].order() for i in range(len(CL_K.gens()))])
I3 = prod([CL_K.gens()[i]^cl_sage[i] for i in range(len(CL_K.gens()))])
assert CL_K(I_sage) == CL_K(I3)
