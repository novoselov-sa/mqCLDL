import nprofile
import trees
import prime_decomp
import polynomial_ring
import sys

# Note: don't use labels d,e,f simultaniously. 
# food = {'a': 5, 'b': 79, 'c': 89} # |Cl_K| = 96 (Sage)
#food = {'a': 5, 'b': 79, 'c': 89, 'd': 97, 'e': 101, 'f': 109}
#food = {'a': 5, 'b': 79, 'c': 89, 'd': 97, 'e': 101} # result: C2 x C2 x C2 x C2 x C2 x C2 x C2 x C2 x C4 x C4 x C4 x C24 x C48 x C48 x C96 x C96
#food = {"a": -3, "b": -11, "c": -71} # |Cl_K| = 56 (Sage), norm_bound = 37 (GM)
#food = {'a': -7, 'b': -11, 'c': -43} # |Cl_K| = 54, norm_bound = 25 (GM)
#food = {'a': -7, 'b': -11, 'c': -67} # |Cl_K| = 66
#food = {'a': -7, 'b': -19, 'c': -43} # |Cl_K| = 70
food = {'a': -7, 'b': -23, 'c': -31} # |Cl_K| = 207
#food = {'a': -7, 'b': -31, 'c': -47} # |Cl_K| = 510
#food = {'a': -7, 'b': -11, 'c': -19, 'd': -23} # |Cl_K| = 23040
#food = {'a': -7, 'b': -11, 'c': -19, 'd': -31}
#food = {'a':17, 'b': 29, 'c': 37, 'd': 41}
#food = {'a':-11, 'b':-19, 'c':-23, 'd':-31} #Cl_K: C3 x C3 x C12 x C2520, norm_bound = 577 (GM)
#food = {'a':-11, 'b':-19, 'c':-23, 'd':-31, 'e':-43} #Cl_K: C2 x C2 x C4 x C4 x C12 x C48 x C96 x C3360 x C3360 x C443520
#food = {'a':-11, 'b':-19, 'c':-23, 'd':-31, 'e':-43, 'g':-47} # HAS NOT TERMINATED
#food = {'a': -11, 'b': -19, 'c': -31}
#food ={'a':-11, 'b':-19}

#food = {'a': 7*11, 'b': 23*43}
#food = {'a': 3*5, 'b': 7*13, 'c': 3*13}

#food = {'a':7, 'b': 79, 'c': 223, 'd': 229}
#food = {'a': 1297, 'b': 4759, 'c': 7057, 'd': 8761}
#food = {'a': 1297, 'b': 7057, 'c': 8761} # h_1297 = 11, h_4759 = 13, h_7057 = 21
#food = {'a': 7057, 'b': 8101, 'c': 8761}
#food = {'a': 229, 'b': 257, 'c': 359, 'd': 401}
#food = {'a': -11, 'b': -19, 'c': -23}
#food = {'a': 5, 'b': 229, 'c': 257}

#food = {"a": 5, "b": 13, "c": 17}
#food = {"a": 5, "b": 13, "c": 17, "d": 29}
#food = {"a": 5, "b": 13, "c": 17, "d": 29, "e": 37} #norm_bound = 2209 (GM)
#food = {"a": 5, "b": 13, "c": 17, "d": 29, "e": 37, "g": 41}
#food = {"a": 5, "b": 13, "c": 17, "d": 29, "e": 37, "g": 41, "h": 53}

#food = {"a": -3, "b": -5, "c": -7, "d": -11} # GRHBound(K) = 455 (Magma), GM bound = 169
#food = {"a": -3, "b": -5, "c": -7, "d": -11, "e": -13, "g": -17, "h": -23}
#food = {"a": -3, "b": -5, "c": -7, "d": -11, "e": -13, "g": -17}

#food = {"a": -3, "b": -7, "c": -11}
#food = {"a": -3, "b": -7, "c": -11, "d": -19}
#food = {"a": -3, "b": -7, "c": -11, "d": -19, "e": -23}
#food = {"a": -3, "b": -7, "c": -11, "d": -19, "e": -23, "g": -31}

#food = {'a': -7, 'b': -11, 'c': -19}
#food = {'a': -7, 'b': -11, 'c': -19, 'd': -23}
#food = {"a": -3, "b": -7, "c": -11, "d": -19, "e": -23, "g": -31}

if len(sys.argv) > 1:
    l = "abcdeghijklmnopqrstuvwxyz"
    d = [ZZ(di) for di in sys.argv[1:]]
    d = tuple(sorted(d, key=lambda di: abs(di)))
    food = {}
    for i in range(len(d)):
        food[l[i]] = d[i]
#print(f"food = {food}")

# Additional primes to factor base. For DLOG computation add here primes from the factorization of norm of target ideal.
PRIMES = []

d = tuple(sorted(food.values(), key=lambda di: abs(di)))
gens = list(sorted(food.keys()))
working_folder="trees/"

# set this variable to True if you want to compute norm bound using Grenié-Molteni algorithm.
COMPUTE_GM_NBOUND = False #(len(d) <= 5)
COMPUTE_GM_NBOUND_QUAD = True

print("Generating trees for field: {})".format(d))

norm_bound = trees.compute_bound(d)
print("Bach bound =", norm_bound)

# compute norm bound using Grenié-Molteni algorithm
if COMPUTE_GM_NBOUND: # computation may be very slow 
    gp.read("norm_bounds/GM_bounds.gp")
    K = NumberField([x^2 - di for di in food.values()])
    GM_bound = ZZ(gp.GRHoptimize(K.pari_nf()))
    print("GM bound =", GM_bound)

# compute maximal bound for quadratic subfields.
if COMPUTE_GM_NBOUND_QUAD:
    gp.read("norm_bounds/GM_bounds.gp")
    I = [1]
    for di in d:
        I += [di*i for i in I]
    GM_bound_quad = 0
    for di in I:
        if di == 1:
            continue
        K = QuadraticField(di)
        b = ZZ(gp.GRHoptimize(K.pari_nf()))
        print(f"GM bound for D = {di}: {b}")
        GM_bound_quad = max(b, GM_bound_quad)
    print("GM bound (quadr.subf.) =", GM_bound_quad)

# comment out the following line if you want to use Bach bound
norm_bound = GM_bound_quad + 130
# norm_bound = GM_bound # bound computed using algorithm of Grenié-Molteni

print("Used bound =", norm_bound)
tr = trees.build_tree(d, {}, norm_bound=norm_bound, names=gens, append_primes=PRIMES)
print("-> done")

print("Removing existing trees ...")
trees.clean_trees(folder=working_folder)
print("-> done")

print("Saving new trees ...")
trees.save_tree(tr, "_".join(gens), names=gens, folder=working_folder)
trees.save_food(food, folder=working_folder)
print("-> done")

nprofile.output([trees, prime_decomp, polynomial_ring])
