import ideal
import relations
import fb
import trees
import clgp

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

proof.all(False)

#pari.allocatemem(1024*1024*1024)
#print(f"Pari stack size: {pari.stacksize()}")

trees_food = trees.get_food()
if food != None:
    assert food == trees_food, f"Run trees generation and relation computation first! Trees are generated for {trees_food} != {food}."

bound = prod(di.abs() for di in d)
I = ideal.random(d, bound=bound)
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
I_hnf = I.hnf(K)
print(f"Target ideal (hnf): {I_hnf}")
I_sage = I.to_sage(K)
print(f"I_sage.norm = {I_sage.absolute_norm()}")
assert I_sage.absolute_norm() == I.q

file = relations.convert_letters(d, seed=1, food=food)
FB = fb.load_primes(tuple(d), file, food=food)
M = relations.load_matrix(d)
B,U,V=M.smith_form()
V_inv = V^(-1)

assert fb.is_smooth(I.q, FB), f"Ideal is not smooth. Run class group generation including all factors of in the norm of ideal!"

print(f"CL_K: {[B[i,i] for i in range(B.nrows())]}")

w = walltime()
#print(f"FB = {FB}")
e = fb.factor(K, I_sage, FB)
print(f"e = {e}")

e_g = clgp.primes_to_gens(e, V, B)
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

e_rev = clgp.gens_to_primes(e_g, V_inv)
print(f"e_rev = {e_rev}")
I2 = fb.prod_ideal(K, FB, e_rev)
CL_K = K.class_group(proof=False)
assert CL_K(I_sage) == CL_K(I2)

#print([CL_K.gens()[i].order() for i in range(len(CL_K.gens()))])
I3 = prod([CL_K.gens()[i]^cl_sage[i] for i in range(len(CL_K.gens()))])
assert CL_K(I_sage) == CL_K(I3)
