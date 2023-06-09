import relations
import trees
import idealprod
import ring
import field
import field_compact
import ideal
import clgp

print(pari.allocatemem(610*1024^3))

SAVE_ALL_RESULTS = "experiments/dlogs"

food = trees.get_food()
d = list(food.values())
d.sort(key=lambda v: abs(v))
d = tuple(d)
d_parent = ()

# Heuristic check for equality I == g * prod P_i^(e_i)], i = 1, ..., #FB.
def ideal_eq_test(K, FB, I, e, g, norms=False, valuations=True):
    assert len(FB) == len(e), f"Wrong exponent vector e, len(e) = {len(e)} != {len(FB)}"
    if norms:
        known_factors = [I.q] + [FB[i].prime for i in range(len(FB))]
        nrm = prod([K.ideal(FB[i].prime, K(FB[i].elts)).absolute_norm()^ZZ(e[i]) for i in range(len(e))])
        if abs(I.q) != abs(g.absnorm(known_factors=known_factors) * nrm): # compare upto units
            print("\t[fail: norms]")
            return False

    if valuations:
        print("[to_sage:", end="", flush=True)
        t = walltime()
        I_sage = I.to_sage(K)
        print(f"{walltime(t)} sec.]\n", end="", flush=True)
        for i in range(len(e)):
            print(f"[primeid_{i}/{len(e)-1}:", end="", flush=True)
            t = walltime()
            P = K.ideal(FB[i].prime, K(FB[i].elts))
            print(f"{walltime(t)} sec.]\n", end="", flush=True)

            print(f"[val_{i}/{len(e)-1}:", end="", flush=True)        
            if g.valuation(P) + e[i] != I_sage.valuation(P):
                print(f"fail, {walltime(t)} sec.]\n", end="", flush=True)
                return False
            print(f"ok, {walltime(t)} sec.]\n", end="", flush=True)
    return True

def dlog_verify(d, d_parent=()):
    print(f"d = {d}, parent = {d_parent}")
    K = field.field(d)
    CL = clgp.clgp(d, d_parent=d_parent)
    fld = "_".join([str(d[i]) for i in range(len(d))])
    fn = f"{SAVE_ALL_RESULTS}/{fld}_dlog"

    print("Loading input ideal ...", flush=True, end="")
    t = walltime()
    try:
        with open(f"{fn}_input", 'r') as f:
            data = f.read()
        J = eval(preparse(data))
        print(f"[ok] {walltime(t)} sec.")
    except FileNotFoundError:
        J = None
        print(f"[error]")

    if J != None:
        print(f"J = {J}")

        print("Loading ideal decomposition ... ", flush=True, end="")
        t = walltime()
        with open(fn, 'r') as f:
            data = f.read()
        G = eval(preparse(data))
        print(f"{walltime(t)} sec.")
        print(f"G = {G.to_sage(K.sage(food))}")

        print("Loading factor base ... ", flush=True, end="")
        t = walltime()
        FB = CL.factor_base()
        print(f"{walltime(t)} sec.")
        
        print("Verifying [valuations] ... ", flush=True, end="")
        #if ideal_eq_test(K.sage(food), FB, J, G.powers, G.element):
        #    print("[ok]")
        #else:
        #    print("[error]")
        assert ideal_eq_test(K.sage(food), FB, J, G.powers, G.element)        
    print()

    d_s  = d[:-1]
    d_t = d[:-2] + d[-1:]
    d_st = d[:-2] + (d[-2]*d[-1],)

    #J_s = ideal.ideals(d_s)(J.q, J.s[:-1])
    #J_t = ideal.ideals(d_t)(J.q, J.s[:-2] + J.s[-1:])
    #J_st = ideal.ideals(d_st)(J.q, J.s[:-2] + (J.s[-2]*J.s[-1],))

    if len(d_s) != 1:
        dlog_verify(d_s, d_parent=d)
    if len(d_t) != 1:
        dlog_verify(d_t, d_parent=d)
    if len(d_st) != 1:
        dlog_verify(d_st, d_parent=d)

dlog_verify(d, d_parent=d_parent)

