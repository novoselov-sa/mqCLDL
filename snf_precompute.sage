import relations
import trees

print(pari.allocatemem(610*1024^3))

food = trees.get_food()
d = list(food.values())
d.sort(key=lambda v: abs(v))
d = tuple(d)

def snf_precompute(d):
    print(f"d = {d}")
    print("loading matrix of relations ... ", end="", flush=True)
    try:
        M = relations.load_matrix(d)
        print("[ok]")
    except:
        M = None
        print("[error]")

    if M != None:
        print("computing SNF ... ", end="", flush=True)
        t = walltime()
        try:
            B = relations.load_matrix(d, suffix="_snf_B")
            U = relations.load_matrix(d, suffix="_snf_U")
            V = relations.load_matrix(d, suffix="_snf_V")
            print(f"[ok] [cached] {walltime(t)} sec.")
        except:
            B,U,V = M.smith_form()
            relations.save_matrix(d, B, suffix="_snf_B")
            relations.save_matrix(d, U, suffix="_snf_U")
            relations.save_matrix(d, V, suffix="_snf_V")
            print(f"[ok] {walltime(t)} sec.")
        
        print("computing inverses (V) ... ", end="", flush=True)
        t = walltime()
        try:
            V_inv = relations.load_matrix(d, suffix="_snf_V_inv")
            print(f"[ok] [cached] {walltime(t)} sec.")
        except:
            V_inv = V^(-1)
            relations.save_matrix(d, V_inv, suffix="_snf_V_inv")
            print(f"[ok] {walltime(t)} sec.")
            
        print("computing inverses (U) ... ", end="", flush=True)
        t = walltime()
        try:
            U_inv = relations.load_matrix(d, suffix="_snf_U_inv")
            print(f"[ok] [cached] {walltime(t)} sec.")
        except:
            U_inv = U^(-1)
            relations.save_matrix(d, U_inv, suffix="_snf_U_inv")
            print(f"[ok] {walltime(t)} sec.")
        print()
    
    if len(d) == 1:
        return

    d_s  = d[:-1]
    d_t = d[:-2] + d[-1:]
    d_st = d[:-2] + (d[-2]*d[-1],)
    snf_precompute(d_s)
    snf_precompute(d_t)
    snf_precompute(d_st)

snf_precompute(d)