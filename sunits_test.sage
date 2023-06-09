import trees
import clgp
import fb
import sunits
import field
import field_compact

food = trees.get_food()
assert food != None, "Factor base not found!"

d_base = sorted(food.values(), key=lambda di: abs(di))

for n in range(len(d_base)):
    d = tuple(d_base[:n+1])
    d_parent = tuple(d_base[:n+2])
    if d == d_parent:
        d_parent = ()

    print(f"d = {d}")
    print(f"d_parent = {d_parent}")

    KC = field_compact.field_compact(d)

    K = field.field(d)
    K_sage = K.sage(food)

    SU = sunits.sunits(d, d_parent=d_parent)
    FB = SU.factor_base()
    M = SU.relations_matrix()

    print("S-units:")
    for i in range(len(SU.gens())):
        su = SU.gens()[i]
        print(su.to_sage(K_sage), "norm:", su.absnorm())
        assert fb.is_smooth(su.absnorm(), FB), f"Wrong S-unit: {su.to_sage(K_sage)}, non-smooth!"

        for j in range(len(FB)):
            I = FB[j].to_sage()
            assert M[i][j] == su.valuation(I), f"I = {(FB[j].prime, FB[j].elts)}, v_I(su) = {su.valuation(I)} != {M[i][j]}"

        try:
            su_e = su.compute_roots(ideal_sqrt=True)
            print(f"-> explicit form: {su_e.to_sage(K_sage)}")
        except Exception as err:
            print(f"-> can't compute explicit form with compute_roots: {err}\n")
            raise Exception("fail")

        print("-> testing reduction mod S-units for field elements ...")
        a = K.random()
        su_e = su_e.evaluate()
        b = a * su_e
        c,s = b.reduce_mod_sunits(SU)
        I1 = K_sage.ideal(b.to_sage(K_sage))
        I2 = K_sage.ideal((c * s.evaluate()).to_sage(K_sage))
        assert I1 == I2, f"{b.to_sage(K_sage)} O_K != {(c * s.evaluate()).to_sage(K_sage)} O_K"
        assert c.bit_size() <= b.bit_size()

        print("-> testing reduction mod S-units for field elements (power products repr.) ...")
        a = KC.random()
        b = a * su_e
        c,s = b.reduce_mod_sunits(SU)
        I1 = K_sage.ideal(b.evaluate().to_sage(K_sage))
        I2 = K_sage.ideal((c * s).evaluate().to_sage(K_sage))
        assert I1 == I2, f"{b.to_sage(K_sage)} O_K != {(c * s).evaluate().to_sage(K_sage)} O_K"
        assert c.bit_size() <= b.bit_size()

        print("-> testing multiplication mod S-units ...")
        a = K.random()
        b = K.random() * su_e
        c,s = a.mul_mod_sunits(b, SU)
        I1 = K_sage.ideal((a*b).to_sage(K_sage))
        I2 = K_sage.ideal((c*s.evaluate()).to_sage(K_sage))
        assert I1 == I2, f"{(a*b).to_sage(K_sage)} O_K != {(c*s.evaluate()).to_sage(K_sage)} O_K"

        print("-> testing exponentiation mod S-units ...")
        a = K.random() * su_e
        e = ZZ.random_element(-4,4)
        b,s = a.pow_mod_sunits(e, SU)
        print(f"e = {e}")

        I1 = K_sage.ideal((a^e).to_sage(K_sage))
        I2 = K_sage.ideal((b*s.evaluate()).to_sage(K_sage))
        assert I1 == I2, f"{(a^e).to_sage(K_sage)} O_K != {(b*s.evaluate()).to_sage(K_sage)} O_K"

        print("-> testing removing of linear dependence mod S-units ...")
        a = K.random() * su_e
        b = K.random() * su_e^2
        e_a = ZZ.random_element(-4,4)
        e_b = ZZ.random_element(-4,4)
        c = KC((a,b), (e_a, e_b))
        r,s = c.reduce_lin_dep(mod_units=True, mod_sunits=True, SU=SU, rollback=True)

        #print(f"c = {c.to_sage(K_sage)}")
        #print(f"c.evaluate() = {c.evaluate().to_sage(K_sage)}")
        #print(f"r = {r.to_sage(K_sage)}")
        #print(f"s = {s.to_sage(K_sage)}")
        #print(f"(r*s).evaluate() = {(r*s).evaluate().to_sage(K_sage)}")

        I1 = K_sage.ideal(c.evaluate().to_sage(K_sage))
        I2 = K_sage.ideal((r*s).evaluate().to_sage(K_sage))
        assert I1 == I2, f"{c.evaluate().to_sage(K_sage)} O_K != {(r*s).evaluate().to_sage(K_sage)} O_K"

        print()
