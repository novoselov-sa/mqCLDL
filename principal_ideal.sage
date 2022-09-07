food = {'a':-11, 'b':-19}
d = list(food.values())
gens = list(sorted(food.keys()))
K = NumberField(x^2 - d[0], names=gens[0])
field_counter = 1
for gen in gens[1:]:
	L = K.extension(x^2 - d[field_counter], names = gen)
	K = L
	field_counter+=1

#
#	generate principal ideal
#
gen = K.random_element()
I = K.ideal(gen*gen.denominator())
print('I:', I)
I_hnf = I.pari_hnf()

SubFields = K.subfields()
print(SubFields)
print(type(SubFields[1]))

Automs = K.automorphisms()
print(Automs)

