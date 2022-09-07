
#d = [17, 29, 37, 41]
#food = {'a':-11, 'b':-19, 'c':-23, 'd':-31}
food = [-11, -19, -23, -31, -43]
#dparent = [-11, -19, -31]
L.<a,b,c,d,e> = NumberField([x^2 - di for di in food])
#Lparent.<a,b,c> = NumberField([x^2 - di for di in dparent])
#Lext.<c> = L.extension(x^2-dparent[2])
Labs.<z> = L.absolute_field()
min_poly = z.absolute_minpoly()




#food ={'a':-11, 'c':-19}
#gens_list = list(sorted(food.keys()))
#variables = [ var(str(gens_list[i])) for i in range(len(gens_list)) ]



#print('gens_list', gens_list, 'food', food, 'variables', variables, type(variables[1]))
#print(gens_list[0], type(var(gens_list[0])))
#res = [1,['a']]
#K = NumberField(x^2 - d[0], names=gens_list[0])

#print('K(gens[0])', K(gens_list[0]), type(gens_list[0]))

#Field_tmp = K
#field_counter = 1
#for gen in gens_list[1:]:
#	print(gen)
#	L = K.extension(x^2 - d[field_counter], names = gen)
#	K = L
#	field_counter+=1
#print(L)


prime = 3
prime_bound = 2000

while prime < prime_bound:
	F = GF(prime)
	Fx = PolynomialRing(F, 'x')
	decomp = (Fx(min_poly)).factor()
	if (len(decomp) == 2**len(food)):
		print(prime)
	prime = next_prime(prime)


