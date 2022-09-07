

food = {'a':-11, 'b':-19, 'c':-23, 'd':-31}
d = list(food.values())
gens = list(sorted(food.keys()))
K = NumberField(x^2 - d[0], names=gens[0])
field_counter = 1
for gen in gens[1:]:
	L = K.extension(x^2 - d[field_counter], names = gen)
	K = L
	field_counter+=1

NormBound = 10**6
randBound = 10**11
INorm = 0
while INorm<NormBound:
#OK = K.maximal_order()
	gen1 = K(ZZ.random_element(-randBound,randBound))
	#gen2 = K(ZZ.random_element(-randBound,randBound))
	for el in gens:
		gen1+=K(ZZ.random_element(-randBound,randBound))*K(el)
	#gen2+=K(ZZ.random_element(-randBound,randBound))*K(el)

	I = K.ideal(gen1, ZZ.random_element(-randBound,randBound))
	I_hnf_pari = I.pari_hnf()

	I_hnf = matrix(ZZ,[I_hnf_pari[i] for i in range(I_hnf_pari.nrows())])
	det_Ihnf = I_hnf.determinant()

	INorm = I.absolute_norm()

print(det_Ihnf, INorm)

Kabs.<x> = K.absolute_field()
minpoly = x.absolute_minpoly()
pariK = pari.nfinit(minpoly)
print('running bnfinit...')
pariK_bnf = pari.bnfinit(minpoly)

for p in range(3,400):
	if not is_prime(p): continue
	if (INorm % p == 0):
		print('prime:', p)
		#vals = []
		#for p_i in K.primes_above(p):
		#	v = I.valuation(p_i)
		#	vals.append([v])
		vals_pari = []
		pabove = pari.idealprimedec(pariK, p)
		for i in range(len(pabove)):
			v = pari.idealval(pariK, I_hnf_pari, pabove[i])
			while v>0:
				I_hnf_pari = pari.idealmul(pariK, I_hnf_pari, pari.idealinv(pariK, pabove[i]))
				v-=v
print('running bnfisprincipal...')
isprincipal = pari.bnfisprincipal(pariK_bnf,I_hnf_pari, flag=1)
print(isprincipal)


