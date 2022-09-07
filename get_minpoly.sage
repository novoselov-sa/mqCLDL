
def mquad_disc(d):
	n = len(d)
	m4 = True
	m3 = False
	m2 = False
	P = []
	res = 1
	for di in d:
		if Mod(di, 4) != 1:
			m4 = False
		if Mod(di, 4) == 3:
			m3 = True
		for p,e in factor(di):
			if p == -1 or p in P:
				continue
			if p == 2:
				m2 = True
			res = res * p
	if m4:
		r = 0
	elif ((not m2) and m3) or (m2 and m4):
		r = 2
	else:
		r = 3
	res = (2^r * res)^(2^(n-1))
	return res

def glob_idx(d_K):
	names = [chr(ord('a')+i) for i in range(len(d_K))]
	K = NumberField([x^2 - d_i for d_i in d_K], names=names)
	disc_poly = K.absolute_polynomial().discriminant()
	mquad_disc(d_K)
	idx =  disc_poly / disc_field
	if len(d_K) == 1:
		return idx
	d_s = d_K[:-1]
	d_t = d_K[:-2] + (d_K[-1],)
	d_st = d_K[:-2] + (d_K[-1] * d_K[-2],)
	return lcm(lcm(lcm(idx, glob_idx(d_s)), glob_idx(d_t)), glob_idx(d_st))



food = {'a': -23, 'c': -43, 'b': -31, 'd': -47}
d = tuple(sorted(food.values(), key=lambda di: abs(di)))
gens = list(sorted(food.keys()))
names = [chr(ord('a')+i) for i in range(len(d))]

K = NumberField([x^2 - d_i for d_i in d], names=names)
minpoly= K.absolute_polynomial()
print(minpoly)
disc_poly = minpoly.discriminant()
disc_field = mquad_disc(d)
idx =  disc_poly / disc_field
g_idx = glob_idx(d)

P = Primes()
p = 11
while true:
	if (Mod(g_idx, p) != 0) and  (not Mod(disc_field, p) != 0):
		print(p)
		break
	p = P.next(p)




