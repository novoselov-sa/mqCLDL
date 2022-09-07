import sys
import nprofile
#profiling = ['get_qstuff','lift_ais','get_symbols','kernel','squares','adjoin_sqrt','shortening_matrix','shorten','get_cl_gens_optimized']
import field
import fb
import units
from memoized import memoized
import char
import copy

STORE_RELATIONS = True
SHORTENNING_ENABLED = True # use shorten relations for subfields using HNF or LLL, or not.
SHORTENNING_USE_LLL = True # use LLL to shorten relations for subfields

COMPUTE_CLGP_FOR_SUBFIELDS = True

aism = 0
sigmas = {}
pointerdic = {}
ais = []
aisc = []
sqrtapprox = {}
finald = []

def set_vars(ais, sigs, pars, d):
  global aism 
  aism = ais
  global finald 
  finald = list(d)
  global sigmas
  sigmas = sigs
  global pointerdic
  pointerdic = pars
  global ais
  ais = aism*[0]
  global aisc
  aisc = aism*[0]
  global sqrtapprox
  #XXXX: If stuff goes wrong, try upping this
  R = RealField(10037)
  I = [1]
  for di in d:
    I += [di*i for i in I]
  for i in I:
    sqrtapprox[i] = (sqrt(i))

#
#	Intput: currect subfields, ais, aisc are global
# Output: lift all three pairs (ais, aisc) to the upper field
#
@memoized
def lift_ais_internal(d, seed):
	#print('seed:', seed)
  K = field.field(d)
  N = 2^len(d)
  I = [1] #WHY DO WE NEED 1 IN I ???
  for di in d:
    I += [di*i for i in I]
  k0 = K(N*[0],1) # TODO: this is 0 in K. Why do we need it???
  A = aism*[k0]
  Ac = aism*[k0]
#print 'pointerdic:', pointerdic
  for k in pointerdic.keys():
    if k in I:
      targ = I.index(k)
			#
			# there ais contains S-units for each of the three subfield
			# aisc contains their conjuagates
			# pointerdic point to the rage where these units are stored in ais (aisc)
			# for each one of the subfields
			#
      B = ais[pointerdic[k][0]:pointerdic[k][0] + pointerdic[k][1]]
      Bc = aisc[pointerdic[k][0]:pointerdic[k][0] + pointerdic[k][1]]

      for b in range(len(B)):
				
				#
				# lift the found S-units (stored as coeff-vectors [free-coeff, coeff_in_front_of d_i's]) to the larger field
				# by creating a vector (twice as long as the coeff-vector), where [coeff_in_front_of d_i's] are shifted accoring
				# to the d_i's position in the larger field
				# i.e., [c_1, c_2] from Q(sqrt(d1*d2)) is lifted to [c_1, 0, 0, c_2]
				# targ is responsible for the new position
				#
        lvec = [B[b].numer.c[0]] + (N-1)*[0]
        lvec[targ] = B[b].numer.c[1]
				#
				# there could be denomiators in the S-units (???)
				#
        aa =  K(lvec, B[b].denom)
        #try: evaluatereal(aa,d) > 0
        #except: print aa,d
        A[pointerdic[k][0] + b] = aa
				
				#
				# same as the above but for the cojuagates
				#
        lvec = [Bc[b].numer.c[0]] + (N-1)*[0]
        lvec[targ] = Bc[b].numer.c[1]
        aac =  K(lvec, Bc[b].denom)
        #try: evaluatereal(aac,d) > 0
        #except: print aac,d
        Ac[pointerdic[k][0] + b] = aac
  return A, Ac

@nprofile.profile
def lift_ais(d, seed):
  return lift_ais_internal(d, seed)

#
# Inputs:
#	-- d describes the field we're currecnly lifting to
# -- high is the size of the factor-based of the field we're currenly lifting to (not sure though)
#
# Outputs:
#
@memoized
def get_qstuff_internal(d, high, seed):
  K = field.field(d)
  N = 2^len(d)
  k0 = K(N*[0],1) # TODO: this is 0 in K. Why do we need it???
  A, Ac = lift_ais(d, seed)
  Aq = []
  Aqc = []
  Adenoms = []
  qs = []
  for a in range(len(A)):
    if A[a] == k0:
      Aq += [(high + 64)*[0]]
      Aqc += [(high + 64)*[0]]
      Adenoms += [1]
    else:
      Adenoms += [A[a].denom]
			#
			# get_qs is in field.sage
			# 
      aq = A[a].get_qs(0,(high + 64))
      aqc = Ac[a].get_qs(0,(high + 64))
			#print('A[a]:', A[a])
      #print "getting qs",a,walltime(t)
      aqn = len(aq)
      Aq += [[aq[i][0] for i in range(aqn)]]
      Aqc += [[aqc[i][0] for i in range(aqn)]]
			#print('Aq:', Aq)
      if qs == []:
            qs =  [aq[i][1] for i in range(len(aq))]
  Aq = matrix(Aq)
  Aqc = matrix(Aqc)
  return A, Ac, Aq, Aqc, Adenoms, qs

@nprofile.profile
def get_qstuff(d,high, seed):
	#print('get_qstuff seed:', seed)
  return get_qstuff_internal(d,high,seed)

@memoized
def quadfield(D, name='a'):
  return QuadraticField(D, name)

#
# returns the classgroup, the elements and factorizations of the relations and the factorbase
# input: d - the generator of quadratic extension
#				FB - factor-base of the type 'tuple' of the form (fb.prime((D))(3,1/2*a + 1/2,[1, 0])
#				aism - factor-base size
#
@memoized
def quadrelations(D, FB, aism, food = None):
  assert len(FB) != 0, "factor base should not be empty"
  #search replacement generator
  #gens =  ['a','b','c','d','e','f','g']
  if food == None:
    gens =  ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
  else:
    gens = list(sorted(food.keys()))

  K = quadfield(D, gens[0])
  OK= K.ring_of_integers()
	#print K, OK
  
  gen = []
  #FB should not be empty and have a generator hidden somewhere
  i = 0
  non_intert_exists = False
  genfinal = FB[i].elts.split(' ')[0]
  while gen == []:
    for ge in genfinal:
      if ge in gens: 
        gen = genfinal.index(ge)
        break
    i += 1
    genfinal = FB[i].elts.split(' ')[0]
  genfinal = genfinal[gen:]
#TODO: FIX THE WHILE LOOP: it's mere purpose is to get the letter of the generator in the ideal's representation
  FB2 = []
  for pr in FB:
    #replace generator
    #pre = pr.elts.replace(genfinal,'a')
    pre = pr.elts.replace(genfinal, gens[0])
    FB2 += [K.ideal(pr.prime, K(pre))]
	#
	# FB2 now is a list with elements of the form Fractional ideal (3, 1/2*a + 1/2)
	#
  
  m = len(FB2)
  SUK = K.S_units(FB2) # outputs a list of generators of the S-unit group



#print 'SUK:',  SUK, len(SUK)
  v = []

  av = []
  avc = []
  R = RealField(53)
  C = ComplexField(53)
  #for su in SUK[:-2] + [SUK[-1]]:
  for su in SUK:
    if su == -1:
      continue
    su0, su1 = list(su) #su0 - free coeff, su1 - gen's coeff
    if su0 + su1*C(sqrt(D)) > 0:
      av += [su] #FIXME: make compatible with imaginary
    else:
      av += [-su]
    if su0 - su1*C(sqrt(D)) > 0:
      avc += [K([su0,-su1])] #FIXME: make compatible with imaginary
    else:
      avc += [-K([su0,-su1])]
    vtemp = m*[0]

#
# each element su from the S-units is factorized and the resulting factors suf are searched insied the factor-base FB2
# vtemp stores the degree of a prime p from FB2 that appears in the factorization of suf
# this degree is stoted in the vtemp at the position of p in FB2
# vector v gives the relations btw. the S-units and the factor-base as
# SU[i] = \prod_j p_j^v[i,j]
#
    for suf in K.ideal(su).factor(): #TODO: exucute factor() several times???
      i = FB2.index(suf[0]) #TODO: OPTIMIZE, index may be costly
      vtemp[i] = suf[1]
    v += [vtemp]
 
  VMat = matrix(ZZ, [vvec for vvec in v])
  #print('VMat:', VMat.elementary_divisors())
  VMat_hnf = VMat.hermite_form(include_zero_rows=False)
  print(f'for D = {D} VMat_hnf:')
  print(VMat_hnf.elementary_divisors())

  if STORE_RELATIONS:
    print(f"Saving relations for field {(D,)}")
    save_matrix([D], VMat_hnf)

  #optional check that the relations hold
  #for i in range(len(av)):
  #    if not check(av[i]*OK, v[i], FB2, m):
  #        print av[i]
  #        print v[i]
  #        raise AssertionError

  #optional class group computation
  #clK = tuple([x for x in matrix(ZZ,v).hermite_form(include_zero_rows = False)[:m].elementary_divisors() if x != 1])
  clK = 0
  #Optional class group checking: keep optional
  #if not clK == K.class_group().elementary_divisors():
  #  raise AssertionError

  #put the elements into our context
	#removes the generators, stores elements from SU as pairs
  K = field.field(tuple([D]))
  for i in range(len(av)):
    #print 'av[i]', av[i]
    d = av[i].denominator()
    av[i] = K((d*av[i]).list(), d)
    avc[i] = K((d*avc[i]).list(), d)
	#print 'av[i].denominator():', d, 'av[i]:', av[i], 'avc[i]', avc[i]

  #initialize elements as quad representation
#print 'len(av):', len(av), 'pointerdic[D][1]:', pointerdic[D][1]

  assert len(av) == pointerdic[D][1]
  L = zero_matrix(len(av), aism) #len(av) >= aism - 1 (since we removed -1 from SU)
  L[:,pointerdic[D][0]:pointerdic[D][0] + len(av)] = identity_matrix(len(av))

  ais[pointerdic[D][0]:pointerdic[D][0] + len(av)] = av
  aisc[pointerdic[D][0]:pointerdic[D][0] + len(av)] = avc
	#print(ais, aisc)
  #identity matrix is for the powers of the base relations
  return clK, av, L, zero_matrix(len(av), aism), v


#Checks if stuff still factorizes correctly
def check(I, vec, FB, m):
    I2 = prod([FB[i]^vec[i] for i in range(m)])
    return I2 == I

@memoized
def relations(d, seed):
	#print('relations seed:', seed)
  n = len(d)
  N = 2^n
  K = field.field(d)
  class R:
    def __init__(f,*args):
      if len(args) == 0: raise Exception('not known how to initialize empty relations')
      # one relation as input of the form (elt, factorization vector)
      # XXX: what is the len of this?
			#print(f, len(args))
      if len(args) == 1:
        print('len(args)==1')
        g = args[0]
        # return the trivial relation
        if g in [-1,1]: # U(1) is 1; U(-1) is -1
          f.element = K((g,)+(0,)*(N-1),1)
          f.powers = (0,)*N
          f.factors = (0,)*N
          return
        if n >= 1 and g.__class__ == relations(d[:-1], seed):
          print "heuj"
          f.element = K(g.element)
          f.powers = g.powers
	  # XXX: need to do something here with the change factor base matrix (primes over other primes and all that, same for the two following)
          # XXX: maybe load something from file here
          f.factors = g.factors + g.factors
          # XXX: might be able to do without qchars,otherwise need len(factors)
          f.qchars = f.element.symbols(0, len(f.factors) + 64)
          return
        if n >= 2 and g.__class__ == relations(d[:-2] + d[-1:], seed):
          f.element = K(g.element)
          f.powers = g.powers[:N//4] + (0,)*(N//4) + g.powers[N//4:] + (0,)*(N//4)
          f.approxlog = g.approxlog[:N//4]*2 + g.approxlog[N//4:]*2
          return
        if n >= 2 and g.__class__ == relations(d[:-2] + (d[-2]*d[-1],), seed):
          f.elements = K(g.element)
          f.powers = g.powers[:N//4] + (0,)*(N//2) + g.powers[N//4:]
          f.approxlog = g.approxlog + g.approxlog[N//4:] + g.approxlog[:N//4]
          return
      if len(args) == 3: # XXX: trusting caller to provide suitable values
				#print('len(args)==3')
        f.powers, f.conj, f.factors = args
        return
      raise Exception('not known how to initialize this format')
    def __repr__(f):
#      return 'relations.relations(%s)(%s,%s,%s,%s)' % (d,f.element,f.powers,f.factors,f.qchars)
      return 'relations.relations(%s)(%s,%s,%s)' % (d,f.powers,f.conj,f.factors)
    def __eq__(self, other):
#      return self.element==other.element and self.factors==other.factors and self.powers==other.powers and self.conj==other.conj
      if vector(self.factors).is_zero() and vector(other.factors).is_zero(): return False
      else: return self.factors==other.factors
    def __hash__(self):
#      return hash(tuple(self.factors) + tuple(self.powers) + tuple(self.conj))
      return hash(tuple(self.factors))
    def __mul__(f,g):
      hpowers = tuple(f.powers[j] + g.powers[j] for j in range(len(f.powers)))
      hconj = tuple(f.conj[j] + g.conj[j] for j in range(len(f.conj)))
      hfactors = tuple(f.factors[j] + g.factors[j] for j in range(len(f.factors)))
      return R(hpowers,hconj,hfactors)
    def __div__(f,g):
      hpowers = tuple(f.powers[j] - g.powers[j] for j in range(len(f.powers)))
      hconj = tuple(f.conj[j] - g.conj[j] for j in range(len(f.conj)))
      hfactors = tuple(f.factors[j] - g.factors[j] for j in range(len(f.factors)))
      return R(hpowers,hconj,hfactors)
    def sqrt(f):
      hpowers = tuple(f.powers[j]/2 for j in range(len(f.powers)))
      hconjs = tuple(f.conj[j]/2 for j in range(len(f.conj)))
      hfactors = tuple(f.factors[j]/2 for j in range(len(f.factors)))
      return R(hpowers,hconjs,hfactors)
    def symbols_log(f, Aq, Aqc, Adenoms, qs):
      s = char.altfinal(*get_qinfo(f, Aq, Aqc, Adenoms, qs))
      if 0 in s:
        raise Exception('zero symbol')
      return tuple(0 if si == 1 else 1 for si in s)
  return R

def subsetprod(d,g,e,seed):
  R = relations(d,seed)
  for gj in g:
    if gj.__class__ != R:
      raise Exception('%s not in relations(%s)' % (gj,d))
  M = len(g[0].factors)
  hpowers = [tuple(sum(gj.powers[k] for gj,eij in zip(g,ei) if eij) for k in range(aism)) for ei in e]
  hconj = [tuple(sum(gj.conj[k] for gj,eij in zip(g,ei) if eij) for k in range(aism)) for ei in e]
  hfactors = [tuple(sum(gj.factors[k] for gj,eij in zip(g,ei) if eij) for k in range(M)) for ei in e]
  return tuple(R(p,c,f) for p,c,f in zip(hpowers,hconj,hfactors))

def powerprod(d,g,e, seed):
  R = relations(d, seed)
  for gj in g:
    if gj.__class__ != R:
      raise Exception('%s not in relations(%s)' % (gj,d))
  n = len(d)
  N = 2^n
  hfactors = [tuple(sum(eij*gj.factors[k] for gj,eij in zip(g,ei) if eij) for k in range(len(gj.factors))) for ei in e]
  hpowers = [tuple(sum(eij*gj.powers[k] for gj,eij in zip(g,ei) if eij) for k in range(len(gj.powers))) for ei in e]
  hconjs = [tuple(sum(eij*gj.conj[k] for gj,eij in zip(g,ei) if eij) for k in range(len(gj.conj))) for ei in e]
  return tuple(R(p,c,f) for p,c,f in zip(hpowers,hconjs,hfactors))

@nprofile.profile
def kernel(n,M):
  return M.left_kernel().basis_matrix().rows()

@nprofile.profile
def squares(n,d,S,e,seed):
  return subsetprod(d,S,e,seed)

@nprofile.profile
def get_symbols(d, n, S, seed):
	#print('get_symbols seed:', seed)
  A, Ac, Aq, Aqc, Adenoms, qs = get_qstuff(d, n, seed)
  return matrix(GF(2), [u.symbols_log(Aq, Aqc, Adenoms, qs) for u in S])

@nprofile.profile
def adjoin_sqrt(n,S,E):
  return S + [sqrt(Ei) for Ei in E]

# return matrix that turns vectors into short nonzero vectors, this function doesn't work, but wants to use LLL
@nprofile.profile
def shortening_matrix_alt(n,S):
  H = len(S)
  relmat = matrix(ZZ,[[H*si for si in Si.factors] for Si in S])
  I = identity_matrix(H)
  new = I.augment(relmat)
  newLLL= new.LLL()
  #print newLLL.str()
  M = matrix([newi[:H] for newi in newLLL if not newi[H:].is_zero()])
  return M

# return matrix that turns vectors into short nonzero vectors, I am doing it with HNF (will this get too big?)
@nprofile.profile
def shortening_matrix(n,S):
  relmat = matrix(ZZ,[Si.factors for Si in S])
  H,M = relmat.hermite_form(transformation=True, include_zero_rows=False)
  return H,M

@nprofile.profile
def shortening_matrix_yet_another(n,S):
  relmat = matrix(ZZ,[Si.factors for Si in S])
  #H,M = relmat.hermite_form(transformation=True, include_zero_rows=False)
  H,M = relmat.LLL(transformation=True)
  # FIXME: remove zero rows from the matrix H with nessary changes in transformation matrix M
  return H,M

@nprofile.profile
def shorten(n,d,S,L,seed):
  return powerprod(d,S,L,seed)

#for the purpose of loading from files
@memoized
def convert_letters(d, seed, food = None):
  if food == None:
    food = {"a": 5, "b": 13, "c": 17, "d": 29, "e": 37, "f": 41, "g": 53, "h": 61}
  revfood = dict((abs(k),v) for v,k in food.iteritems())
	#print revfood
  #revfood = {5 : "a", 13 : "b", 17 : "c", 29: "d", 37 : "e", 41 : "f", 53 : "g", 61 : "h"}
  d2 = [list(zip(*list(abs(di).factor())))[0] for di in d]
  filecode = '_'.join([''.join([revfood[abs(dii)] for dii in di]) for di in d2])
  return filecode

def get_qinfo(r1, Aq, Aqc, Adenoms, qs):
     vec = []
     vecc = []
     aq = []
     aqd = []
     aqc = []
     aqcd = []
     for i in range(len(r1.powers)):
       if (r1.powers[i] != 0):
           vec += [r1.powers[i]]
           aq += [Aq[i]]
           aqd += [Adenoms[i]]
       if (r1.conj[i] != 0):
           vecc += [r1.conj[i]]
           aqc += [Aqc[i]]
           aqcd += [Adenoms[i]]
     vec = vector(vec)
     vecc = vector(vecc)
     denom = max(vec.denominator(), vecc.denominator())
     return aq, aqd, aqc, aqcd, qs, denom*vec, denom*vecc, denom

def save_matrix(d, M, zero_rows=False):
  fn = "relations/" + "_".join([str(di) for di in sorted(d)]) + "_relations"
  r = M.rows()
  with open(fn, "w") as f:
    f.write("[")
    for i in range(len(r)):
      if r[i].is_zero():
        continue
      f.write(str(r[i]))
      if i != len(r) - 1:
        f.write(",\n")
    f.write("]")

def load_matrix(d):
  fn = "relations/" + "_".join([str(di) for di in sorted(d)]) + "_relations"
  with open(fn, "r") as f:
    data = f.read()
    return matrix(eval(data))

# generators of full unit group mod torsion
def relations_internal(d, file, seed, food=None):
	#print('Seed relations_internal', seed)
  n = len(d)
  if n == 0:
    return ()
  
	#print 'd:', d

  if n == 1:
    R = relations(d, seed)
    P = fb.load_primes(tuple(d), file, food=food)
    #print "loaded primes:", d, P
		#print 'tuple(P):', tuple(P)
		#print 'food:', food
		
		#
    #compute quad relations and factorbase
		#
    RQ = quadrelations(d[0], tuple(P), aism, food=food)
		#
		# now RQ contains clK (possibly 0), av (S-units in the form of pairs), L, zero_matrix(len(av), aism), v (list of relations)
		# the below concatenates the triples (i-th unit vector, 0-vector, i-th relation vector) into one large list
		# then R(*r) calls the constructor of R with the triple r
		#
		# This way S-units are stored via the relation vectors RQ.factors: S[i]  = \prod_j FB[j]^(RQ[i].factors[j])
		# RQ.powers storex indices of the resulting elements in the list SUK = K.S_units(FB)
		#
    RQ = [R(*r) for r in zip(*RQ[2:5])]
		#print('RQ:', RQ)
    return RQ 

  N = 2^n
  K = field.field(d)
  R = relations(d, seed)
	#print('R:', R)
  file = convert_letters(d,seed, food=food)

	#print "d", d
  R1 = relations_internal(d[:-1], file, seed, food=food)
	#print "R1:", R1
  FB1 = fb.load_primes(tuple(d[:-1]), file, food=food)
	#print('FB1 file', d, file, len(FB1) )
  mr1 = len(R1)
  m = len(FB1)

  R2 = relations_internal(d[:-2] + d[-1:], file, seed, food=food)
  mr2 = len(R2)
  FB2 = fb.load_primes(tuple(d[:-2] + d[-1:]), file, food=food)

  R3 = relations_internal(d[:-2] + (d[-2]*d[-1],), file, seed, food=food)
  mr3 = len(R3)
  FB3 = fb.load_primes(tuple(d[:-2]) + (d[-2]*d[-1],), file, food=food)

  print(f"Lifting into field {d}")
  #lift field 1
	#
	# lifting is done via going through the collected relations R1
	# a new vector newv is created. It stores indices of the primes that take part in the relation
	# lifted to the factor-base of the above field
	# Example: if (5) = [1,1] is involved in a simple relation 5 =5 (i.e, 5 is in the units of the subfield)
	# which is described by (0,0,1,0...) relation vector. Such relation will be lifted to (0,0,1,1,0,..) relation
	# Putting the non-zero indicies properly is what is done below
	#
	
  R1new = []
  for r1 in R1:
    currentp = 0
    #print('r1:', r1)
    newv = []
    assert len(FB1)==len(r1.factors), f'{len(FB1)}!={len(r1.factors)}'
    for i in range(len(FB1)):
      #print('newv:', newv)
      #print('current p:', FB1[i].prime, 'FB1[i].powers:', FB1[i].powers, 'r1.factors[i]', r1.factors[i])
      if currentp != FB1[i].prime:
            newv += [r1.factors[i]*j for j in FB1[i].powers]
            currentp = FB1[i].prime
      else:
         c = len(FB1[i].powers)
         newv[-c:] = [newv[-c:][h] + r1.factors[i]*FB1[i].powers[h] for h in range(c)]
    #print('newv:', newv)
    R1new += [R(r1.powers, r1.conj, newv)]

  #lift field 2
  R2new = []
  for r2 in R2:
    assert len(FB2)==len(r2.factors), f'{len(FB2)}!={len(r2.factors)}'
    currentp = 0
    newv = []
    for i in range(len(FB2)):
      if currentp != FB2[i].prime: 
            newv += [r2.factors[i]*j for j in FB2[i].powers]
            currentp = FB2[i].prime
      else:
         c2 = len(FB2[i].powers)
         newv[-c2:] = [newv[-c2:][h] + r2.factors[i]*FB2[i].powers[h] for h in range(c2)]
    R2new += [R(r2.powers, r2.conj, newv)]

  #lift field 3: the lifts are already conjugated in the loaded primes
  R3new = []
  for r3 in R3:
    currentp = 0
    newv = []
    for i in range(len(FB3)):
      if currentp != FB3[i].prime:
            newv += [r3.factors[i]*j for j in FB3[i].powers]
            currentp = FB3[i].prime
      else:
         c = len(FB3[i].powers)
         newv[-c:] = [newv[-c:][h] + r3.factors[i]*FB3[i].powers[h] for h in range(c)]
    conjugate =  zip(*[(ZZ((1+c3)/2)*a - ZZ((-1+c3)/2)*b,
                        ZZ((1+c3)/2)*b - ZZ((-1+c3)/2)*a) for a,b,c3 in zip(r3.powers,r3.conj,sigmas[d[-1]])])
    conjugate = list(conjugate)
    R3new += [R(conjugate[0], conjugate[1], newv)]

  #let's do a cheat check
  #for i in range(len(R3new)):
  #   #print "praying", r3.element
  #   for i in range(len(r3.powers)):
  #     if (r3.powers[i] != 0) or (r3.conj[i] != 0):
  #       print r3.powers[i],r3.conj[i], ais[i]

  #S: check relations
  #for r3 in R3new:
  #  for i in range(len(r3.powers)):
  #    if not ((r3.powers[i] != 0) or (r3.conj[i] != 0) or ((r3.conj[i] == 0) and (r3.powers[i] == 0))):
  #      print("[R3new check failed] ", r3.powers[i], r3.conj[i], ais[i])
  #      assert False, 'R3new check failed'

  #S = uniq(sorted(R1new + R2new + R3new))
  S = list(set(R1new + R2new + R3new))
  #print('R1new.factors:', [Ri.factors for Ri in R1new])
	#print('R2new.factors:', [Ri.factors for Ri in R2new])
	#print('R3new.factors:', [Ri.factors for Ri in R3new])
  print("Preparing ais")
  sys.stdout.flush()
	#
	#	len(newv) is the legnth of the last newv-vector coming from the third lift
	# it should (?) be of length |factor-base of the field above|
 	#	TODO: EXPLAIN THIS!!!
  A, Ac, Aq, Aqc, Adenoms, qs = get_qstuff(d, len(newv), seed)

  print("Creating Symbols")
  #create symbols matrix
  #M = get_symbols(S, Aq, Aqc, Adenoms, qs)
	#
	# This should be the A-matrix from Alg. 2 from the original paper
	# https://msp.org/obs/2019/2-1/obs-v2-n1-p07-p.pdf
	#
  M = get_symbols(d, len(newv), S, seed)
	#print M.str()

  print("Computing Kernel")
  sys.stdout.flush()
  #compute kernel
  e = kernel(n, M)
  #print("kernel", e)

  print("Computing group of squares")
  sys.stdout.flush()
  #compute the group of R^2
  E = squares(n,d,S,e,seed)
	#print("squares2", E)

  #S: check for non-squares in squares2
  for rel in E:
    for fac in rel.factors:
      if fac % 2 != 0:
        print("Error in checking for squares:", rel)
        break

  print("Adjoining square roots")
  sys.stdout.flush()
  #compute the square roots and adjoin to the original
  S = adjoin_sqrt(n,S,E)

  M = None
  T = S
  m = len(T[0].factors)
  if list(d) != finald:
    if SHORTENNING_ENABLED:
      print("shortening relations")
      if not SHORTENNING_USE_LLL:
        print("-> using HNF ...")
        M,L = shortening_matrix(n,S) # this thing uses HNF and requires a lot of memory
      else:
        print("-> using LLL ...")
        #L = shortening_matrix_alt(n,S) # does not works, but wants to use LLL in some
        M,L = shortening_matrix_yet_another(n,S) # direct application of LLL

    print("shorten")
    T = shorten(n,d,S,L,seed)
    #T = uniq(sorted(T))
    T = list(set(T))

    if STORE_RELATIONS:
      print(f"Saving relations for field {d}")
      if M == None:
        M = matrix(ZZ,[Ti.factors for Ti in S])
      save_matrix(d, M)

    if COMPUTE_CLGP_FOR_SUBFIELDS:
      print(f"Computing class group structure ...")
      if M == None:
        M = matrix(ZZ,[Ti.factors for Ti in S])
      el_divisors =  M.elementary_divisors()
      print(f'-> class group for d = {d}: {[e for e in el_divisors if e not in [0,1]]}')

    #cheating for filtered out units
    U = units.generators(tuple(d))
    #print('U:', U)
    # add the other units (not -1)
    I = [1]
    for di in d:
      I += [di*i for i in I]
      unitpointers =  [pointerdic[i][0]+pointerdic[i][1]-1 for i in I[1:]]
    for u in U[1:]:
      if evaluatereal(u.element,d) < 0: u.element = -u.element
      assert len(unitpointers) == len(u.powers) - 1
      #assert len(unitpointers) == N - 1
      unitpow = aism*[0]
      unitconj = aism*[0]
      for i in range(N-1):
        if ais[unitpointers[i]].numer.c[0] < 0 or ais[unitpointers[i]].numer.c[1] < 0: 
          if u.powers[i+1] < 0:
            unitpow[unitpointers[i]] = abs(u.powers[i+1])
          else:
            unitconj[unitpointers[i]] = abs(u.powers[i+1])
        else: 
          if u.powers[i+1] < 0:
            unitconj[unitpointers[i]] = abs(u.powers[i+1])
          else:
            unitpow[unitpointers[i]] = abs(u.powers[i+1])
      T += [R(unitpow,unitconj, m*[0])]
  else:
    #print matrix(ZZ,[Ti.factors for Ti in T]).str()
		#print "class group for field", d
    M = matrix(ZZ,[Ti.factors for Ti in T])
    M_HNF = M.hermite_form(include_zero_rows = False)
    el_divisors =  M_HNF.elementary_divisors()
    print(f'Class group: {el_divisors}')
    if STORE_RELATIONS:
      print(f"Saving relations for field {d}")
      save_matrix(d, M_HNF)
	#
	# get_cl_gens() is rather slow
	# TODO: optimize
	#
    #generators = get_cl_gens_optimized(d, T, food, False)
    #print(f'Class group generators: {generators}')
  if False:
    print("powers")
    print(matrix(ZZ,[Ti.powers for Ti in T]).str())
    print("conj")
    print(matrix(ZZ,[Ti.conj for Ti in T]).str())
    print("factors")
    print(matrix(ZZ,[Ti.factors for Ti in T]).str())

  print(f"done with field: {d}\n")
  sys.stdout.flush()

  return tuple(T)

#
#	This function returns the generators of the class group
# It follows the explanation in Biasse-V. Chap. 3A
# ATTENTION: in BV, Chap 3A there is a bug:
# we should construct the generators g using the elems of V^-1, not V!
# Computations can be constly for large relation-set T!
#
@nprofile.profile
def get_cl_gens_optimized(d, T, food, genetors_explicit):
  hnf = matrix(ZZ,[Ti.factors for Ti in T])
  M_HNF = hnf.hermite_form(include_zero_rows = False)
  t = walltime()
  D,U,V=hnf.smith_form()
  el_divisors =  [D[i,i] for i in range(D.ncols())]
  print(f'el_divisors: {el_divisors}')
  Vinv = matrix(ZZ, V.inverse())
  assert len(el_divisors)== V.ncols()
	#print(D) # D is a diagonal with the elem. divisor on the diagonal
	#print(U*hnf*V) # this product gives elementary_divisors and zeros at the bottom
	# Load primes for output field
  generators = []
  if genetors_explicit: # This is VERY slow
    file = convert_letters(d, seed, food=food)
    P_fin=fb.load_primes(tuple(d), file, food=food)
    gens = list(sorted(food.keys()))
    K = NumberField(x^2 - d[0], names=gens[0])
    field_counter = 1
    for gen in gens[1:]:
       L = K.extension(x^2 - d[field_counter], names = gen)
       K = L
       field_counter+=1
    generators = []
    for j in range(V.ncols()):
       if el_divisors[j]==1: continue
       g = K.ideal(1)
       for i in range(len(P_fin)):
           pi = K.ideal([P_fin[i].prime, P_fin[i].elts])
           pi_power = pi^(Vinv[j,i])
           g = g*pi_power
       generators.append(g)
  return generators

def evaluatereal(f, d):
   I = [1]
   for di in d:
     I += [di*i for i in I]
   val = sum([f.numer.c[i]*sqrtapprox[I[i]] for i in range(len(I))])
   return val/f.denom
