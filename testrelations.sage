import os
import sys
import nprofile
import units
import goodprime
import goodprimecheat
import mult
import div
import subsetprod
import powerprod
import sqrt
import field
import relations
import timeit
import fb
import trees
from sage.combinat.subset import SubsetsSorted
import argparse

def parse_files(d, dick, food = None, gens = None):
  if food == None:
    food = {"a": 5, "b": 13, "c": 17, "d": 29, "e": 37, "f": 41, "g": 53, "h": 61}
  if gens == None:
    gens = list(food.keys())
  n = len(d)
  file = gens[0] + '_'
  aDone = False
  for i in SubsetsSorted(gens[1:n])[1:]:
    #print(i, file)
    count = 0
    counting = False
    counter = 0
    currentfile = file + ''.join(i)
    f = open("trees/" + currentfile + ".christine","r")
    #if d[0]>0:
    #  f = open("trees/" + currentfile + ".christine", "r")
    #else:
    #  f = open("trees/imaginary/" + currentfile + ".christine", "r")
    #print('openning ', "trees/" + currentfile + ".christine")
    for line in f:
      if line[0] == 'F':
        line = line[:-1]
        spl = line.split(' ')
        lst = spl[2]
        field = []
        for elt in lst.split('_'):
          #print elt
          #field += [prod(int(food[k]^elt.count(k)) for k,_ in food.iteritems())]
          field += [prod(int(food[k]^elt.count(k)) for k,_ in food.items())]
        if len(field) == 1:
          #if field[0] == 5 and aDone == True:
          #if field[0] == food.values()[0] and aDone == True: 
          if field[0] == d[0] and aDone == True: 
            counting = False
            continue 
          #if field[0] == 5 and aDone == False: aDone = True 
          #if field[0] == food.values()[0] and aDone == False: aDone = True
          if field[0] == d[0] and aDone == False: aDone = True
          counting = True
          counter = field[0]
        else:
          counting = False
      else:
        if counting == True: 
          dick[counter] += 1/3
    #print 'dick:', dick
    f.close()
  #for k, v in dick.iteritems():
  for k, v in dick.items():
    dick[k] = ZZ(dick[k] + 1)
    #print k, v, dick
  return dick

def make_dictionary(d):
   dick = {}
   for i in Subsets(d)[1:]:
      dick[prod(i)] = 0
   return dick

#returns (# quad relations, dictionary of dick[gennorm] = (#rels Q(gen,), where in vec start), list of pointers, list of gens, dictionary of conjugation vectors)
def get_info(d, pars):
   I = [1]
   for di in d:
     I += [di*j for j in I]
   aism = 0
   point = [] #S: offsets of S-units of Q(sqrt(di)) in ais
   for di in I[1:]:
      point += [aism]
      #aism += pars[di]
      if di > 0:
        #S: quadratic real case, the number of S-unit group generators is |S| + 1
        pars[di] = (point[-1], pars[di])
      else:
        #S: quadratic imaginary case, the number of S-unit group generators is |S|
        K = NumberField(x^2 - di, names="z") 
        if K.unit_group(proof=False).order() == 2: #S: TODO, replace this by explicit conditions
          pars[di] = (point[-1], pars[di] - 1)
        else:
          print("[info] Quadratic subfield with di = {} has non-trivial units".format(di))
          pars[di] = (point[-1], pars[di]) #S: there are non-trivial units
      aism += pars[di][1]
			#print('pars interm:', pars)
   sigmas = {}
   for di in I[1:]:
     sigmas[di] = get_sigma(I, di, aism, pars)
   return aism, pars, point, I[1:], sigmas

def get_sigma(I, di, aism, pars):
  vec = aism*[1]
  for i in I:
    if i % di == 0:
      vec[pars[i][0]:pars[i][0] + pars[i][1]] = pars[i][1]*[-1]
  return vec

#d = (5,13,17,29,37,41,53,61)
#d = (5,13,17,29,37,41,53)
#d = (5,13,17,29,37,41)
#d = (5,13,17,29,37)
#d = (5, 13, 17, 29)
#d = (5,13,17)
#d = (5,13)
#d = (5,377)
#d = (5,37)
#d = (5,)
#d = (5,29)
#d = (5,17,53)

# Note: don't use labels d,e,f simultaneously.
#food = {"a": 5, "b": 13, "c": 17}
#food = {"a": -3, "b": -11, "c": -71} #|Cl_K| = 56, norm_bound = 37 (GM)
#food = {'a': -7, 'b': -11, 'c': -43} # |Cl_K| = 54, norm_bound = 25 (GM)
#food = {'a': -7, 'b': -11, 'c': -67} # |Cl_K| = 66
#food = {'a': -7, 'b': -19, 'c': -43}
#food = {'a': -7, 'b': -23, 'c': -31} # |Cl_K| = 207
#food = {'a': -7, 'b': -31, 'c': -47} # |Cl_K| = 510
#food = {'a':17, 'b':29, 'c': 37, 'd':41}
#food = {'a':-7, 'b':-11, 'c':-19, 'd':-23}
#food = {'a':-11, 'b':-19, 'c':-23, 'd':-31, 'e':-47} 
#food = {'a': -11, 'b': -19, 'c': -23, 'd': -31, 'e': -43, 'g': -47}
#food = {'a': -11, 'b': -19, 'c': -23, 'd': -31, 'e': -43} #Cl_K: C2 x C2 x C4 x C4 x C12 x C48 x C96 x C3360 x C3360 x C443520
#food = {'a':-11, 'b':-19, 'c':-31}
#food = {'a': -11, 'b': -19, 'c': -23}
#food ={'a':-11, 'c':-31}
#food = {"c": 17, "f": 41, "g": 53} # |Cl_K| = 36
#food = {"a": -7, "b": -23}  # |Cl_K| = 3
#food = {'a': 5, 'b': 79, 'c': 89} # Sage: |Cl_K| = 96, Result: 192
#food = {'a': 19, 'b': 23}
#food = {'a':-19, 'b':-71}
#food = {'a': 5, 'b': 79, 'c': 89, 'd': 97, 'e': 101} # result: C2 x C2 x C2 x C2 x C2 x C2 x C2 x C2 x C4 x C4 x C4 x C24 x C48 x C48 x C96 x C96
#food = {'a': -11, 'b': -19, 'c': -31}
#food = {'a': -7, 'b': -11, 'c': -19}
#d = tuple(sorted(food.values()))

#food = {'a':-11, 'b':-19, 'c':-23, 'd':-31} #Cl_K: C3 x C3 x C12 x C2520
#food = {'a': 7*11, 'b': 11*43}
#food = {'a': 4759, 'b': 8761}
#food = {'a': 7057, 'b': 8101, 'c': 8761}
#food = {'a': 229, 'b': 257, 'c': 359, 'd': 401}
#food = {'a': 5, 'b': 229, 'c': 257}

#food = {'a': -7, 'b': -11, 'c': -19, 'd': -23} # |Cl_K| = 23040, norm_bound = 361 (GM)
#food = {"a": 5, "b": 13, "c": 17}
#food = {"a": 5, "b": 13, "c": 17, "d": 29}
#food = {"a": 5, "b": 13, "c": 17, "d": 29, "e": 37} #norm_bound = 2209 (GM)
#food = {"a": 5, "b": 13, "c": 17, "d": 29, "e": 37, "g": 41}
#food = {"a": 5, "b": 13, "c": 17, "d": 29,  "e": 37, "g": 41, "h": 53}
#food = {"a": -2777106107, "b": -2777106487, "c": -2777106991}
#food = {'a': 1297, 'b': 7057, 'c': 8761} # h_1297 = 11, h_4759 = 13, h_7057 = 21

#food = {"a": -3, "b": -7, "c": -11, "d": -19}
#food = {"a": -3, "b": -7, "c": -11, "d": -19, "e": -23}
food = trees.get_food() # load food from trees

seed = 1
def testrelations(food, seed):
	d = tuple(sorted(food.values(), key=lambda di: abs(di)))
	gens = list(sorted(food.keys()))
	#print('d:', d, 'gens:', gens)

	K = field.field(d)
	print("Regulator: ", units.approxregulator(d))
	dick = make_dictionary(d)
	print("dick:", dick)
	pars = parse_files(tuple(d), dick, food = food)
	print('pars:', pars)
	aism, pars, points, I, sigmas = get_info(d, pars)
	print('pars:', pars)
	relations.set_vars(aism, sigmas, pars, d)
	t = walltime()
	file = "_".join(gens)
	res = relations.relations_internal(tuple(d), file, food = food, seed = seed)
	print('len(relations)', len(res))
	#cProfile.run('relations.relations_internal(tuple(d), file, food = food)')
	print("total computation: ", walltime(t))
	# f = open("relations/" + str(d)+'_relations','w')
	# for i in range(len(res)):
	# 	f.write(str(list(res[i].factors))+ "\n" )
	# f.close()

	return 1

parser = argparse.ArgumentParser(description='Generation of trees describing splitting of primes over subfields of multiquadratic field.')
parser.add_argument("food", type=str, nargs="*", help="dictionary of the form {'a': d_1, ..., 'e': d_n} containing information of multiquadratic field and corresponding labels for generators")
parser.add_argument("--store", action=argparse.BooleanOptionalAction, default=True, help="store (or not) relations matrices")

parser.add_argument('--shorten', '-s', choices=['hnf', 'lll', 'bkz'], default="lll", help="algorithm for shortening of intermediate relation matrices.")
parser.add_argument('--shorten-final', '-sf', choices=['hnf', 'lll', 'bkz'], default="lll", help="algorithm for shortening of final relation matrix.")
parser.add_argument('--bkz-bsize', type=Integer, help="Block size for BKZ (default = log(d_1 * ... * d_n) / 2)")

args = parser.parse_args()

relations.set_shortening_alg(args.shorten)
relations.set_final_shortening_alg(args.shorten_final)

if args.food != []:
  food = eval(preparse(" ".join(args.food)))

trees_food = trees.get_food()

if trees_food != None:
  assert food == trees_food, f"Run trees generation first! {food} != {trees_food}."

if args.shorten == "bkz" or args.shorten_final == 'bkz':
  if args.bkz_bsize != None:
    relations.set_bkz_block_size(args.bkz_bsize)
  else:
    dd = prod(food.values()).abs()
    #bs = floor(log(dd)/2)
    bs = floor((log(dd))^(1/2) * log(log(dd)) * len(food.values()))
    # for small dd we use default value of block size
    if bs < 40:
      bs = 40
    relations.set_bkz_block_size(bs)

trees.clear_dir("relations")
testrelations(food, seed)
#nprofile.output([goodprime,units,relations])
