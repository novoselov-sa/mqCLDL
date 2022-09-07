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
from sage.combinat.subset import SubsetsSorted


def parse_files(d, dick, food = None, gens = None):
  if food == None:
    food = {"a": 5, "b": 13, "c": 17, "d": 29, "e": 37, "f": 41, "g": 53, "h": 61}
  if gens == None:
    gens = food.keys()
  n = len(d)
  file = gens[0] + '_'
  aDone = False
  for i in SubsetsSorted(gens[1:n])[1:]:
    #print i, file
    count = 0
    counting = False
    counter = 0
    currentfile = file + ''.join(i)
    f = open("trees/" + currentfile + ".christine","r")
    for line in f:
      if line[0] == 'F':
        line = line[:-1]
        spl = line.split(' ')
        lst = spl[2]
        field = []
        for elt in lst.split('_'):
          field += [prod(int(food[k]^elt.count(k)) for k,_ in food.iteritems())]
        if len(field) == 1:
          if field[0] == d[0] and aDone == True: 
            counting = False
            continue
          if field[0] == d[0] and aDone == False: aDone = True
          counting = True
          counter = field[0]
        else:
          counting = False
      else:
        if counting == True: 
          dick[counter] += 1/3
    f.close()
  for k, v in dick.iteritems():
    dick[k] = ZZ(dick[k] + 1)
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
        if K.unit_group().order() == 2: #S: TODO, replace this by explicit conditions
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


def testrelations(food, seed):
   d = tuple(sorted(food.values(), key=lambda di: abs(di)))
   gens = list(sorted(food.keys()))

   K = field.field(d)
	 #print('d:', d)
   dick = make_dictionary(d)
	 #print "dick:", dick
   pars = parse_files(tuple(d), dick, food = food)
	 #print 'pars:', pars
   aism, pars, points, I, sigmas = get_info(d, pars)
	 #print 'pars:', pars
   relations.set_vars(aism, sigmas, pars, d)
   t = walltime()
   file = "_".join(gens)
   print('file:', file)
   res = relations.relations_internal(tuple(d), file, seed, food = food)
	 #print('len(relations)', len(res))
	 #print "total computation: ", walltime(t)
	 #f = open('relations/'+str(d)+'_relations','w')
		 #for i in range(len(res)):
		 # f.write(str(list(res[i].factors))+ "\n" )
		 #f.close()
   return res

#food = {'a': -19, 'c': -31, 'b': -23, 'd': -43}
#seed = 0
#testrelations(food,seed)


