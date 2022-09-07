import sys
import units

nrange = range(3,11)
loops = 100

def gs_smallest(S): # S is a tuple of units
  M = matrix(QQ,[u.approxlog for u in S])
  G,T = M.gram_schmidt()
  return min(sqrt(RR(r*r)) for r in G.rows())

totalsqUgs = {}
totallgUgs = {}
totallgQgs = {}
totallgindex = {}

for n in nrange:
  B = prime_range(2*n*n)
  totalsqUgs[n] = 0
  totallgUgs[n] = 0
  totallgQgs[n] = 0
  totallgindex[n] = 0
  for loop in range(loops):
    d = tuple(sorted(Subsets(B,n).random_element()))
    print('% d',d)
    n = len(d)
    N = 2^n

    U = units.generators_mod_torsion(d)
    Ugs = gs_smallest(U)
    print('% Ugs',n,Ugs)
    lgUgs = log(Ugs)/log(2.0)
    totalsqUgs[n] += Ugs*Ugs
    totallgUgs[n] += lgUgs

    Q = units.mqgenerators_mod_torsion(d)
    Q = sorted(Q,key=lambda u:RR(u.approxlog[0]^2))
    Qgs = gs_smallest(Q)
    print('% Qgs',n,Qgs)
    lgQgs = log(Qgs)/log(2.0)
    totallgQgs[n] += lgQgs

    lgindex = units.mqindex(d).log(2)
    print('% lgindex',n,lgindex)
    totallgindex[n] += lgindex

    sys.stdout.flush()

avgsqUgs = {}
avglgUgs = {}
avglgQgs = {}
avglgindex = {}
for n in nrange:
  avgsqUgs[n] = totalsqUgs[n] / loops
  avglgUgs[n] = totallgUgs[n] / loops
  avglgQgs[n] = totallgQgs[n] / loops
  avglgindex[n] = totallgindex[n] / loops

print('$n$&' + '&'.join('%d' % n for n in nrange) + '\\\\')
print('\\noalign{\\hrule}')
print('average $\\log_2||u^*||$ for $U_L$&' + '&'.join(str('%.3f' % avglgQgs[n]) for n in nrange) + '\\\\')
print('average $\\log_2||u^*||$ for $\\units\\LL$&' + '&'.join(str('%.3f' % avglgUgs[n]) for n in nrange) + '\\\\')
print('% average $||u^*||^2$ for $\\units\\LL$&' + '&'.join(str('%.3f' % avgsqUgs[n]) for n in nrange) + '\\\\')
print('average $\\log_2(\\#(\\units\\LL/U_L))$&' + '&'.join(str('%.3f' % avglgindex[n]) for n in nrange) + '\\\\')
print('% loops',loops)
