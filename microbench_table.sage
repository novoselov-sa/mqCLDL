from sage.doctest.util import Timer
timer = Timer()

import goodprime
import mult
import norm
import div
import sqrt

nrange = range(3,11)

x = {}

for n in nrange:
  N = 2^n
  d = (2,3,5,7,11,13,17,19,23,29)[:n]

  bits = 1000
  loops = 101

  timings = []
  for loop in range(loops):
    f = tuple(randrange(1-2^bits,2^bits) for j in range(N))
    g = tuple(randrange(1-2^bits,2^bits) for j in range(N))
    timer.start()
    h = mult.mult(n,N,d,f,g)
    timings += [timer.stop().cputime]
  x['mult',n] = median(timings)

  timings = []
  for loop in range(loops):
    f = tuple(randrange(1-2^bits,2^bits) for j in range(N))
    timer.start()
    h = mult.square(n,N,d,f)
    timings += [timer.stop().cputime]
  x['square',n] = median(timings)
  
  timings = []
  for loop in range(loops):
    f = tuple(randrange(1-2^bits,2^bits) for j in range(N))
    timer.start()
    h = norm.normonestep(n,N,d,f,n-1)
    timings += [timer.stop().cputime]
  x['relnorm',n] = median(timings)
  
  timings = []
  for loop in range(loops):
    f = tuple(randrange(1-2^bits,2^bits) for j in range(N))
    timer.start()
    h = norm.absnorm(n,N,d,f)
    timings += [timer.stop().cputime]
  x['absnorm',n] = median(timings)
  
  timings = []
  for loop in range(loops):
    f = tuple(randrange(1-2^bits,2^bits) for j in range(N))
    g = tuple(randrange(1-2^bits,2^bits) for j in range(N))
    h = mult.mult(n,N,d,f,g)
    timer.start()
    f2 = div.div(n,N,d,h,g)
    timings += [timer.stop().cputime]
    assert f2 == f
  x['div',n] = median(timings)

  timings = []
  for loop in range(loops):
    f = tuple(randrange(1-2^bits,2^bits) for j in range(N))
    h = mult.square(n,N,d,f)
    timer.start()
    f2,S = sqrt.squareroot(n,d,n,h)
    timings += [timer.stop().cputime]
    assert S == set()
    assert f2 == f or f2 == tuple(-fj for fj in f)
  x['sqrt',n] = median(timings)

ops = ['mult','square','relnorm','absnorm','div','sqrt']

print '\\begin{tabular}{r|r|%s}' % ('|r'*len(ops))
print '$n$&$2^n$&%s\\\\' % '&'.join(ops)
print '\\noalign{\\hrule}'
for n in nrange:
  print '%d&%d&%s\\\\' % (n,2^n,'&'.join('%.4f' % x[op,n] for op in ops))
print '\\end{tabular}'
