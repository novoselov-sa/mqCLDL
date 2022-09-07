import sys
from sage.doctest.util import Timer

class profile(object):
  def __init__(self, func):
    self.func = func
    self.time = {}
    # self.time[n] is (c,t)
    # where c is number of calls with first arg n
    # and t is time spent in those calls
  def __call__(self, *args, **kwargs):
    #n = args[0]
    n = str(args[0])
    timer = Timer()
    timer.start()
    try:
      result = self.func(*args, **kwargs)
    finally:
      s = timer.stop().cputime
      c,t = 0,0
      if n in self.time:
        c,t = self.time[n]
      self.time[n] = c+1,t+s
      c,t = 0,0
      if None in self.time:
        c,t = self.time[None]
      self.time[None] = c+1,t+s
    return result

def output(modules):
  result = ''
  result2 = ''

  for m in modules:
    for x in m.profiling:
      try:
        f = m.__dict__[x]
        for n in sorted(f.time):
          y = '%13.5f = %10s * %8.5f: %s.%s' % (f.time[n][1],f.time[n][0],f.time[n][1]/f.time[n][0],m.__name__,x)
          z = '%s(%s)' % (y,n)
          result += z + '\n'
          if n == None: result2 += y + '\n'
      except:
        pass

  sys.stdout.write(result + result2)
