import collections

class memoized(object):
  def __init__(self,func):
    self.func = func
    self.cache = {}
    self.__name__ = 'memoized:' + func.__name__
  def __call__(self, *args, **kwargs):
#    if not isinstance(args, collections.Hashable):
    if not isinstance(args, collections.abc.Hashable):
      return self.func(*args, **kwargs)
    #if not args in self.cache:
      #self.cache[args] = self.func(*args, **kwargs)
    #return self.cache[args]

    kw = [el for el in kwargs.items() if isinstance(el[1], collections.abc.Hashable)]

    k = frozenset([("args", args)] + kw)
    if not k in self.cache:
      self.cache[k] = self.func(*args, **kwargs)
    return self.cache[k]
