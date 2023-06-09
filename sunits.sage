import ring
import field
import field_compact
import relations
import trees
import clgp
from memoized import memoized

@memoized
def sunits(d, d_parent=()):
  K = field.field(d)
  KC = field_compact.field_compact(d)

  class SU:

    def __init__(f,*args, **kwargs):
      if len(args) == 1:
        f.element = args[0]
        return
      raise Exception(f'not known how to initialize sunits.sunits({d})({args})')

    def __repr__(f):
      return 'sunits.sunits(%s)(%s)' % (d,f.element)
  
    def __hash__(f):
      return f.element.__hash__()

    @staticmethod
    def relations_matrix():
      return load(d)[1]

    @staticmethod
    def generators():
      return load(d)[0]
    
    @staticmethod
    def gens():
      return SU.generators()

    @staticmethod
    def factor_base(sage=False):
      return clgp.get_FB(d, d_parent=d_parent, sage=sage)

    @staticmethod
    def preload():
      return SU.generators()
    
    @staticmethod
    def one():
      return KC.one()
    
    @staticmethod
    def from_prime_prod(h):
      '''
      Constructs an S-unit from prime product prod_i FB[i]^h[i]. Assumes that the prime product is principal ideal.
      '''
      FB = SU.factor_base(sage=False)
      assert len(h) == len(FB)
      M = SU.relations_matrix()
      v = M.solve_left(vector(h))
      hh = prod([SU.gens()[i]^v[i] for i in range(len(v))], KC.one())
      return hh

  return SU

def save(d, rels, A, Ac):
  facts = matrix(ZZ,[Ti.factors for Ti in rels])
  relations.save_matrix(d, facts)
  pows = matrix(QQ,[Ti.powers for Ti in rels])
  relations.save_matrix(d, pows, suffix="_powers")
  conj = matrix(QQ,[Ti.conj for Ti in rels])
  relations.save_matrix(d, conj, suffix="_conj")
  relations.save_matrix(d, A, suffix="_base")
  relations.save_matrix(d, Ac, suffix="_base_conj")

@memoized
def load(d):
  KC = field_compact.field_compact(d)
  t = walltime()
  M = relations.load_matrix(d)
  print(f"[relations|{walltime(t)} sec.] ", end="", flush=True)
  t = walltime()
  M_pows = relations.load_matrix(d, suffix="_powers")
  print(f"[powers|{walltime(t)} sec.] ", end="", flush=True)
  t = walltime()
  M_conj = relations.load_matrix(d, suffix="_conj")
  print(f"[conj|{walltime(t)} sec.] ", end="", flush=True)
  t = walltime()
  A = relations.load_matrix(d, suffix="_base", format="list")
  print(f"[base|{walltime(t)} sec.] ", end="", flush=True)
  t = walltime()
  Ac = relations.load_matrix(d, suffix="_base_conj", format="list")
  print(f"[base_conj|{walltime(t)} sec.] ", end="", flush=True)
  #print(f"M = {M}")
  #print(f"A = {A}")
  assert M.nrows() == M_pows.nrows()
  assert M_pows.ncols() == len(A)
  SU = []
  t = walltime()
  for i in range(M_pows.nrows()):
    #su = KC.one()
    #for j in range(M_pows.ncols()):
    #   su *= KC(A[j]) ^ M_pows[i,j] *  KC(Ac[j]) ^ M_conj[i,j]
    su = KC(tuple(A) + tuple(Ac), tuple(M_pows[i]) + tuple(M_conj[i])).trim()
    SU.append(su)
  print(f"[init|{walltime(t)} sec.] ", end="", flush=True)
  return SU,M
