# Computation of square roots for S-smooth ideal I using saturation.
# Requires precomputed class group.

import field
import field_compact

def find_squares(d, T, r):
    '''
    Recognizing squares among products of elements from a set T.
    '''
    low,high = 0,r + 64
    # The procedure works with symbols_log if T doesn't contain units.
    t0_s = T[0].symbols2_log(low,high) 
    #t0_s = T[0].symbols_log(low,high)
    A = matrix(GF(2), len(T), len(t0_s))
    A[0] = t0_s
    for i in range(1,len(T)):
        A[i] = T[i].symbols2_log(low,high)
        #A[i] = T[i].symbols_log(low,high)
    V = A.left_kernel().basis_matrix().rows()
    K = field.field(d)
    KC = field_compact.field_compact(d)
    squares = []
    for v in V:
        sq = KC(K.one())
        for j in range(len(v)):
            sq *= T[j] ^ ZZ(v[j])
        squares.append(sq)
    return squares

def find_square_exp(h, T, r = 0):
    '''
    Given a list T of field elements in and element of the field h the method finds binary vector e such that h*T[i]^e[i] is a square.

    Important: this method doesn't support rational exponents in compact representation of field elements.
    '''
    low,high = 0,r + 64
    assert len(T) != 0
    v = vector(GF(2), h.symbols2_log(low, high))
    
    A = matrix(GF(2), len(T), len(v))
    for i in range(len(T)):
        A[i] = T[i].symbols2_log(low,high)
    e = A.solve_left(v)
    return e

def find_square_rat_exp(h, T, r = 0):
    '''
    Given a list T of field elements in and element of the field h the method finds binary vector e such that h*T[i]^e[i] is a square.

    This method supports rational exponents in compact representation of field elements.
    '''
    low,high = 0,r + 64
    assert len(T) != 0
    v = vector(GF(2), h.altsymbols_log(low, high))
    
    A = matrix(GF(2), len(T), len(v))
    for i in range(len(T)):
        A[i] = T[i].altsymbols_log(low,high)
    e = A.solve_left(v)
    return e
