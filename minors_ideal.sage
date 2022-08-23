from IPython import embed

def minors_ideal(m,r,I=None):
    from collections import Counter

    S = m.base_ring()
    if I is None: I = S.ideal()

    R = PolynomialRing(S.base_ring(),S.gens()+('d',))
    dv = R.gens()[-1]

    def rec(m,I,d,r):
        print (r)
        if m.is_zero():
            return I.elimination_ideal(dv).change_ring(S)
        if r == 1:
            return (I + m.coefficients()).elimination_ideal(dv).change_ring(S)

        # select the most common nonzero element of 
        f = Counter(m.dict().values()).most_common(1)[0][0]

        # Case 1: pass to the closed set where f == 0
        Icur = I + f
        mcur,lo = elimination_by_units(m,Icur)
        I1 = rec(mcur[lo:,lo:],Icur,d,r-lo) if lo < r else S.ideal(1)

        # Case 2: pass to the open set where f not identically zero, e.g., localize at f
        mcur = m.subs({dv : dv*f})
        Icur = I.elimination_ideal(dv) + [dv*d*f-1]
        mcur,lo = elimination_by_units(mcur,Icur)
        I2 = rec(mcur[lo:,lo:],Icur,d*f,r-lo) if lo < r else S.ideal(1)

        # Combine the result of both cases: $I_1\cap I_2$ corresponds to the
        # union of the corresponding parameter sets
        return I1.intersection(I2)

    return rec(m.change_ring(R), I.change_ring(R), R.one(), r)

# M : matrix with base ring a polynomial ring over QQ or a quotient of such by an ideal
# nminors: a parameter determining when to revert to the base case in the
#    recursive algorithm. Output should be independent of its value, but
#    performance may be tuned
#
# returns the radical of the ideal of r by r minors of m
def minors_ideal1(m,r,I=None):
    from collections import Counter

    S = m.base_ring()
    if I is None: I = S.ideal()

    R = PolynomialRing(S.base_ring(),S.gens()+('d',))
    dv = R.gens()[-1]

    def rec(m,I,d,r):
        if m.is_zero():
            yield I
            return
        if r == 1:
            yield I + m.coefficients()
            return

        f = Counter(m.dict().values()).most_common(1)[0][0]

        # pass to the open set where f not identically zero, e.g., localize at f
        mcur = m.subs({dv : dv*f})
        Icur = I.elimination_ideal(dv) + [dv*d*f-1]
        mcur,lo = elimination_by_units(mcur,Icur)
        if lo < r:
            for Iout in rec(mcur[lo:,lo:],Icur,d*f,r-lo):
                yield Iout

        # pass to the closed set where f == 0
        Icur = I + f
        mcur,lo = elimination_by_units(m,Icur)
        if lo < r:
            for Iout in rec(mcur[lo:,lo:],Icur,d,r-lo):
                yield Iout

    return S.ideal(1).intersection(*(I.elimination_ideal(dv).change_ring(S)
            for I in rec(m.change_ring(R), I.change_ring(R), R.one(), r)))

# M : matrix with base ring a polynomial ring over QQ or a quotient of such by an ideal
# nminors: a parameter determining when to revert to the base case in the
#    recursive algorithm. Output should be independent of its value, but
#    performance may be tuned
#
# returns the radical of the ideal of r by r minors of m
def minors_ideal2(m,r,I=None):
    from collections import Counter

    S = m.base_ring()
    if I is None: I = S.ideal()

    def getR(m,I,needed_ds):
        if S.ngens() + needed_ds > m.base_ring().ngens():
            numds = next(fibonacci(k) for k in range(2,1000) if fibonacci(k) >= needed_ds)
            # print (numds)
            R = PolynomialRing(S.base_ring(),S.gens() +
                tuple('d%d' % i for i in range(numds)))
            return m.change_ring(R),I.change_ring(R)
        return m,I

    def rec(m,I,di,r):
        if m.is_zero():
            yield I.elimination_ideal(I.ring().gens()[S.ngens():]).change_ring(S)
            return
        if r <= 2:
            I = I + m.minors(r)
            yield I.elimination_ideal(I.ring().gens()[S.ngens():]).change_ring(S)
            return

        f = Counter(m.dict().values()).most_common(1)[0][0]
        m,I = getR(m,I,di+1)

        # pass to the open set where f not identically zero, e.g., localize at f
        Icur = I + [I.ring().gen(S.ngens()+di)*f - 1]
        mcur,lo = elimination_by_units(m,Icur)
        # print ('one',m.dimensions(),lo)
        if lo < r:
            for Iout in rec(mcur[lo:,lo:],Icur,di+1,r-lo):
                yield Iout

        # pass to the closed set where f == 0
        Icur = I + f
        mcur,lo = elimination_by_units(m,Icur)
        # print ('two',m.dimensions(),lo)
        # print([t for ee in mnext.list() for t in ee.monomials()])
        if lo < r:
            for Iout in rec(mcur[lo:,lo:],Icur,di,r-lo):
                yield Iout

    return S.ideal(1).intersection(*rec(m,I,0,r))

# m : matrix
# 
# row reduces m as far as possible pivoting using unit entries. Column swaps are
# also performed. m is modified in place and the number of pivots r is returned. 
# After this operation, m has the form
# [ T A ]
# [ 0 B ]
# where T is r by r and  upper triangular with units along the diagonal, B 
# contains no units
def elimination_by_units(m,I):
    m = m.apply_map(lambda e: I.reduce(e))

    from sage.libs.singular.function import singular_function
    lift = singular_function('lift')

    # computes the inverse mod $I$ if it exists, otherwise None
    def try_inverse(e):
        try:
            return lift(I+e,1).list()[-1]
        except RuntimeError:
            return None

    r = 0
    while True:
        if r == min(m.nrows(),m.ncols()):
            break

        einv = None
        for (i,j),e in m[r:,r:].dict().items():
            einv = try_inverse(e)
            if einv is not None:
                break
        if einv is None:
            break
        i += r
        j += r

        m.swap_rows(r,i)
        m.swap_columns(r,j)
        m[r,:] *= einv
        for k in m.nonzero_positions_in_row(r):
            m[r,k] = I.reduce(m[r,k])
        assert m[r,r] == 1
        for i in m.column(r)[r+1:].nonzero_positions():
            i += r+1
            m.add_multiple_of_row(i,r,-m[i,r])
            for k in m.nonzero_positions_in_row(i):
                m[i,k] = I.reduce(m[i,k])
        r += 1

    return m,r
  

# vim: ft=python

