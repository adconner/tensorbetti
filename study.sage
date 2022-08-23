from itertools import *
from functools import cache
from IPython import embed

load('tensorbetti.sage')
load('minors_ideal.sage')

a = 4
n = 4
r = 7

# F = GF(5)
# F = GF(101)
# F = GF(1031)
F = GF(32003)
# F = GF(65537)

h = Tinv(random_tensor(F,n,r))
I = h.ideal_to(2)
jac = matrix([[p.derivative(x) for x in I.ring().gens()] for p in I.gens()[::-1]])

minorsize = 12
jacm = random_matrix(F,minorsize,jac.nrows())*jac*random_matrix(F,jac.ncols(),minorsize)


def reduce_fn_memo(I):
    ltI = ideal([p.lt() for p in I.groebner_basis()])
    maxgbgen = max(p.degree() for p in ltI.gens())
    # assume m = n*x, n not in ltI and x a variable
    @cache
    def reduce1(m):
        print ('reduce',m)
        if m.degree() <= maxgbgen:
            return I.reduce(m)

        if m not in ltI:
            return m

        # exists by assumption on m
        x = next(x for x in m.variables() if m.quo_rem(x)[0] not in ltI)
        # exists since there are no generators of ltI in this degree
        y = next(y for y in m.variables() if m.quo_rem(y)[0] in ltI)

        ml = m.quo_rem(x*y)[0]
        # ml not in ltI
        # ml*y not in ltI
        # ml*x in ltI
        # ml*x*y == m in ltI

        # call is smaller because
        mlx = reduce1(ml*x)

        # all calls smaller because
        # ml*x > mlx
        # m = ml*x*y > mlx*y
        res = sum([mlx.monomial_coefficient(n) * reduce1(n*y) for n in mlx.monomials()])
        assert all(m not in ltI for m in res.monomials())
        return res
    def reduce(p): # reduction of lin*pnormal mod I
        print ('reduce_poly')
        return sum([p.monomial_coefficient(m)*reduce1(m) for m in p.monomials()])
    return reduce

def det_red(m,I,reduce=None):
    if reduce is None:
        reduce = reduce_fn_memo(I)

    @cache
    def minor(cols):
        print ('minor',cols)
        nonlocal I
        if len(cols) == 0:
            return I.ring().one()
        i = m.nrows() - len(cols)
        p = sum([(-1)**ji * m[i,j]*minor(cols[:ji]+cols[ji+1:]) for ji,j in enumerate(cols)])
        p = reduce(p)
        print (cols,len(p.monomials()))
        return p
    return minor(tuple(range(m.ncols())))


def mult_maps(I,dupto):
    allmons = { 0 : [I.ring().one()] }
    allmonsix = { 0: { I.ring().one() : 0 } }
    allmats = { }
            
    for d in range(1,dupto+1):
        print(d)
        mons = []
        monsix = {}
        def get(m):
            if m not in monsix:
                monsix[m] = len(mons)
                mons.append(m)
            return monsix[m]
        mats = []
        for x in I.ring().gens():
            print(x)
            mat = {}
            for mi,m in enumerate(allmons[d-1]):
                mr = m*x
                if mr in monsix:
                    mat[(mi,monsix[mr])] = 1
                else:
                    mrred = reduce(mr)
                    if mr == mrred:
                        mat[(mi,get(mr))] = 1
                    else:
                        for n,c in mrred.dict().items():
                            mat[(mi,get(I.ring().monomial(*n)))] = c
            mats.append(mat)
        mats = [matrix(I.base_ring(),len(allmons[d-1]),len(mons),mat) for mat in mats]
        allmats[d] = mats
        allmons[d] = mons
        allmonsix[d] = monsix
    return allmons,allmonsix,allmats

def upper_tri_assume_all_generic(m,I):
    F = I.ring().base_ring()
    S = PolynomialRing(F,I.ring().gens()+('d',))
    dv = S.gens()[-1]
    I = I.change_ring(S)
    m = m.change_ring(S)
    d = 1
    for i in range(m.nrows()):
        f = m[i,i]
        m = m.subs({dv : dv * f})
        I = I.elimination_ideal(dv) + [dv*d*f-1]
        dold = d
        d *= f

        m[i] *= dv*dold
        m[i,i] = 1
        for j in range(i+1,m.nrows()):
            print(i,j)
            m.add_multiple_of_row(j,i,-m[j,i])
        I.groebner_basis()
        print ('found gb')
        m = m.apply_map(lambda e: I.reduce(e))
        print ('simplified')

    return (I,m)

def upper_tri_assume_all_generic_mult(m,I):
    F = I.ring().base_ring()
    S = PolynomialRing(F,I.ring().gens()+tuple('d%d' % i for i in range(minorsize)))
    dvs = S.gens()[-minorsize:]
    I = I.change_ring(S)
    m = m.change_ring(S)
    for i,d in enumerate(dvs):
        I = I + (d*m[i,i] - 1)
        m[i] *= d
        m[i,i] = 1
        for j in range(i+1,m.nrows()):
            print(i,j)
            m.add_multiple_of_row(j,i,-m[j,i])
        I.groebner_basis()
        print ('found gb')
        m = m.apply_map(lambda e: I.reduce(e))
        print ('simplified')

    return (I,m)

# dat = upper_tri_assume_all_generic(jacm,I)

# vim: ft=python
