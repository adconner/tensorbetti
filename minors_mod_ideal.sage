from itertools import *
from functools import cache
import numpy as np
import scipy.sparse as sp

load('tensorbetti.sage')
load('minors_ideal.sage')

def get_mat(ps):
    mons = sorted(set(m for p in ps for m in p.monomials()),reverse=True)
    print (1)
    return matrix([[p.monomial_coefficient(m) for m in mons] for p in ps])

def generic_minor_samp(M,minorsize):
    F = M.base_ring().base_ring()
    return random_matrix(F,minorsize,M.nrows())*M*random_matrix(F,M.ncols(),minorsize)

def generic_sparse_minor_samp(M,minorsize):
    mons = sorted(set([m for p in M.list() for m in p.monomials()]),reverse=True)
    def to_vec(p):
        return [p.coefficient(m) for m in mons]
    F = M.base_ring().base_ring()
    M = random_matrix(F,minorsize,M.nrows())*M
    js = [0]
    for i in range(M.nrows()):
        j = js[-1]
        print (i,j)
        if j == M.ncols():
            js.extend([M.ncols()]*(M.nrows()-i))
            break
        Mrow = matrix([to_vec(p) for p in M[i,j:].list()])
        r = Mrow.rank()
        Mrow = Mrow.augment(identity_matrix(F,Mrow.nrows()),subdivide=True)
        Mrow.echelonize()
        act = Mrow.subdivision(0,1)
        M = M*block_diagonal_matrix([identity_matrix(F,j),act.T])
        js.append(j+r)
    rtrans = random_matrix(F,M.ncols(),minorsize)
    for j in range(minorsize-1,-1,-1):
        # in order that the lower right j:,j: block has nontrivial determinant, necessary to choose 
        # js[k] so that minorsize-k >= minorsize-j and (js[k] < M.ncols() and js[k] < js[k+1])
        # and M.ncols() - js[k] >= minorsize-j
        # for now, choose js[k] largest satisfying this property (giving most
        # sparse, least general jacobian)
        k = j
        cond = lambda k: js[k] < M.ncols() and js[k] < js[k+1] and \
            M.ncols() - js[k] >= minorsize - j
        while k >= 0 and not cond(k):
            k -= 1
        assert cond(k)
        rtrans[:js[k],j] = 0
    return M*rtrans

# M from simplify_polynomial_matrix
def generic_sparse_minor_samp2(M,minorsize):
    js = [0]
    j = 0
    for r in M:
        while j < len(r) and not r[j].is_zero():
            j += 1
        js.append(j)
    print (js)
    ltrans = random_matrix(F,minorsize,M.nrows())
    rtrans = random_matrix(F,M.ncols(),minorsize)
    for j in range(minorsize-1,-1,-1):
        # k: number of zero entries beginning jth column
        # in order that the lower right j:,j: block has nontrivial determinant, necessary to choose 
        # js[k] so that minorsize-k >= minorsize-j and (js[k] < M.ncols() and js[k] < js[k+1])
        # and M.ncols() - js[k] >= minorsize-j
        # for now, choose js[k] largest satisfying this property (giving most
        # sparse, least general jacobian)
        k = j
        cond = lambda k: js[k] < M.ncols() and js[k] < js[k+1] and \
            M.ncols() - js[k] >= minorsize - j
        while k >= 0 and not cond(k):
            k -= 1
        assert cond(k)
        rtrans[:js[k],j] = 0
        ltrans[:k,k:] = 0
    return ltrans*M*rtrans
    

# looks for row permutation and column GL putting m in approximately upper triangular form
def simplify_polynomial_matrix(M):
    mons = sorted(set([m for p in M.list() for m in p.monomials()]))
    def to_vec(p):
        return [p.coefficient(m) for m in mons]
    F = M.base_ring().base_ring()
    j = 0
    for i in range(M.nrows()):
        print (i,j)
        if j == M.ncols():
            break
        i2 = min(range(i,M.nrows()),key=
                 lambda i2: (matrix([to_vec(p) for p in M[i2,j:].list()]).rank(),i2))
        M.swap_rows(i,i2)
        Mrow = matrix([to_vec(p) for p in M[i,j:].list()])
        r = Mrow.rank()
        Mrow = Mrow.augment(identity_matrix(F,Mrow.nrows()),subdivide=True)
        Mrow.echelonize()
        act = Mrow.subdivision(0,1)
        M = M*block_diagonal_matrix([identity_matrix(F,j),act.T])
        j += r
    return M
        

def reduce_fn_memo(I):
    ltI = ideal([p.lt() for p in I.groebner_basis()])
    maxgbgen = max(p.degree() for p in ltI.gens())
    # assume m = n*x, n not in ltI and x a variable
    @cache
    def reduce1(m):
        if m.degree() <= maxgbgen:
            return I.reduce(m)

        if m not in ltI:
            return m

        # exists since there are no generators of ltI in this degree
        y = next(y for y in m.variables()[::-1] if m.quo_rem(y)[0] in ltI)

        ml = m.quo_rem(y)[0]
        ml = reduce1(ml)

        res = sum([ml.monomial_coefficient(n) * reduce1(n*y) for n in ml.monomials()])
        print ('reduce',m)
        return res
    def reduce(p): # reduction of lin*pnormal mod I
        return p.parent().sum([p.monomial_coefficient(m)*reduce1(m) for m in p.monomials()])
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
        p = sum([(-1)**ji * m[i,j]*minor(cols[:ji]+cols[ji+1:]) for ji,j in enumerate(cols) if not m[i,j].is_zero()])
        p = reduce(p)
        print (cols,len(p.monomials()))
        return p
    return minor(tuple(range(m.ncols())))

def det_red_batch(ms,I,mmaps=None):
    if mmaps is None:
        mmaps = mult_maps(I)
    _,reducemap,mult_maps_noreduce_stacked = mmaps

    n = ms[0].nrows()
    R = I.ring()
    F = R.base_ring()
    char = int(F.characteristic())
    m = np.array([[ [[int(m[i,j].monomial_coefficient(x)) 
                     for m in ms] for x in R.gens()] for j in range(n)] for i in range(n)],
                 dtype=np.int32)

    def multiply(x,p,d):
        intoout, intoin = mult_maps_noreduce_stacked(d)
        rmap = reducemap(d+1)
        prod = (x.astype(np.int32).reshape(x.shape[0],1,x.shape[1]) *\
                  p.reshape(1,p.shape[0],p.shape[1])).reshape(-1,p.shape[1])
        prod %= char
        prod = prod.astype(np.double)
        res = intoin @ prod
        res = rmap @ res
        res += intoout @ prod
        res = res.astype(np.int64)
        if (res > 2**53).any():
            raise ValueError("floating point precision loss")
        res %= char
        return res.astype(np.int32)
        
    @cache
    def minor(cols):
        print ('minor',cols)
        if len(cols) == 0:
            return np.ones((1,len(ms)),dtype=np.int32)
        i = n - len(cols)
        p = sum([(-1)**ji * multiply(m[i,j],minor(cols[:ji]+cols[ji+1:]),len(cols)-1) for ji,j in enumerate(cols) if (m[i,j] != 0).any()])
        p %= char
        print (cols,p.shape[0])
        return p
    return matrix(F,minor(tuple(range(n))))

def mult_maps(I):
    R = I.ring()
    F = R.base_ring()
    ltI = ideal([p.lt() for p in I.groebner_basis()])
    gbmaxd = max([p.degree() for p in I.groebner_basis()])
    xsix = {x:i for i,x in enumerate(R.gens())}
    char = int(F.characteristic())

    @cache
    def mons(d):
        print('mons',d)
        if d == 0:
            return ((R.one(),),())
        monsout = set()
        monsin = set()
        for m in mons(d-1)[0]:
            for x in R.gens():
                n = m*x
                if n in ltI:
                    monsin.add(n)
                else:
                    monsout.add(n)
        return (tuple(sorted(monsout)),tuple(sorted(monsin)))
    @cache
    def monsix(d):
        print('monsix',d)
        monsout,monsin = mons(d)
        return ({m:i for i,m in enumerate(monsout)}, {m:i for i,m in enumerate(monsin)})
    @cache
    def reducemap(d):
        print('reducemap',d)
        monsout,monsin = mons(d)
        if d <= gbmaxd:
            mat = []
            for m in monsin:
                p = I.reduce(m)
                mat.append([p.monomial_coefficient(n) for n in monsout])
            # do this to get shape correct when len(monsin) == 0
            return np.array(matrix(len(monsin),len(monsout),mat).T,dtype=np.int16)

        # each monomial, divide by y staying in ltI, reduce there, multiply by y, then reduce here (using only lower monomials)

        seen = set()
        divmaps = [None for yi in range(R.ngens())]
        for yi in range(R.ngens()-1,-1,-1):
            y = R.gen(yi)
            cur = {}
            for mli,ml in enumerate(mons(d-1)[1]):
                n = ml*y
                if n in monsix(d)[1] and n not in seen:
                    seen.add(n)
                    cur[(mli,monsix(d)[1][n])] = 1
            divmaps[yi] = matrix(F,len(mons(d-1)[1]),len(mons(d)[1]),cur)

        # mapiot = matrix(F,len(monsin),len(monsout))
        # mapiit = matrix(F,len(monsin),len(monsin),sparse=True)
        mapiot = np.zeros((len(monsin),len(monsout)),dtype=np.int32)
        mapiit = sp.lil_array((len(monsin),len(monsin)),dtype=np.int16)
        rprevt = reducemap(d-1).T
        for yi,(divmap,(intoout,intoin)) in enumerate(zip(divmaps,mult_maps_noreduce(d-1))):
            if len(divmap.dict()) == 0:
                continue
            print(yi,end=" ",flush=True)
            intoin = np.array(sorted(intoin.dict().keys())).T
            intoout = np.array(sorted(intoout.dict().keys())).T
            divmap = np.array(sorted(divmap.dict().keys())).T

            if len(intoin) > 0:
                mapiit[divmap[1].reshape(-1,1),intoin[0]] = rprevt[divmap[0].reshape(-1,1),intoin[1]]
            mapiot[divmap[1].reshape(-1,1),intoout[0]] = rprevt[divmap[0].reshape(-1,1),intoout[1]]
        print("solving %d %d %d.. "% (mapiit.nnz, len(monsin),len(monsout)),end="",flush=True)
        mapiit = mapiit.tocsc()
        for i,(a,b) in enumerate(zip(mapiit.indptr,mapiit.indptr[1:])):
            J = mapiit.indices[a:b]
            E = mapiit.data[a:b]
            mapiot[J] += mapiot[i]*E.reshape(-1,1)
            mapiot[J] %= char
        print('finished')
        return mapiot.T.astype(np.int16)

    @cache 
    def mult_maps_noreduce(d): # d -> d+1
        print('mult_maps_noreduce',d)
        mlo, _ = mons(d)
        monsout, monsin = mons(d+1)
        monsoutix, monsinix = monsix(d+1)
        maps = []
        for x in R.gens():
            intoin = {}
            intoout = {}
            # map = matrix(F,len(monsout),len(mlo))
            for mi,m in enumerate(mlo):
                n = m*x
                if n in monsinix:
                    intoin[(monsinix[n],mi)] = 1
                else:
                    intoout[(monsoutix[n],mi)] = 1
            intoin = matrix(F,len(monsin),len(mlo),intoin)
            intoout = matrix(F,len(monsout),len(mlo),intoout)
            maps.append((intoout,intoin))
        return maps
    @cache
    def mult_maps_noreduce_stacked(d):
        return (block_matrix([[a for a,_ in mult_maps_noreduce(d)]]),
                block_matrix([[b for _,b in mult_maps_noreduce(d)]]))

    @cache
    def mult_maps_noreduce_stacked_numpy(d):
        intoout, intoin = mult_maps_noreduce_stacked(d)
        dat = list(intoin.dict().keys())
        I = np.array([i for i,_ in dat])
        J = np.array([j for _,j in dat])
        E = np.ones(len(dat),dtype=np.double)
        intoin = sp.coo_array((E,(I,J)),shape=intoin.dimensions(),dtype=np.int8)
        intoin = intoin.tocsr()

        dat = list(intoout.dict().keys())
        I = np.array([i for i,_ in dat])
        J = np.array([j for _,j in dat])
        E = np.ones(len(dat),dtype=np.double)
        intoout = sp.coo_array((E,(I,J)),shape=intoout.dimensions(),dtype=np.int8)
        intoout = intoout.tocsr()
        return intoout,intoin

    return mons,reducemap,mult_maps_noreduce_stacked_numpy


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


# vim: ft=python
