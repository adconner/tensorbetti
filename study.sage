from itertools import *
from functools import cache
from IPython import embed

load_attach_mode(load_debug=True)

load('tensorbetti.sage')
load('minors_ideal.sage')

def get_mat(ps):
    mons = sorted(set(m for p in ps for m in p.monomials()),reverse=True)
    print (1)
    return matrix([[p.monomial_coefficient(m) for m in mons] for p in ps])

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

def det_red_batch(ms,I,mult_map=None):
    if mult_map is None:
        mult_map = mult_maps(I)[0]
    n = ms[0].nrows()
    R = I.ring()
    F = R.base_ring()
    m = [[ matrix(F,[[m[i,j].monomial_coefficient(x) for x in R.gens()] for m in ms]).T 
        for j in range(n)] for i in range(n)]

    def multiply(x,p,d):
        assert x.ncols() == p.ncols()
        return mult_map(d) * matrix(F,[a.outer_product(b).list() 
                for a,b in zip(x.columns(),p.columns())]).T
        
    @cache
    def minor(cols):
        print ('minor',cols)
        if len(cols) == 0:
            return matrix(F,1,len(ms),[R.one()]*len(ms))
        i = n - len(cols)
        p = sum([(-1)**ji * multiply(m[i][j],minor(cols[:ji]+cols[ji+1:]),len(cols)-1) for ji,j in enumerate(cols) if not m[i][j].is_zero()])
        print (cols,p.nrows())
        return p
    return minor(tuple(range(n)))

def mult_maps(I):
    R = I.ring()
    F = R.base_ring()
    ltI = ideal([p.lt() for p in I.groebner_basis()])
    gbmaxd = max([p.degree() for p in I.groebner_basis()])
    xsix = {x:i for i,x in enumerate(R.gens())}

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
        return (tuple(sorted(monsout,reverse=True)),tuple(sorted(monsin,reverse=True)))
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
            return matrix(F,len(monsin),len(monsout),mat).T

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

        rels = block_matrix([[identity_matrix(F,len(monsin),sparse=True),
                                matrix(F,len(monsin),len(monsout),sparse=True)]])
        for i in range(len(monsin)):
            rels[i,i] = 1
        for yi,(divmap,(intoout,intoin)) in enumerate(zip(divmaps,mult_maps_noreduce(d-1))):
            print(yi,end=" ",flush=True)
            multmap = block_matrix([[-intoin],[intoout]])
            ixs = [i for i,j in multmap.dict().keys()]
            ixs.sort()
            jxs = [j for i,j in divmap.dict().keys()]
            jxs.sort()
            rhs = multmap[ixs,:]*reducemap(d-1)*divmap[:,jxs]
            print(rhs.dimensions(),RDF(rhs.density()))
            rels[jxs, ixs] = rhs.T
        print("solving",rels.dimensions())
        rels.echelonize()
        return rels[:,len(monsin):].T.dense_matrix()

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
    def mult_maps_reduce(d): # d -> d+1
        print('mult_maps_reduce',d)
        return [reducemap(d+1)*intoin + intoout for intoout,intoin in mult_maps_noreduce(d)]
    @cache
    def mult_map(d):
        return block_matrix([[m for m in mult_maps_reduce(d)]])
    return mult_map,mult_maps_reduce,mult_maps_noreduce,reducemap,mons


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

x0 = h.samp()
jact = jac.T
jact = random_matrix(F,16,16)*jact
jact=jact.apply_map(lambda e: e(x0)).augment(identity_matrix(16)).echelon_form()[:,-16:]*jact
jact = simplify_polynomial_matrix(jact)

my_reduce = reduce_fn_memo(I)
mult_map,mult_maps_reduce,mult_maps_noreduce,reducemap,mons = mult_maps(I)

# dat = upper_tri_assume_all_generic(jacm,I)

# vim: ft=python
