from itertools import *
from collections import Counter
from copy import copy
from IPython import embed
from sage.libs.singular.function_factory import ff
from sage.libs.singular.option import opt,opt_ctx
load_attach_mode(load_debug=True)

# opt['prot'] = True
# opt['deg_bound'] = 3

F = GF(32003)
# params = 6+2
params = 18
R = PolynomialRing(F,['x%d' % i for i in range(9)]+['t%d' % i for i in range(params)],
                   order = TermOrder('degrevlex',9) + TermOrder('degrevlex',params))
# R = PolynomialRing(F,['x%d' % i for i in range(9)]+['t%d' % i for i in range(params)])
xs = R.gens()[:9]
ts = R.gens()[9:]
# S = FractionField(PolynomialRing(F,'t',params))
# ts = S.gens()
# R = PolynomialRing(S,'x',9)
# xs = R.gens()
qs = matrix(R,3,3,xs)[::-1,::-1].T.adjugate().list()


TM = matrix(R,6,9)
TM[:,:6] = identity_matrix(R,6)

# TM[:3,6:] = identity_matrix(R,3)
# TM[3:,6:] = identity_matrix(R,3)

TM[:,6:] = matrix(R,6,3,ts)
for i in range(3):
    for j in range(3):
        if i > j:
            TM[i,6+j] = 0
TM[0,6] = 1
TM[3,6] = 1
TM[0,-1] = 0

TMkeep = copy(TM)

TMg = copy(TM)
TMg[0,7] = 0
TMg[1,8] = 0

TM1 = copy(TMg)
TM1[4,6] = 1
TM1[5,6] = 1
TMs = [copy(TM1)]
TM1[4,6] = 0
TMs.append(copy(TM1))
TM1[5,6] = 0
TMs.append(copy(TM1))
TM1[4,6] = 1
TMs.append(copy(TM1))

## MORE needed 

# to test, should get basic one
for i in range(3):
    for j in range(3):
        if i != j:
            TMkeep[i,6+j] = 0
            TMkeep[3+i,6+j] = 0

# TM[:3,6:] = diagonal_matrix((0,)+ts[:2])

# TM[3,7] = 1
# TM[3,8] = 1
# TM[4:,6:] = matrix(R,2,3,ts[2:])

# TM[3,7] = 1
# TM[3,8] = 1
# TM[4:,6:] = matrix(R,2,3,ts[3:]+(0,))
# TM[3,6] = ts[2]

# TM = TMkeep
TM = TMs[2]
ps = [sum(a*b for a,b in zip(qs,r)) for r in TM]


xdeg = lambda p: -1 if p == 0 else sum(p.exponents()[0][:9])

def idealmatrix(ps,qs=[]):
    return vector(ps+qs).column().augment(
        block_matrix([[identity_matrix(len(ps)) ],[matrix(len(qs),len(ps))]]))
def allmonomials(xs,d):
    return [prod(ys) for ys in combinations_with_replacement(xs,d)]

# dbound = 6
# MM = idealmatrix(ps,allmonomials(xs,dbound+1))
# M = matrix(sorted([r for r in ff.groebner(MM.T) if not r[1:].is_zero()
#                       and xdeg(next(e for e in r[1:] if e!=0))+2 <= dbound
#                   ],reverse=True))

# M = idealmatrix(ps)
# M = matrix(ff.groebner(idealmatrix(ps).T))

J = R.ideal(ps)

def lmc(p):
    m = R.monomial(*tuple(p.lm().exponents()[0][:9])+(0,)*params)
    return (m,p.coefficient(m))

syz = [sum(a*b for a,b in zip(ps,r)) 
    for r in ff.syz(matrix([qs[:6]]))]
    # for r in ff.syz(matrix([[p(xs+(0,)*len(ts) ) for p in ps]]))]
syzcs = [[lmc(m)[1] for m in r.monomials()] for r in syz]

S = PolynomialRing(F,['x%d' % i for i in range(9)])

Jtarget=ideal([p(xs+(1,)*len(ts)) for p in J.gens()])
Jlt = ideal([p.lm() for p  in Jtarget.groebner_basis()])



def walk_ltIs(J,degbound=None):
    if degbound is not None:
        J = J + [prod(m) for m in combinations_with_replacement(xs,degbound+1)]
    J = R.ideal(ff.groebner(J))
    complexity = lambda J,m: (m.degree(), max(p.degree() for p in J.gens()),sum(p.degree() for p in J.gens()))
    from heapq import heappop, heappush
    pkey = lambda p: (-p.degree(), p.lm()) # use with descending order

    q = [(None,R.ideal(),[],R.one(),J)]
    # assumes ss not in I and ts not in I and I prime
    def walk_choices(I, ss, ts):
        if len(ts) == 0:
            yield (I,ss)
            return
        t = ts[0]
        tf = [f for f,_ in t.factor() if f not in ss]
        if len(tf) == 0:
            for Iout,ssout in walk_choices(I, ss, ts[1:]):
                yield (Iout,ssout)
            return
        for Iout,ssout in walk_choices(I, ss + tf, ts[1:]):
            yield (Iout,ssout)
        Is= (I + t).minimal_associated_primes()
        # print(len(Is),end=" ",flush=True)
        for I2 in Is:
            if all(s not in I2 for s in ss):
                for Iout,ssout in walk_choices(I2, ss, ts[1:]):
                    yield (Iout,ssout)
        
    Jcache = {}
    hilbertseries = set()
    while len(q) > 0:
        _, I, ss, mprev, J = heappop(q)
        if I.groebner_basis() not in Jcache:
            print(tuple(I.groebner_basis()),mprev,end=" ",flush=True)
            J = R.ideal(ff.std(J,I))
            lmcs = {}
            ltJ = []
            for p in J.gens():
                m,c = lmc(p)
                if m != 1:
                    ltJ.append(S.monomial(*m.exponents()[0][:S.ngens()]))
                lmcs.setdefault(m,[]).append(c)
            lmcs = sorted(lmcs.items(),key=lambda p: pkey(p[0]), reverse=True)
            ltJ = S.ideal(ltJ)
            hilbertseries.add(ltJ.hilbert_numerator())
            print(len(hilbertseries))
            Jcache[I.groebner_basis()] = (J,lmcs,ltJ)

        J,lmcs,_ = Jcache[I.groebner_basis()]
        try:
            m, cs = next( (m,cs) for m,cs in lmcs if pkey(m) < pkey(mprev) )
            if degbound is not None and m.degree() > degbound:
                continue
        except StopIteration:
            continue
        print (len(cs),end=" ",flush=True)
        for I2,ss2 in walk_choices(I,ss,cs):
            heappush(q,(complexity(I2,m), I2, ss2, m, J))

    return Jcache,hilbertseries

def walk_ltIs2(J):
    complexity = lambda J: (max(p.degree() for p in J.gens()),sum(p.degree() for p in J.gens()))
    q = [(None,R.ideal(),J)]
    from heapq import heappush, heappop
    checking = set()
    res = []
    hilbertseries = {}
    while len(q) > 0:
        _, I, J = heappop(q)
        J = R.ideal( (J+I).groebner_basis() )
        # J = R.ideal( ff.groebner(J+I) )
        res.append(J)
        lmcs = {}
        for r in J.gens():
            m,c = lmc(r)
            if m != 1:
                lmcs.setdefault(m,[]).append(c)
        ltJ = S.ideal([S.monomial(*m.exponents()[0][:S.ngens()]) for m in lmcs.keys()])
        hn = ltJ.hilbert_numerator()
        hilbertseries[hn] = J
        print(len(hilbertseries),list(Counter(map(xdeg,J.gens())).values()),tuple(I.gens()),hn)
        for cs in lmcs.values():
            for I2 in (I + cs).minimal_associated_primes():
                if I2 not in checking:
                    checking.add(I2)
                    heappush(q, (complexity(I2), I2, J))
    return res

def lmc_row(r):
    i,p = next((i,p) for i,p in enumerate(r) if p != 0)
    return (i,*lmc(p))

def walk_ltIs_module(M):
    complexity = lambda J: (max(p.degree() for p in J.gens()),sum(p.degree() for p in J.gens()))
    q = [(None,R.ideal(),M)]
    from heapq import heappush, heappop
    checking = set()
    res = []
    while len(q) > 0:
        _, I, M = heappop(q)
        M = block_matrix([[M], [vector(I.gens()).column().augment(
                                 matrix(R,I.ngens(),M.ncols()-1))] ])
        M = matrix(ff.groebner(M.T))
        print(M.dimensions(),tuple(I.gens()))
        res.append((M,I))
        lmcs = {}
        for r in M:
            i,m,c = lmc(r)
            lmcs.setdefault((i,m),[]).append(c)
        for cs in lmcs.values():
            for J in (I + cs).minimal_associated_primes():
                if J not in checking:
                    checking.add(J)
                    heappush(q, (complexity(J), J, M))
                    # heappush(q, (-J.dimension(), J, M))
    return res


# opt['deg_bound'] = 3

# M = matrix(R,2,3)
# M[:,1:] = identity_matrix(R,2)
# M[0,0] = xs[0]*xs[2]*ts[0]
# M[1,0] = xs[1]*xs[2]*ts[1]

# F = GF(32003)
# S = PolynomialRing(F,'t',3*9)
# R = PolynomialRing(S,'x',9)
# ps = matrix(R,3,3,R.gens()).adjugate().list()
# qs = [p*x for p in ps for x in R.gens()]
# ms = sorted([m for q in qs for m in q.monomials()],reverse=True)
# M = matrix(S,[[q.monomial_coefficient(m) for m in ms] for q in qs])

def echelon_forms(m,I=None):
    from collections import Counter

    S = m.base_ring()
    if I is None: I = S.ideal()

    R = PolynomialRing(S.base_ring(),S.gens()+('d',))
    dv = R.gens()[-1]

    def rec(m,I,d,r,jxs):
        if m[r:,r:].is_zero():
            yield (m,I,jxs,d)
            return

        # select the most common nonzero element of 
        f = Counter(m[r:,r:].dict().values()).most_common(1)[0][0]
        print (r,f)

        # Case 1: pass to the closed set where f == 0
        Icur = I + f
        if 1 not in Icur:
            mcur,rnext,jxsnext = elimination_by_units_jxs(m,Icur,r,jxs)
            for res in rec(mcur,Icur,d,rnext,jxsnext):
                yield res

        # Case 2: pass to the open set where f not identically zero, e.g., localize at f
        mcur = m.subs({dv : dv*f})
        Icur = I.elimination_ideal(dv) + [dv*d*f-1]
        if 1 not in Icur:
            mcur,rnext,jxsnext = elimination_by_units_jxs(mcur,Icur,r,jxs)
            for res in rec(mcur,Icur,d*f,rnext,jxsnext):
                yield res

    jxs = list(range(m.ncols()))
    m,r,jxs = elimination_by_units_jxs(m,I,0,jxs)
    return m,r
    print('begin',r)
    return rec(m.change_ring(R), I.change_ring(R), R.one(), r, jxs)


# m : matrix
# 
# row reduces m as far as possible pivoting using unit entries. Column swaps are
# also performed. m is modified in place and the number of pivots r is returned. 
# After this operation, m has the form
# [ T A ]
# [ 0 B ]
# where T is r by r and  upper triangular with units along the diagonal, B 
# contains no units
def elimination_by_units_jxs(m,I,r,jxs):
    m = m.apply_map(lambda e: I.reduce(e))

    jxs = copy(jxs)

    # computes the inverse mod $I$ if it exists, otherwise None
    def try_inverse(e):
        try:
            return ff.lift(I+e,1).list()[-1]
        except RuntimeError:
            return None

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
        jxs[r], jxs[j] = jxs[j], jxs[r]
        m[r,:] *= einv
        for k in m.nonzero_positions_in_row(r):
            m[r,k] = I.reduce(m[r,k])
        assert m[r,r] == 1
        for i in m.column(r).nonzero_positions():
            if i == r:
                continue
        # for i in m.column(r)[r+1:].nonzero_positions():
        #     i += r+1
            m.add_multiple_of_row(i,r,-m[i,r])
            for k in m.nonzero_positions_in_row(i):
                m[i,k] = I.reduce(m[i,k])
        r += 1

    return m,r,jxs
  

# vim: ft=python

