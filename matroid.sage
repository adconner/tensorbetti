# study betti decomposition of C3 C3 C3

# For T, equations of Tinv are given by Tdual after inv, saturated at the locus
# of indeterminacy of inv. View T as Tdual evaluated at codimension n-1 x n-1 minors.

# when n = 3, we are considering subspaces of dimension 6 of the 9 dimensional
# space spanned by 2x2 minors of a symbolic matrix X. In the span of n-1 by n-1 minors, 
# the orbit of an ordinary minor has projective dimension 2(n-1) = 4, and so
# intersects the 6 dimensional space Tinv (projective 5 dimensional) in a projective
# curve. Assume the projective curve is not contained in a hyperplane of Tinv
# (is this always true?), so we can pick a basis of Tinv consisting of ordinary
# minors X mapsto det(AXB), where A is 2x3 and B is 3x2. Representing such A
# and B by their perps, we can give the data of Tinv by 3 x 6 matrices A and B,
# each subject to the (independent) left actions of GL_3. Lets conjecture that
# the betti table depends only on the colinearity relations (representable
# matroid structure) among the columns of A and B

from IPython import embed
from itertools import *
from functools import cache
from sage.libs.singular.function_factory import ff

def is_concise(T):
    F = T[0].base_ring()
    import numpy as np
    T = np.array(T)
    for i in range(3):
        flat = matrix(F,T.reshape(T.shape[0],-1))
        if flat.rank() != flat.nrows():
            return False
        T = T.transpose((1,2,0))
    return True

def flattening_ranks(T):
    F = T[0].base_ring()
    import numpy as np
    T = np.array(T)
    rs = []
    for i in range(3):
        rs.append(matrix(F,T.reshape(T.shape[0],-1)).rank())
        T = T.transpose((1,2,0))
    return rs

def ABtops(A,B):
    n = A.nrows()
    F = A.base_ring()
    R = PolynomialRing(F,'x',n**2)
    X = matrix(R,n,n,R.gens())
    return [ (e.row().right_kernel_matrix()*X*f.row().right_kernel_matrix().T).det() for e,f in zip(A.columns(),B.columns()) ]

def pstoT(ps):
    R = ps[0].parent()
    n = sqrt(R.ngens())
    X = matrix(R,n,n,R.gens())
    ms = X.adjugate().list()
    T = ff.lift(ideal(ms), ideal(ps)).change_ring(R.base_ring()).T.right_kernel_matrix()
    T = [matrix(F,n,n,r.list()) for r in T]
    return T

def Ttops(T):
    n = T[0].nrows()
    R = PolynomialRing(T[0].base_ring(),'x',n**2)
    X = matrix(R,n,n,R.gens())
    Tinv = matrix([m.list() for m in T]).right_kernel_matrix()
    ps = [sum(a*b for a,b in zip(X.adjugate().list(),r)) for r in Tinv]
    return ps

def sat(ps,minors=False):
    R = ps[0].parent()
    n = sqrt(R.ngens())
    X = matrix(R,n,n,R.gens())
    if not minors:
        return ideal(ps).saturation(X.det())[0]
    else:
        return ideal(ps).saturation(ideal(X.adjugate().list()))[0]

def quotient_by_det(T):
    n = T[0].nrows()
    R = PolynomialRing(T[0].base_ring(),'x',n**2)
    X = matrix(R,n,n,R.gens())
    ps = X.adjugate().list()
    M2 = matrix(ff.syz(ideal(ps))).T
    T = matrix([m.list() for m in T])
    x = vector([R.gen(i*n) for i in range(n) ] + [0]*(n**2-n))
    K = (T*x.column()).augment(-T*M2)
    ss = ff.syz(K)
    return ideal([r[0] for r in ss])


def unsat(I):
    R = I.ring()
    n = sqrt(R.ngens())
    X = matrix(R,n,n,R.gens())
    ps = [p for p in I.intersection(ideal(X.adjugate().list())).gens() if p.degree() == n-1]
    return ps

def getAs(m,F=GF(32003),withrep=True):
    if m == 1:
        yield matrix(F,3,1,{(0,0):1})
        return
    def relsf(B):
        M = Matroid(B)
        return frozenset(M.dependent_r_sets(2)) | frozenset(M.dependent_r_sets(3))
    seen = set()
    for A in getAs(m-1,F,withrep):
        if withrep:
            for j in range(A.ncols()):
                B = A.augment(A[:,j],subdivide=False)
                rels = relsf(B)
                if rels not in seen:
                    seen.add(rels)
                    yield B
        # dont need to check seen here
        yield normalize(A.augment(random_vector(F,3),subdivide=False))
        for js in combinations(range(m-1),2):
            if A[:,js[0]] == A[:,js[1]]:
                continue
            B = normalize(A.augment(A[:,js]*random_vector(F,2),subdivide=False))
            rels = relsf(B)
            if rels not in seen:
                seen.add(rels)
                yield B
        for js in combinations(range(m-1),2):
            if A[:,js[0]] == A[:,js[1]]:
                continue
            for ks in combinations( [j for j in range(js[0]+1,m-1) if j != js[1]],2 ):
                if A[:,ks[0]] == A[:,ks[1]]:
                    continue
                K = block_matrix([ [A[:,js].T.right_kernel_matrix()],
                    [A[:,ks].T.right_kernel_matrix()] ]).right_kernel_matrix()
                B = normalize(A.augment(K.row(0),subdivide=False))
                if B.column(m-1) in B.columns()[:m-1]:
                    continue
                rels = relsf(B)
                if rels not in seen:
                    seen.add(rels)
                    yield B

def getAs_iso(m,F=GF(32003),withrep=True):
    As = list(getAs(m,F,withrep))
    def key(A):
        M = Matroid(A)
        out = [tuple([e+1 for e in sorted(s)]) 
            for s in chain(M.dependent_r_sets(2),M.dependent_r_sets(3))]
        return tuple(sorted(out))
    As = {key(A): A for A in As}
    os = gap.OrbitsDomain(gap.SymmetricGroup(m), list(As.keys()), gap.OnSetsSets)
    # return gap.Orbits(gap.SymmetricGroup(m), list(As.keys()), gap.OnSetsSets).sage(),As
    return [As[tuple(map(tuple,o[0]))] for o in os.sage()]

def normalize(A):
    A.echelonize()
    jxs = A.pivots()
    try:
        jnorm = max(range(A.ncols()), key=lambda j: len(A[:,j].nonzero_positions()))
    except ValueError:
        return A
    for i,j in enumerate(jxs):
        if A[i,jnorm] != 0:
            A[i] /= A[i,jnorm]
    for j in range(A.ncols()):
        A[:,j] /= next(e for e in A[:,j].list() if e != 0)
    return A

def getABs(m,F=GF(32003)):
    R = PolynomialRing(F,'x',9)
    X = matrix(R,3,3,R.gens())
    ms = ideal(X.adjugate().list())

    def key(A):
        M = Matroid(A)
        out = [tuple([e+1 for e in sorted(s)]) 
            for s in chain(M.dependent_r_sets(2),M.dependent_r_sets(3))]
        return tuple(sorted(out))

    out = []
    As = { key(A) : A for A in getAs(m,F)}
    Bs = { key(A) : A for A in getAs(m,F)}
    for o in gap.OrbitsDomain(gap.SymmetricGroup(m),list(As.keys()),gap.OnSetsSets).sage():
        Akey = o[0]
        Astab = gap.Stabilizer(gap.SymmetricGroup(m),Akey,gap.OnSetsSets)
        print(Astab.Size())
        for o in gap.OrbitsDomain(Astab,list(As.keys()),gap.OnSetsSets).sage():
            Bkey = o[0]

            A = As[tuple(map(tuple,Akey))]
            B = Bs[tuple(map(tuple,Bkey))]

            ps = ABtops(A,B)
            Tinv = ff.lift(ms,ideal(ps)).change_ring(F)
            if Tinv.rank() < m:
                continue
            out.append((A,B))
    # # gap.OrbitsDomain( list(res.keys())
    # return res
    return out




F = GF(32003)
A = random_matrix(F,3,6)
B = random_matrix(F,3,6)

# vim: ft=python
