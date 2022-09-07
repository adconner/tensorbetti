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

def sat(ps):
    R = ps[0].parent()
    n = sqrt(R.ngens())
    X = matrix(R,n,n,R.gens())
    return ideal(ps).saturation(ideal(X.adjugate().list()))[0]

def unsat(I):
    R = I.ring()
    n = sqrt(R.ngens())
    X = matrix(R,n,n,R.gens())
    ps = [p for p in I.intersection(ideal(X.adjugate().list())).gens() if p.degree() == n-1]
    return ps

def getAs_norep(m,F=GF(32003)):
    if m == 1:
        yield random_matrix(F,3,1)
        return
    seen = set()
    for A in getAs_norep(m-1,F):
        yield A.augment(random_vector(F,3),subdivide=False)
        for js in combinations(range(m-1),2):
            B = A.augment(A[:,js]*random_vector(F,2),subdivide=False)
            rels = frozenset(Matroid(B).dependent_r_sets(3))
            if rels not in seen:
                seen.add(rels)
                yield B
        for js in combinations(range(m-1),2):
            for ks in combinations( [j for j in range(js[0]+1,m-1) if j != js[1]],2 ):
                K = block_matrix([ [A[:,js].T.right_kernel_matrix()],
                    [A[:,ks].T.right_kernel_matrix()] ]).right_kernel_matrix()
                B = A.augment(F.random_element()*K.row(0),subdivide=False)
                C = copy(B)
                for j in range(C.ncols()):
                    C[:,j] /= next(e for e in C[:,j].list() if e!=0)
                if C.column(m-1) in C.columns()[:m-1]:
                    continue
                rels = frozenset(Matroid(B).dependent_r_sets(3))
                if rels not in seen:
                    seen.add(rels)
                    yield B

def getAs_norep_iso(m,F=GF(32003)):
    As = list(getAs_norep(m))
    def key(A):
        return tuple(sorted([tuple([e+1 for e in sorted(s)]) 
            for s in Matroid(A).dependent_r_sets(3)]))
    As = {key(A): A for A in As}
    os = gap.OrbitsDomain(gap.SymmetricGroup(m), list(As.keys()), gap.OnSetsSets)
    return [As[tuple(map(tuple,o[1].sage()))] for o in os]

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

F = GF(32003)
A = random_matrix(F,3,6)
B = random_matrix(F,3,6)

# vim: ft=python
