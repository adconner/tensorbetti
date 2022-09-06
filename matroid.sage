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

    
F = GF(32003)
A = random_matrix(F,3,6)
B = random_matrix(F,3,6)

# vim: ft=python
