
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

# vim: ft=python
