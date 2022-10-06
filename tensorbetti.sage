
def Tinv(T):
    F = T[0].base_ring()
    n = len(T)
    def samp():
        return sum(e * m for e, m in zip(random_vector(F, n), T)).adjugate().list()
    return ParameterizedVariety(samp)

def Tprimal(T):
    F = T[0].base_ring()
    n = len(T)
    def samp():
        return sum(e * m for e, m in zip(random_vector(F, n), T)).list()
    return ParameterizedVariety(samp)

def random_tensor(F, a, n, r):
    return random_tensor_gen(F, a, n, n, r)

def random_tensor_gen(F, a, b, c, r):
    return sum_rank_ones([(random_vector(F, a), random_vector(F,b), random_vector(F,c)) 
            for i in range(r)])

# returns the parameterized variety of the gauss map of p applied to the linear
# subspace
def gauss_map(T,p):
    F = T[0].base_ring()
    a = len(T)
    partials = [p.derivative(x) for x in p.parent().gens()]
    def samp():
        m = sum(e * m for e, m in zip(random_vector(F, a), T)).list()
        return [q(m) for q in partials]
    return ParameterizedVariety(samp)

def generic_project(V,dim):
    P = random_matrix(V.F,dim,V.n)
    def samp():
        return (P*vector(V.samp())).list()
    return ParameterizedVariety(samp)

def Fn_poly(F,dim,n):
    R = PolynomialRing(F,'x',dim)
    return sum(prod(m) for m in combinations_with_replacement(R.gens(),n-1))

def sum_rank_ones(rank1s, sparse=True):
    from operator import add

    def rank_one(a, b, c):
        return [e * b.column() * c.row() for i, e in enumerate(a)]

    res = list(reduce(lambda a, b: map(add, a, b), [rank_one(*m) for m in rank1s]))
    return res

def dual_tensor(T):
    K = matrix([m.list() for m in T]).right_kernel_matrix()
    return [matrix(K.base_ring(),T[0].nrows(),T[0].ncols(),list(r)) for r in K] 

# vim: ft=python
