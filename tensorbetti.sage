from itertools import *
from sage.libs.singular.option import opt, opt_ctx
import sage.libs.singular.function_factory
from sys import stdout

load('ParameterizedVariety.sage')

from sage.libs.singular.function_factory import ff
sing = ff
std = sing.std
syz = sing.syz
res = sing.res
betti = sing.betti

# see libSingular: Options sagemath

# 4 6 11
# [  1   0   0]
# [  0   0   0]
# [  0  66 180]
# [  0   0 705]
# [  0   0 224]
# [  0   0  64]

# 4 6 13
# [   1    0    0]
# [   0    7    0]
# [   0  144 1153]
# [   0    0  678]
# [   0    0  406]

# 4 6 15
# [  1   0   0]
# [  0  36  80]
# [  0   0 893]
# [  0   0 565]
# [  0   0  78]

# 4 6 16
# [  1   0   0]
# [  0  52 236]
# [  0   0 778]
# [  0   0 313]
# [  0   0  45]

# 4 6 noproj
# [  1   0   0]
# [  0  52 236]
# [  0   0 801]
# [  0   0 147]
# [  0   0   2]

# [1 0 0 0]
# [0 3 0 0]
# [0 0 3 0]
# [0 0 0 1]


# [ 1  0  0  0  0]
# [ 0  0  0  0  0]
# [ 0  4  0  0  0]
# [ 0  1  0  0  0]
# [ 0  0  6  0  0]
# [ 0  0 16 30 12]

# n=4; h=Tinv(dual_tensor([matrix(GF(32003),n,n,{(i,i):1}) for i in range(n)]+[matrix(GF(32003),
# n,n,lambda i,j: 1)])); I=h.ideal_to(4); betti(res(I,0))
# [  1   0   0   0   0   0   0   0   0]
# [  0   0   0   0   0   0   0   0   0]
# [  0   5   0   0   0   0   0   0   0]
# [  0   5   0   0   0   0   0   0   0]
# [  0   0  26   0   0   0   0   0   0]
# [  0   0  50 134  70   0   0   0   0]
# [  0   0   0 100 288 316 160  36   0]
# [  0   0   0   0   0   0   0   0   1]


def main():
    a = 3
    n = 3
    r = 7
    dim = 16

    # F = GF(5)
    # F = GF(101)
    # F = GF(1031)
    F = GF(32003)
    # F = GF(65537)

    # opt['prot'] = True

    # return dual_tensor(random_tensor(F,n,r))
    return Tinv(random_tensor(F,n,r))
    # return Tinv(dual_tensor(random_tensor(F,n,r)))

    # F4 = Fn_poly(F,n**2,4)
    # return gauss_map(random_tensor(F,n,r),F4)

    # return generic_project(Tinv(dual_tensor(random_tensor(F,n,r))),dim)

    # for n in range(2,6):
    #     for a in range(max(n**2-3,1),n**2+1):
    #         print ("tensor size %d %d %d" % (a,n,n))
    #         I = Tinv(random_tensor_gen(F,a,n,n,100)).ideal_to(4)
    #         # print ('finding betti table')
    #         print (betti(res(I,0)))

    # I = Tinv(random_tensor_gen(F,n**2-a,n,n,100)).ideal_to(4)
    # print (betti(res(I,0)))

    # for r in range(4,10):
    #     print (n,r)
    #     I = Tinv(dual_tensor(random_tensor_gen(F,a,n,n,r))).ideal_to(4)
    #     print (betti(res(I,0)))


# Given a tensor, returns the parameterized variety Tinv
def Tinv(T):
    F = T[0].base_ring()
    n = len(T)
    def samp():
        return sum(e * m for e, m in zip(random_vector(F, n), T)).adjugate().list()
    return ParameterizedVariety(samp)

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


def random_tensor(F, n, r):
    return sum_rank_ones(
        [tuple(random_vector(F, n) for j in range(3)) for i in range(r)]
    )

def random_tensor_gen(F, a, b, c, r):
    return sum_rank_ones([(random_vector(F, a), random_vector(F,b), random_vector(F,c)) 
            for i in range(r)])

def dual_tensor(T):
    K = matrix([m.list() for m in T]).right_kernel_matrix()
    return [matrix(K.base_ring(),T[0].nrows(),T[0].ncols(),list(r)) for r in K] 

if __name__ == "__main__":
    h = main()

# vim: ft=python
