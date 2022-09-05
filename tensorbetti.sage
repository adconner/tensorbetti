from itertools import *
from sage.libs.singular.option import opt, opt_ctx
import sage.libs.singular.function_factory
from sys import stdout

sing = sage.libs.singular.function_factory.ff
std = sing.std
syz = sing.syz
res = sing.res
betti = sing.betti

# see libSingular: Options sagemath

def main():
    n = 4
    r = 4
    dim = 16

    # F = GF(5)
    # F = GF(101)
    # F = GF(1031)
    F = GF(32003)
    # F = GF(65537)

    # opt['prot'] = True

    # return dual_tensor(random_tensor(F,n,r))
    # return Tinv(random_tensor(F,n,r))
    # return Tinv(dual_tensor(random_tensor(F,n,r)))

    F4 = Fn_poly(F,n**2,4)
    return gauss_map(random_tensor(F,n,r),F4)

    # return generic_project(Tinv(dual_tensor(random_tensor(F,n,r))),dim)

    # for r in range(4,10):
    #     print (n,r)
    #     I = Tinv(dual_tensor(random_tensor(F,n,r))).ideal_to(4)
    #     print (betti(syz(I)))


def memoize(obj):
    cache = obj.cache = {}

    import functools

    @functools.wraps(obj)
    def memoizer(*args, **kwargs):
        if args not in cache:
            cache[args] = obj(*args, **kwargs)
        return cache[args]

    return memoizer

class ParameterizedVariety:
    def __init__(self, samp):
        self.samp = samp
        v = samp()
        self.n = len(v)
        self.F = v[0].parent()
        self.R = PolynomialRing(self.F,'x',self.n)

    def sampm(self, mis):
        s = self.samp()
        return [prod(s[i] for i in mi) for mi in mis]

    @memoize
    def ideal_to(self, d):
        if d < 0:
            return self.R.ideal()
        Ilower = self.ideal_to(d - 1)
        print("getting component of ideal in degree %d" % d)
        ltI = ideal([p.lm() for p in Ilower.groebner_basis(deg_bound=d)])
        ms = [
            (mi, m)
            for mi in combinations_with_replacement(range(self.R.ngens()), d)
            for m in [prod(self.R.gen(i) for i in mi)]
            if m not in ltI
        ]
        print("%d monomials undetermined, " % len(ms),end="")
        stdout.flush()
        mis, ms = [mi for mi, _ in ms], [m for _, m in ms]
        eqs = matrix(self.F, [self.sampm(mis) for i in range(len(mis))])
        pscur = [sum(a*m for a,m in zip(r,ms)) for r in eqs.right_kernel_matrix()]
        print ("%d relations found" % len(pscur))
        gbto = (Ilower + pscur).groebner_basis(deg_bound=d)
        return self.R.ideal(gbto)

    @memoize
    def Id_mod_lower_basis(self, d):
        I = self.ideal_to(d - 1)
        ltlower = ideal([p.lm() for p in I.groebner_basis(deg_bound=d)])
        return [p for p in self.ideal_to(d).gens() if p.lm() not in ltlower]

    def write_to(self, deg, pre="examples/"):
        # open("%sT.txt" % pre, "w").write(str([map(list, m) for m in self.T]) + "\n")
        for d in range(deg + 1):
            ps = self.Id_mod_lower_basis(d)
            if len(ps) == 0:
                continue
            open("%sI%d.txt" % (pre, d), "w").write(",\n".join(map(str, ps)) + "\n")

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
