from itertools import *
from sage.libs.singular.option import opt, opt_ctx
import sage.libs.singular.function_factory

sing = sage.libs.singular.function_factory.ff
std = sing.std
syz = sing.syz
res = sing.res
betti = sing.betti

# see libSingular: Options sagemath

def main():
    n = 5
    r = 10
    dim = 16

    # F = GF(5)
    # F = GF(101)
    # F = GF(1031)
    F = GF(32003)
    # F = GF(65537)

    # opt['prot'] = True
    return generic_project(Tinv(random_tensor(F,n,r)),dim)

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
        print("%d monomials undetermined" % len(ms))
        mis, ms = [mi for mi, _ in ms], [m for _, m in ms]
        eqs = matrix(self.F, [self.sampm(mis) for i in range(len(mis))])
        pscur = [sum(a*m for a,m in zip(r,ms)) for r in eqs.right_kernel_matrix()]
        gbto = (Ilower + pscur).groebner_basis(deg_bound=d)
        return self.R.ideal(gbto)

    @memoize
    def Id_mod_lower_basis(self, d):
        I = self.ideal_to(d - 1)
        ltlower = ideal([p.lm() for p in I.groebner_basis(deg_bound=d)])
        return [p for p in self.ideal_to(d).gens() if p.lm() not in ltlower]

    def write_to(self, deg, pre="examples/"):
        open("%sT.txt" % pre, "w").write(str([map(list, m) for m in self.T]) + "\n")
        for d in range(deg + 1):
            ps = self.Id_mod_lower_basis(d)
            if len(ps) == 0:
                continue
            open("%sI%d.txt" % (pre, d), "w").write(",\n".join(map(str, ps)) + "\n")

def Tinv(T):
    F = T[0].base_ring()
    n = T[0].nrows()
    def samp():
        return sum(e * m for e, m in zip(random_vector(F, n), T)).adjugate().list()
    return ParameterizedVariety(samp)

def generic_project(V,dim):
    P = random_matrix(V.F,dim,V.n)
    def samp():
        return (P*vector(V.samp())).list()
    return ParameterizedVariety(samp)

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


if __name__ == "__main__":
    h = main()

# vim: ft=python
