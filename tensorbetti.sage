

def main():
    n = 5
    r = 8

    F = GF(5)
    # F = GF(101)
    # F = GF(1031)
    # F = GF(32003)
    # F = GF(65537)

    
    for r in [5..9]:
        for run in range(3):
            T = random_tensor(F, n, r)
            # T=[random_matrix(F,n,n) for i in range(n)]
            # T=[random_matrix(ZZ,n,n,x=-100,y=100) for i in range(n)]

            # T = [m.change_ring(F) for m in matrixmult(2,2,2)]

            h = TensorBetti(T)
            h.write_to(3,'examples/r%drun%d_' % (r,run))

from itertools import *
from sage.libs.singular.option import opt, opt_ctx
import sage.libs.singular.function_factory

syz = sage.libs.singular.function_factory.ff.syz
std = sage.libs.singular.function_factory.ff.std

# see libSingular: Options sagemath


def memoize(obj):
    cache = obj.cache = {}

    import functools

    @functools.wraps(obj)
    def memoizer(*args, **kwargs):
        if args not in cache:
            cache[args] = obj(*args, **kwargs)
        return cache[args]

    return memoizer


class TensorBetti:
    def __init__(self, T):
        self.T = T
        self.n = T[0].nrows()
        self.F = T[0].base_ring()
        self.R = PolynomialRing(
            FractionField(self.F),
            ["m%d%d" % (i, j) for i, j in product(range(self.n), range(self.n))],
        )

    def samp(self):
        return (
            sum(e * m for e, m in zip(random_vector(self.F, self.n), self.T))
            .adjugate()
            .list()
        )

    def sampm(self, mis):
        s = self.samp()
        return [prod(s[i] for i in mi) for mi in mis]

    def mat_to_ps(self, mat, ms):
        return [sum(a * m for a, m in zip(r, ms)) for r in mat]

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
        pscur = self.mat_to_ps(eqs.right_kernel_matrix(), ms)
        gbto = (Ilower + pscur).groebner_basis(deg_bound=d)
        return ideal(gbto)

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
    main()

# vim: ft=python
