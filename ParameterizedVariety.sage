from functools import cache

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

    @cache
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

    @cache
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

# vim: ft=python
