from functools import cache
import numpy as np

class ParameterizedVariety:
    def __init__(self, samp):
        self.samp = samp
        v = samp()
        self.n = len(v)
        self.F = v[0].parent()
        self.p = self.F.characteristic()
        self.R = PolynomialRing(self.F,'x',self.n)

    def sampm(self, mis, cnt=1):
        if len(mis) == 0:
            return matrix(self.F,cnt,0)
        samps = []
        for _ in range(cnt):
            s = self.samp()
            samps.append([prod(s[i] for i in mi) for mi in mis])
        return matrix(self.F,samps)

    def sampm_numpy(self, mis, cnt=1):
        assert self.p > 0 and self.p < 2**31 and self.F == GF(self.p)
        dtype = np.int16 if self.p < 2**7 else np.int32 if self.p < 2**15 else np.int64
        if len(mis) == 0:
            return np.zeros((cnt,0),dtype=dtype)
        samps = np.array([self.samp() for _ in range(cnt)],dtype=dtype)
        mis = np.array(mis)
        res = samps[:,mis[:,0]]
        for i in range(1,mis.shape[1]):
            # print(i,end=' ',flush=True)
            res *= samps[:,mis[:,i]]
            res %= self.p
        return res

    def sampp(self, ps, cnt=1):
        mis = []
        misix = {}
        P = [[] for _ in ps]
        for pi,p in enumerate(ps):
            for m in p.monomials():
                mi = tuple(i for i,e in m.exponents()[0].sparse_iter() for _ in range(e))
                if mi not in misix:
                    misix[mi] = len(mis)
                    mis.append(mi)
                P[pi].append((misix[mi],int(p.monomial_coefficient(m))))
        try:
            samp = self.sampm_numpy(mis,cnt)
            P = [np.array(p) for p in P]
            eqs = matrix(self.F,[ np.sum( (samp[:,p[:,0]]*p[:,1]) % self.p, axis = 1) %
                self.p for p in P]).T
            return eqs
        except AssertionError:
            samp = self.sampm(mis,cnt)
            P = matrix(self.F,len(ps),len(mis),{(i,j) : e  for i,p in enumerate(P) for j,e in p})
            return samp * P.T

    @cache
    def ideal_to(self, d):
        if d <= 0:
            return self.R.ideal()
        Ilower = self.ideal_to(d - 1)
        print("getting component of ideal in degree %d" % d)
        Ilower = self.R.ideal(Ilower.groebner_basis(deg_bound=d))
        ltI = [tuple(p.lm().exponents()[0]) for p in Ilower.gens() if not p.is_zero()]
        ms = [(tuple(i for i,e in enumerate(m) for _ in range(e)) ,self.R.monomial(*m)) 
                for m in monomial_ideal_complement(self.R.ngens(), d, ltI)]
        # ltI = self.R.ideal([p.lm() for p in Ilower.gens() if not p.is_zero()])
        # ms = [
        #     (mi, m)
        #     for mi in combinations_with_replacement(range(self.R.ngens()), d)
        #     for m in [prod(self.R.gen(i) for i in mi)]
        #     if m not in ltI
        # ]
        print("%d monomials undetermined, " % len(ms),end="",flush=True)
        mis, ms = [mi for mi, _ in ms], [m for _, m in ms]
        eqs = matrix(self.F,self.sampm_numpy(mis, len(mis)))
        pscur = [sum(a*m for a,m in zip(r,ms)) for r in eqs.right_kernel_matrix()]
        print ("%d relations found" % len(pscur))
        if len(pscur) == 0:
            return Ilower
        else:
            return self.R.ideal((Ilower + pscur).groebner_basis(deg_bound=d))
    
    @cache
    def ideal_to_intersect(self, d, J=None):
        if J is None:
            J = self.R.ideal(1)
        if d <= 0:
            return self.R.ideal()
        Ilower = self.ideal_to_intersect(d - 1, J)
        print("getting component of ideal in degree %d" % d)
        ltI = [tuple(p.lm().exponents()[0]) for p in Ilower.gens() if not p.is_zero()]
        ltJ = [tuple(p.lm().exponents()[0]) for p in J.groebner_basis() if not p.is_zero()]
        ltJd = [self.R.monomial(*m) for m in monomial_ideal_component(ltJ, d, ltI)]
        ps = [m - J.reduce(m) for m in ltJd]
        ps.sort(reverse=True)
        print("%d monomials undetermined, " % len(ps),end="",flush=True)
        eqs = self.sampp(ps,len(ps))
        pscur = [sum(a*p for a,p in zip(r,ps)) for r in eqs.right_kernel_matrix()]
        print ("%d relations found" % len(pscur))
        return Ilower + pscur

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
