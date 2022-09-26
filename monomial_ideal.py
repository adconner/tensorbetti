
def monomial_ideal_complement(n, deg, ms):
    if deg < 0:
        return []
    def search(lb,ub,slb,i):
        if i == len(ms):
            for pt in IntegerVectors(deg-slb,n,outer=[u-l for l,u in zip(lb,ub)]):
                yield tuple(e+f for e,f in zip(lb,pt))
            return
        # print(lb,ub,i)
        m = ms[i]
        lb = copy(lb)
        for j,e in enumerate(m):
            # this is the first position smaller than m
            if e-1 >= lb[j]:
                ublast = ub[j]
                ub[j] = min(ub[j],e-1)
                for pt in search(lb,ub,slb,i+1):
                    yield pt
                ub[j] = ublast
            if e > lb[j]:
                slb += e - lb[j]
                lb[j] = e
                if lb[j] > ub[j]:
                    return
                if slb > deg:
                    return
    return search([0]*n,[deg]*n,0,0)

def monomial_ideal_component(ms, deg, less=[]):
    n = len(ms[0])
    less = copy(less)
    for m in ms:
        for pt in monomial_ideal_complement(n, deg - sum(m), 
                [[max(e-f,0) for e,f in zip(mp,m)] for mp in less]):
            yield tuple(e+f for e,f in zip(m,pt))
        less.append(m)
