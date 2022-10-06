load('all.sage')
F = GF(32003)

n = 4
a = n
r = 100
m = n-1

T = random_tensor(F,a,n,r)
h = Tinv(T)
# h = Tprimal(T)

# I = h.ideal_to(4)
I = h.ideal_to_intersect_minors(6,m)

print(ff.betti(ff.res(I,1)))

# vim: ft=python
