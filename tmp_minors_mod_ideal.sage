load('all.sage')

a = 4
n = 4
r = 5

# F = GF(5)
# F = GF(101)
# F = GF(1031)
F = GF(32003)
# F = GF(65537)

h = Tinv(random_tensor(F,n,r))
I = h.ideal_to(2)
jac = matrix([[p.derivative(x) for x in I.ring().gens()] for p in I.gens()[::-1]])

x0 = h.samp()
jact = jac.T
jact = random_matrix(F,16,16)*jact
jact=jact.apply_map(lambda e: e(x0)).augment(identity_matrix(16)).echelon_form()[:,-16:]*jact
jact = simplify_polynomial_matrix(jact)

my_reduce = reduce_fn_memo(I)
mmaps = mult_maps(I)
mons,reducemap,mult_map = mmaps

# dat = upper_tri_assume_all_generic(jacm,I)

# vim: ft=python
