using Pkg
Pkg.activate(@__DIR__)

using Crystalline, TETB

sg_num = 221
brs = calc_bandreps(sg_num)

coefs = zeros(Int, length(brs))
coefs[6] = 1
cbr = CompositeBandRep(coefs, brs)

Rs = [[0, 0, 0]]

tbs = TETB.tb_hamiltonian(cbr, Rs)

# it only obtains (or at least prints) the onsite terms. Does it only print for the first 
# term in Rs?

hops = TETB.obtain_symmetry_related_hoppings([[0, 0, 0]], cbr.brs[6], cbr.brs[6])