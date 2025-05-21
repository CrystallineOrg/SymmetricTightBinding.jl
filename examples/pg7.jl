using Pkg
Pkg.activate(@__DIR__)

using Crystalline, TETB

sgnum, D = 7, 2
brs = calc_bandreps(sgnum, Val(D))
cbr = @composite brs[end-1] + brs[end] # pick (2a|A) and (2a|B) EBR
tbs = TETB.tb_hamiltonian(cbr)

# ------------------------------------------------- #

hops = obtain_symmetry_related_hoppings(Rs, cbr.brs[end-1], cbr.brs[end-1])
