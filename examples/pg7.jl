using Pkg
Pkg.activate(@__DIR__)

using Crystalline, TETB

pg_num, D = 7, 2
brs = calc_bandreps(pg_num, Val(D))
coefs = zeros(Int, length(brs))
coefs[end-1], coefs[end] = 1, 1  # pick (2a|A) and (2a|B) EBR
cbr = CompositeBandRep(coefs, brs)
Rs = [[0, 0]]
tbs = TETB.tb_hamiltonian(cbr, Rs)

# ------------------------------------------------- #

hops = obtain_symmetry_related_hoppings(Rs, cbr.brs[end-1], cbr.brs[end-1])