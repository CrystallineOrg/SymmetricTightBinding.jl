using Pkg
Pkg.activate(@__DIR__)

using Crystalline, TETB

sg_num, D = 224, 3
brs = calc_bandreps(sg_num, Val(D))
c = zeros(Int, length(brs))
c[13], c[19] = 1, 1
cbr = CompositeBandRep(c, brs)
Rs = [[0, 0, 0]]

tb_model = tb_hamiltonian(cbr, Rs, time_reversal=false)

hops = TETB.obtain_symmetry_related_hoppings(Rs, brs[13], brs[19])