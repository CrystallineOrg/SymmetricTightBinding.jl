using Pkg
Pkg.activate(@__DIR__)

using Crystalline, TETB

sgnum, D = 224, 3
brs = calc_bandreps(sgnum, Val(D))
c = zeros(Int, length(brs))
c[13], c[19] = 1, 1
cbr = CompositeBandRep(c, brs)
Rs = [[0, 0, 0]]

tb_model = tb_hamiltonian(cbr, Rs; timereversal = false)

hops = TETB.obtain_symmetry_related_hoppings(Rs, brs[13], brs[19])
