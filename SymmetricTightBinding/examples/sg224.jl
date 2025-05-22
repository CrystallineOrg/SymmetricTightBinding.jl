using Pkg
Pkg.activate(@__DIR__)

using Crystalline, SymmetricTightBinding

sgnum, D = 224, 3
brs = calc_bandreps(sgnum, Val(D))
cbr = @composite brs[13] + brs[19]
Rs = [[0, 0, 0]]

tb_model = tb_hamiltonian(cbr, Rs)

hops = obtain_symmetry_related_hoppings(Rs, brs[13], brs[19])
