using Pkg
Pkg.activate(@__DIR__)

using Crystalline, SymmetricTightBinding, GLMakie

sgnum = 42
brs = calc_bandreps(sgnum)
cbr = @composite brs[1] + brs[3]

Rs = [[0, 0, 0]]

tb_model = tb_hamiltonian(cbr, Rs)

hops = obtain_symmetry_related_hoppings(Rs, brs[1], brs[3])

plot(hops[1])
