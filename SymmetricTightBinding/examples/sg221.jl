using Pkg
Pkg.activate(@__DIR__)

using Crystalline, SymmetricTightBinding, GLMakie

sgnum = 221
brs = calc_bandreps(sgnum)
cbr = @composite brs[6]

ordering = SymmetricTightBinding.OrbitalOrdering(cbr.brs[6])

Rs = [[0, 0, 0]]

tbs = tb_hamiltonian(cbr, Rs)

# show the hopping orbitals that we are using to build the model

hops = obtain_symmetry_related_hoppings(Rs, cbr.brs[6], cbr.brs[6])

plot(hops[2])
