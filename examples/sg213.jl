using Pkg
Pkg.activate(@__DIR__)

using Crystalline, TETB

sg_num = 213
brs = calc_bandreps(sg_num)
coefs = zeros(Int, length(brs))
coefs[end] = 1
cbr = CompositeBandRep(coefs, brs)

or = TETB.hamiltonian_term_order(cbr.brs[end], cbr.brs[end])

Rs = [[0, 0, 0]]

tb_model = tb_hamiltonian(cbr, Rs)

# show the hopping orbitals that we are using to build the model

hops = obtain_symmetry_related_hoppings(Rs, cbr.brs[end], cbr.brs[end])