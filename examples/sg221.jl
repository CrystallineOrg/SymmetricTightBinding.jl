using Pkg
Pkg.activate(@__DIR__)

using Crystalline, TETB

sg_num = 221
brs = calc_bandreps(sg_num)

coefs = zeros(Int, length(brs))
coefs[6] = 1
cbr = CompositeBandRep(coefs, brs)

or = TETB.hamiltonian_term_order(cbr.brs[6], cbr.brs[6])

Rs = [[0, 0, 0]]

tbs = tb_hamiltonian(cbr, Rs)

# show the hopping orbitals that we are using to build the model

hops = obtain_symmetry_related_hoppings(Rs, cbr.brs[6], cbr.brs[6])

hop_plot(hops[2])