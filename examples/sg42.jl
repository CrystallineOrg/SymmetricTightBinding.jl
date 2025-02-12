using Pkg
Pkg.activate(@__DIR__)

using Crystalline, TETB

sg_num = 42
brs = calc_bandreps(sg_num)
coefs = [1, 0, 1, 0, 0, 0]
cbr = CompositeBandRep(coefs, brs)

Rs = [[0, 0, 0]]

tb_model = tb_hamiltonian(cbr, Rs)

hops = obtain_symmetry_related_hoppings(Rs, brs[1], brs[3])

hop_plot(hops[1])