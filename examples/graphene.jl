using Pkg
Pkg.activate(@__DIR__)
using Crystalline, TETB

pg_num = 17
brs = calc_bandreps(pg_num, Val(2))
coefs = zeros(Int, length(brs))
coefs[5] = 1
cbr = CompositeBandRep(coefs, brs)

tb_model = TETB.tb_hamiltonian(cbr, [[0, 0]])

# visualize graphene in our set-up
a₁ = [1, 0]
a₂ = [-1 / 2, sqrt(3) / 2]
r₁ = [1 / 3, 2 / 3]
r₂ = [2 / 3, 1 / 3]

# hopping vectors can be obtained with the function
hop_vecs = TETB.obtain_symmetry_related_hoppings([[0, 0]], brs[5], brs[5])