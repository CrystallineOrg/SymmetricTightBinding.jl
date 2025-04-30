using Pkg
Pkg.activate(@__DIR__)

using Crystalline, TETB

# This example shows a non-inversion symmetric case with 2D complex site-symmetry representations
pgnum, D = 13, 2
brs = calc_bandreps(pgnum, Val(D))
coefs = zeros(length(brs))
coefs[[1, 3]] .= 1
cbr = CompositeBandRep(coefs, brs)

br₁ = brs[1]
br₂ = brs[3]

# I will study here the term involving br₁ → br₂, which will be a non-diagonal term

Rs = [[0, 0]]

hops = obtain_symmetry_related_hoppings(Rs, br₁, br₂, true)