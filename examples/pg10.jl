using Pkg
Pkg.activate(@__DIR__)

using Crystalline, TETB

pg_num = 10
brs = calc_bandreps(pg_num, Val(2))

# needs to do that to find the WPs properly
coefs = zeros(Int, length(brs))
coefs[1] = 1
coefs[5] = 1

cbr = CompositeBandRep(coefs, brs)

Rs = [[0, 0]]

## -------------------------------------------------------------------------------------- ##

sym_hop = TETB.obtain_symmetry_related_hoppings(Rs, cbr.brs[1], cbr.brs[5])

## -------------------------------------------------------------------------------------- ##