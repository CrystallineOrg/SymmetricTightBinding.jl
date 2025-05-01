using Pkg
Pkg.activate(@__DIR__)

using Crystalline, TETB

sgnum = 13
brs = calc_bandreps(sgnum, Val(2))
coefs = [1, 1, 0, 0, 0, 1]
cbr = CompositeBandRep(coefs, brs)

## debug for the BUG #1

br1 = brs[1]
br2 = brs[2]
br3 = brs[6]

## because the problem can be seen in block [1,1], I will study here the term involving 
## br1 → br1.

Rs = [[0, 0]]

δss = TETB.obtain_symmetry_related_hoppings(Rs, br1, br1)

println(δss[RVec([0, 0])][1])

# this results in a vector v bigger than expected probably the sort of error.
