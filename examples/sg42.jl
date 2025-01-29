using Pkg
Pkg.activate(@__DIR__)

using Crystalline, TETB

sg_num = 42
brs = calc_bandreps(sg_num)
coefs = [1, 0, 1, 0, 0, 0]
cbr = CompositeBandRep(coefs, brs)

## debug for the BUG #3

br1 = brs[1]
br2 = brs[3]

Rs = [[0, 0, 0]]

δss = TETB.obtain_symmetry_related_hoppings(Rs, br1, br2)

wps1 = orbit(group(br1))
wps2 = orbit(group(br2))

cntr = centering(sg_num)
wps1 = primitivize.(wps1, cntr)
wps2 = primitivize.(wps2, cntr)

println("wps1 = ", wps1)
println("wps2 = ", wps2)

println("δs = ", last(δss[RVec([-1 / 2, 1 / 2, -1])][1]))

# Possible problem on primitivize and/or reduce_translation_to_unitrange