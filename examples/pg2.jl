using Pkg
Pkg.activate(@__DIR__)

using Crystalline, TETB

##- Compute the necessary things for obtaining the hoppings

pg_num = 2
brs = calc_bandreps(pg_num, Val(2))
c = [0, 0, 0, 0, 1, 0, 0, 1]
cbr = CompositeBandRep(c, brs)

# needs to do that to find the WPs properly
br1 = brs[end]
br2 = brs[end-3]

gen = generators(num(cbr))
wps1 = orbit(group(br1))
wps2 = orbit(group(br2))

gens, sgrep = sgrep_induced_by_siteir_generators(cbr)

##- Compute the orbits of Δ's taking into considerations the symmetries ------------------##

Rs = [[0, 0]] # vector containing the translations we want to consider

δs = TETB.obtain_symmetry_related_hoppings(Rs, br1, br2)

##----------------------------------------------------------------------------------------##