using Pkg
Pkg.activate(@__DIR__)

using Crystalline, TETB

##- Compute the necessary things for obtaining the hoppings

pg_num = 10
brs = calc_bandreps(pg_num, Val(2))

# needs to do that to find the WPs properly
br = brs[end-3]

ops = spacegroup(num(br), dim(br))
wps = orbit(group(br))

sgrep = sgrep_induced_by_siteir_generators(br)

##- Compute the orbits of Δ's taking into considerations the symmetries ------------------##

Rs = [[0,0], [1, 0]] # vector containing the translations we want to consider

Δs = TETB.find_symmetry_related_hoppings(Rs, br, br)

println(Δs)

##----------------------------------------------------------------------------------------##