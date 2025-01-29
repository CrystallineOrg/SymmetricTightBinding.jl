using Pkg
Pkg.activate(@__DIR__)

# this example is for a _centered_ lattice; the key is to show that the orbits don't
# erroneously go outside the primitive cell (since the primitive and conventional cells)
# don't coincide

using Crystalline, TETB

##- Compute the necessary things for obtaining the hoppings

pg_num = 5
brs = calc_bandreps(pg_num, Val(2))

# needs to do that to find the WPs properly
br = brs[1]

ops = spacegroup(num(br), dim(br))
wps = orbit(group(br))

gens, sgrep = sgrep_induced_by_siteir_generators(br)

##- Compute the orbits of Δ's taking into considerations the symmetries ------------------##

Rs = [[0,0], [1, 0], [1,-1]] # vector containing the translations we want to consider

Δs = TETB.find_symmetry_related_hoppings(Rs, br, br)

##----------------------------------------------------------------------------------------##