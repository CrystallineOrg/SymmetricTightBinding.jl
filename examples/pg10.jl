using Pkg
Pkg.activate(@__DIR__)

using Crystalline, TETB

##- Compute the necessary things for obtaining the hoppings

pg_num = 10
brs = calc_bandreps(pg_num, Val(2))

# needs to do that to find the WPs properly
br1 = brs[1]
br2 = brs[5]

wp1 = orbit(group(br1))
wp2 = orbit(group(br2))

sgrep1 = sgrep_induced_by_siteir_generators(br1)
sgrep2 = sgrep_induced_by_siteir_generators(br2)

##- Compute the orbits of Δ's taking into considerations the symmetries ------------------##

Rs = [[0, 0]] # vector containing the translations we want to consider

δss = TETB.obtain_symmetry_related_hoppings(Rs, br1, br2)

##- Compute the matrix M that will encode the Hamiltonian as a numerical matrix ----------##

or, Mv = TETB.construct_M_matrix(first(values(δss)), br1, br2)

##----------------------------------------------------------------------------------------##
