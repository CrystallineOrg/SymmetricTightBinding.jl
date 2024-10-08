using Pkg
Pkg.activate(@__DIR__)

using Crystalline, TETB

##- Compute the necessary things for obtaining the hoppings

pg_num = 10
brs = calc_bandreps(pg_num, Val(2))
c = [0, 0, 0, 0, 1, 0, 0, 0]
cbr = CompositeBandRep(c, brs)

# needs to do that to find the WPs properly
br = brs[end-3]

gen = generators(num(cbr))
wps = orbit(group(br))

sgrep = sgrep_induced_by_siteir_generators(cbr)

##- ⊕