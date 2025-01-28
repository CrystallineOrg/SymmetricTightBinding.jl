using Pkg
Pkg.activate(@__DIR__)

### necessary packages
using Crystalline, TETB

sg_num = 224
brs = calc_bandreps(sg_num, Val(3))
c = zeros(Int, length(brs))
c[13] = 1
c[19] = 1
cbr = CompositeBandRep(c, brs)

tb_model = TETB.tb_hamiltonian(cbr, [[0, 0, 0]]) # WARNING: we obtain different number of 
# hoppings across orbits. Let me check that.

hops = TETB.obtain_symmetry_related_hoppings([[0, 0, 0]], brs[13], brs[19])