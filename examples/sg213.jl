using Pkg
Pkg.activate(@__DIR__)

using Crystalline, TETB

sgnum = 213
brs = calc_bandreps(sgnum)
coefs = zeros(Int, length(brs))
coefs[end] = 1
cbr = CompositeBandRep(coefs, brs)

ordering1 = TETB.OrbitalOrdering(cbr.brs[end])

Rs = [[0, 0, 0]]

tb_model = tb_hamiltonian(cbr, Rs)

# let me check how the basis vectors look like

h_orbits = obtain_symmetry_related_hoppings(Rs, cbr.brs[end], cbr.brs[end])

Mm, tₐᵦ_basis = TETB.obtain_basis_free_parameters(h_orbits[2], cbr.brs[end], cbr.brs[end])

# checks

c = Int[]
for i in eachindex(tₐᵦ_basis)
    if !isreal(tₐᵦ_basis[i])
        append!(c, i)
    end
end
println(c)

# tₐᵦ_basis[c] will give back the basis vectors with complex entries. Let me
# focus on tₐᵦ_basis[2]

println(tₐᵦ_basis[2])