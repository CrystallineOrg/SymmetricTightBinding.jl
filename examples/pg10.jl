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

δs = TETB.obtain_symmetry_related_hoppings(Rs, br1, br2)

##- Compute the matrix M that will enconde the Hamiltonian as a numerical matrix ---------##

# EBRs: (q|A), (w|B)
# Wyckoff positions: q, w
#   q: q1, ..., qN
#   w: w1, ..., wM
# Site symmetry irreps: A, B
#   A: A1, ..., AJ
#   B: B1, ..., BK
# δs = [δ1, δ2, ..., δn]
#   δ1: qi₁¹ -> wj₁¹, qi₁² -> wj₁², ...
#   δ2: qi₂¹ -> wj₂¹, qi₂² -> wj₂², ...
# v = [exp(ik⋅δ1), exp(ik⋅δ2), ..., exp(ik⋅δn)]
# t = [[t(δ1) ...], [t(δ2) ...], ..., [t(δn) ...]]
#   t(δ1): [t(qi₁ᵅ -> wj₁ᵅ, A_f -> B_g) ...]

# Current example: (1a|E), (2c|A)
#   ___w2__
#  |   x   |
#  |q1 x   x w1
#  |_______|
#   δs = [1/2x, -1/2x, 1/2y, -1/2y]
#      δ1: q1 -> w1 + G1
#      δ2: q1 -> w1 + G2
#      δ3: q1 -> w2 + G3
#      δ4: q1 -> w2 + G4
# t = [t(δ1)..., t(δ2)..., t(δ3)..., t(δ4)...]
#   t(δ1): [t(q1 -> w1, G1, E1 -> A1), t(q1 -> w1, G1, E2 -> A1)]
#   t(δ2): [t(q1 -> w1, G2, E1 -> A1), t(q1 -> w1, G2, E2 -> A1)]
#   t(δ3): [t(q1 -> w2, G3, E1 -> A1), t(q1 -> w2, G3, E2 -> A1)]
#   t(δ4): [t(q1 -> w2, G4, E1 -> A1), t(q1 -> w2, G4, E2 -> A1)]

function construct_M_matrix(
    δs::Dict{RVec{D},Dict{RVec{D},Vector{Tuple{RVec{D},RVec{D},RVec{D}}}}},
    br1::NewBandRep{D},
    br2::NewBandRep{D}
) where {D}

    Ms = Matrix[]
    for Δ in keys(δs)
        δ = δs[Δ]
        M = Matrix{Float64}(undef, 0, 0)

        ## Introduce code that assings a 1 to the correct position

        ##

        push!(Ms, M)
    end
    return Ms
end

##--------------------------------------------------------------------------------------.##