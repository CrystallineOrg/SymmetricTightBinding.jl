#=
** Implementation and theory notes **

TRS can be understood as a spatial symmetry when acting on the Hamiltonian:
   D(𝒯)H(k)D(𝒯)† = H(𝒯k) → Γ(𝒯)H*(k)Γ(𝒯)† = H(-k),
where D is the whole operator and Γ is only its unitary part, so D(𝒯) = Γ(𝒯)𝒯.
For *spinless* band representations, 𝒯² = +1, the basis of the representation is real and
Γ(𝒯) can be chosen as Γ(𝒯) = I so that H*(k) = H(-k).
For *spinful* band representations, 𝒯² = -1 and we can (and do, cf. Crystalline's
`timereversal_unitary`) choose Γ(𝒯) = 𝟙 ⊗ J with J = iσʸ ⊗ 𝟙ₙ per site (𝟙 over the
sites of the Wyckoff position's orbit, J over the orbitals at each site: σʸ over the
Kramers index, 𝟙ₙ over the remaining n). In both cases, Γ is real (a signed permutation).

NB: The "realification" of the representation, i.e., the basis choice that brings Γ to
the above forms, is performed in Crystalline via `physical_realify`, applied when
`timereversal = true`. This is done in `bandreps` itself, i.e., can be assumed here.

We have that Hₛₜ(k) = vᵀ(k) Mₛₜ t, with s and t indexing the orbitals of `brₐ` and `brᵦ`,
and Mₛₜ = `Mm[:,:,s,t]` (cf. `construct_M_matrix`). Remember that we split t into real and
imaginary parts: t = [real(t); i imag(t)], and Mm acts as [Mm Mm]. Thus, we can rewrite the
Hamiltonian as:
   Hₛₜ(k) = vᵀ(k) [Mₛₜ Mₛₜ] [real(t); i imag(t)]

Thus, the action of TRS on the Hamiltonian, H(-k) = ΓH*(k)Γ†, can be decomposed into parts:

1. Hₛₜ(-k) = vᵀ(-k) [Mₛₜ Mₛₜ] [real(t); i*imag(t)]

2. (ΓH*(k)Γ†)ₛₜ = Σₘₙ Γₐ,ₛₘ H*ₘₙ(k) Γᵦ,ₜₙ                  (NB: Γ is real, so Γ† = Γᵀ)
                = Σₘₙ Γₐ,ₛₘ (v*)ᵀ(k) [Mₘₙ Mₘₙ] [real(t); -i*imag(t)] Γᵦ,ₜₙ
                = (v*)ᵀ(k) (Σₘₙ Γₐ,ₛₘ [Mₘₙ Mₘₙ] Γᵦ,ₜₙ) [real(t); -i*imag(t)]
                = (v*)ᵀ(k) [M̃ₛₜ -M̃ₛₜ] [real(t); i*imag(t)]
                = vᵀ(-k) [M̃ₛₜ -M̃ₛₜ] [real(t); i*imag(t)]
   with M̃ₛₜ = Σₘₙ Γₐ,ₛₘ Mₘₙ Γᵦ,ₜₙ, using that Γ is real and that v*(k) = v(-k).
   In practice, we compute M̃ slice by slice over the hopping-vector (i) and free-parameter
   (j) axes of Mm, i.e., as the matrix product M̃[i,j,:,:] = Γₐ Mm[i,j,:,:] Γᵦᵀ.

Imposing the condition H(-k) = ΓH*(k)Γ† we get:
  vᵀ(-k) [Mₛₜ-M̃ₛₜ Mₛₜ+M̃ₛₜ] [real(t); i*imag(t)] = 0   (for all s, t)
For spinless band representations, M̃ = M, and this is simply [0 2Mₛₜ] (we use
[0 Mₛₜ], which has the same null space). For spinful band representations, M̃ permutes
the orbitals of each site among Kramers partners (with signs), so TRS relates hoppings,
e.g., t↑↑ = t↓↓* and t↑↓ = -t↓↑*. Note that Γ enters twice, so its overall sign is
irrelevant.

This way of casting the problem is very convenient (and doesn't require any modification of 
the v vectors to include "reversed" hoppings, which is not in general a necessary feature).
We just need to ensure that we intersect the [real(t); i*imag(t)] basis with the null space
of [Mₛₜ-M̃ₛₜ Mₛₜ+M̃ₛₜ] for all s, t.
=#

"""
    obtain_basis_free_parameters_TRS(
        h_orbit::HoppingOrbit{D}, 
        brₐ::BandRep{D}, 
        brᵦ::BandRep{D}, 
        orderingₐ::OrbitalOrdering{D} = OrbitalOrdering(brₐ),
        orderingᵦ::OrbitalOrdering{D} = OrbitalOrdering(brᵦ),
        Mm::AbstractArray{4, Int} = construct_M_matrix(h_orbit, brₐ, brᵦ, orderingₐ, orderingᵦ)
    )                             --> Tuple{Array{Int,4}, Vector{Vector{ComplexF64}}}

Obtain the basis of free parameters for the hopping terms between `brₐ` and `brᵦ` associated
with the hopping orbit `h_orbit` under time-reversal symmetry.

Real and imaginary parts of the basis vectors are differentiated explicitly.
"""
function obtain_basis_free_parameters_TRS(
    h_orbit::HoppingOrbit{D},
    brₐ::BandRep{D},
    brᵦ::BandRep{D},
    orderingₐ::OrbitalOrdering{D} = OrbitalOrdering(brₐ),
    orderingᵦ::OrbitalOrdering{D} = OrbitalOrdering(brᵦ),
    Mm::AbstractArray{Int, 4} = construct_M_matrix(h_orbit, brₐ, brᵦ, orderingₐ, orderingᵦ),
) where {D}
    S = isspinful(brₐ)
    S == isspinful(brᵦ) || error("both band representations must have the same spin")

    # unitary parts `Γ` of time reversal across the orbitals of `brₐ` and `brᵦ`, in the
    # orbital ordering of `OrbitalOrdering` (site-major, partner function-minor)
    Γₐ = site_induced_timereversal_unitary(brₐ)
    Γᵦ = site_induced_timereversal_unitary(brᵦ)

    # NB: we want to keep `_aggregate_constraints` due to its efficiency in building the
    # constraint matrix. So, although seemingly unnecessary, we stick with its Q & Z tensor
    # structure, implementing the constraint tensor as Q, with Z = 0, with the final
    # constraints being a row-wise aggregation of Q-Z

    # Step 1: the Z tensor (zero-valued tensor)
    Z = 0 # stand-in for a zero-tensor (but no need to allocate it explicitly)

    # Step 2: the Q tensor (the block tensor "[Mᵢⱼ-M̃ᵢⱼ Mᵢⱼ+M̃ᵢⱼ]")
    sMm = size(Mm)
    Q = zeros(Int, (sMm[1], 2sMm[2], sMm[3], sMm[4]))
    if !S # spinless: `M̃ = M`, so the constraint tensor is `[0 M]`
        dst_indices = CartesianIndices((1:sMm[1], (sMm[2]+1):2sMm[2], 1:sMm[3], 1:sMm[4]))
        copyto!(Q, dst_indices, Mm, CartesianIndices(axes(Mm))) # efficient [0 Mm]
    else
        # spinful: the constraint tensor is `[M-M̃ M+M̃]`, with `M̃ = Γₐ M Γᵦᵀ` acting on
        # the orbital axes; since `Γₐ` and `Γᵦ` are signed permutations, `M̃` is integer
        for j in axes(Mm, 2), i in axes(Mm, 1)
            Mᵢⱼ = @view Mm[i, j, :, :]
            M̃ᵢⱼ = Γₐ * Mᵢⱼ * Γᵦ'
            Q[i, j, :, :]        .= Mᵢⱼ .- M̃ᵢⱼ # first block "M-M̃"
            Q[i, sMm[2]+j, :, :] .= Mᵢⱼ .+ M̃ᵢⱼ # second block "M+M̃"
        end
    end

    # Step 3: use `_aggregate_constraints` to build the constraint matrix 
    #         ~(Q-Z)[i, :, s, t] (aggregated over i,s,t)
    constraints = _aggregate_constraints(Q, Z)
    tₐᵦ_basis_matrix = nullspace(constraints; atol = NULLSPACE_ATOL_DEFAULT)

    return tₐᵦ_basis_matrix
end

"""
    site_induced_timereversal_unitary(br::BandRep)           --> AbstractMatrix{<:Real}
    site_induced_timereversal_unitary(cbr::CompositeBandRep) --> Matrix{<:Real}

Return the unitary part `Γ` of time reversal `𝒯 = ΓK` across all orbitals of `br`, i.e.,
across sites in the orbit of the Wyckoff position and site-symmetry orbitals at each site.
This is the time-reversal counterpart of [`sgrep_induced_by_siteir`](@ref): i.e., the
(unitary part of the) action of time-reversal symmetry on the orbitals of a band
representation.

The orbitals are ordered according to `OrbitalOrdering(br)`: i.e., the returned matrix is
`𝟙 ⊗ Γₛ`, with `Γₛ = Crystalline.timereversal_unitary(br.siteir)` the unitary part for a
single site (`𝟙` for spinless, `J = iσʸ ⊗ 𝟙ₙ` for spinful site irreps).

If a `CompositeBandRep` is supplied as input, the matrix is the block-diagonal
generalization, stacking over the contained band representation content.
"""
function site_induced_timereversal_unitary(br::BandRep)
    Γₛ = timereversal_unitary(br.siteir)
    n = length(cosets(group(br))) # number of sites in the orbit of `br`'s Wyckoff position
    return kron(I(n), Γₛ)
end

function site_induced_timereversal_unitary(cbr::CompositeBandRep)
    T = isspinful(cbr) ? Int : Bool # return Matrix{Int} (spinful) / Matrix{Bool} (spinless)
    return _apply_across_matrix_blocks_of_composite_bandrep(
        site_induced_timereversal_unitary, cbr, T
    )
end
