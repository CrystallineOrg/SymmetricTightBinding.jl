function sgrep_induced_by_siteir_excl_phase(
    br::NewBandRep{D},
    op::SymOperation{D},
) where {D}
    # NB: `calc_bandreps` in Crystalline already applies `physical_realify` if
    #     `timereversal` is true, so we don't need to manually redo it for `siteir` below
    siteir = br.siteir
    siteir_dim = irdim(siteir)
    siteg = group(siteir)
    wps = orbit(siteg)
    mult = length(wps)
    g = op

    block_axis = BlockedOneTo(collect(siteir_dim:siteir_dim:mult*siteir_dim))
    ρ = zeros(ComplexF64, block_axis, block_axis) # `BlockedMatrix` backed by a `Matrix`
    for (α, (gₐ, qₐ)) in enumerate(zip(cosets(siteg), wps))
        check = false
        for (β, (gᵦ, qᵦ)) in enumerate(zip(cosets(siteg), wps))
            tᵦₐ = constant(g * parent(qₐ) - parent(qᵦ)) # ignore free parts of the WP
            # compute h = gᵦ⁻¹ tᵦₐ⁻¹ g gₐ
            h = compose(
                compose(compose(inv(gᵦ), SymOperation(-tᵦₐ), false), g, false),
                gₐ,
                false,
            )
            idx_h = findfirst(==(h), siteg)
            if !isnothing(idx_h) # h ∈ siteg and qₐ and qᵦ are connected by `g`
                ρ[Block(β, α)] .= siteir.matrices[idx_h]
                # we build the representation acting as the transpose, i.e., 
                # gΦ(k) = ρᵀ(g)Φ(Rk), where Φ(k) is the site-symmetry function of
                # the bandrep. This yields ρⱼᵦᵢₐ(g) = e(-i(Rk)·v) Γⱼᵢ(g) δ(gqₐ, qᵦ),
                # where Γ is the representation of the site-symmetry group, and 
                # δ(gqₐ, qᵦ) is the Kronecker delta mod τ ∈ T.
                # we are building it as the transpose because we want to keep good
                # composition order: ρ(g₁g₂) = ρ(g₁)ρ(g₂). Check `trs_notes.md`.

                # NB: We do not include the (usually redundant) exponential (k-dependent) 
                #     phases. Note that these phases are NOT REDUNDANT if we mean to use
                #     use the sgrep as the group action on eigenstates, e.g., for
                #     determining the irreps of a tight-binding Hamiltonian; for this, use
                #     `sgrep_induced_by_siteir` instead.
                check = true
                break
            end
        end
        check || error("failed to find any nonzero block")
    end

    return ρ
end

function sgrep_induced_by_siteir_excl_phase(
    cbr::CompositeBandRep{D},
    op::SymOperation{D},
) where {D}
    N = occupation(cbr)
    ρ = zeros(Complex, N, N)
    j = 0
    for (cᵢ, brᵢ) in zip(cbr.coefs, cbr.brs)
        iszero(cᵢ) && continue
        Nᵢ = occupation(brᵢ)
        ρᵢ = sgrep_induced_by_siteir_excl_phase(brᵢ, op)
        for _ in 1:Int(cᵢ)
            ρ[j+1:j+Nᵢ, j+1:j+Nᵢ] .= ρᵢ
            j += Nᵢ
        end
    end
    return ρ
end

# NB: we do not include the (usually redundant) exponential (k-dependent) phases below.
#     Note that these phases are NOT REDUNDANT if we are intending to use the sgrep as 
#     the group action on eigenstates, e.g., for determining the irreps of a tight-binding
#     Hamiltonian
"""
    sgrep_induced_by_siteir_generators(br::NewBandRep{D}) where {D}
    --> Tuple{Vector{SymOperation{D}}, Vector{BlockMatrix{ComplexF64}}}

Induce a representation for the generators of the SG from a representation of the site-symmetry 
group of a particular maximal WP.
"""
function sgrep_induced_by_siteir_generators(br::NewBandRep{D}) where {D}
    gens = generators(num(br), SpaceGroup{D})
    return gens, sgrep_induced_by_siteir.(Ref(br), gens)
end

"""
Direct sum of matrices
"""
function ⊕(As::AbstractMatrix...)
    return cat(As...; dims = Val((1, 2)))
end

function sgrep_induced_by_siteir_generators(cbr::CompositeBandRep{D}) where {D}
    # TODO: primitivize the generators of the space group
    gens = generators(num(cbr), SpaceGroup{D})
    return gens, sgrep_induced_by_siteir.(Ref(cbr), gens)
end
