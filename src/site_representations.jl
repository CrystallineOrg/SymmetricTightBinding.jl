
"""
    sgrep_induced_by_siteir_excl_phase(br::NewBandRep, op::SymOperation) -> Matrix{ComplexF64}
    sgrep_induced_by_siteir_excl_phase(cbr::CompositeBandRep, op::SymOperation) -> Matrix{ComplexF64}

Computes the representation matrix of a symmetry operation `op` induced by the site
symmetry group of a band representation `br` or composite band representation `cbr`,
excluding the global momentum-dependent phase factor.
"""
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
            idx_h = findfirst(h′ -> isapprox(h, h′, nothing, false), siteg)
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
        check || error(lazy"failed to find any nonzero block (br=$br, siteg=$siteg)")
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

# ---------------------------------------------------------------------------------------- #
# Site-induced symmetry representation matrix _with_ phase factors

"""
    SiteInducedSGRepElement{D}(ρ::AbstractMatrix, positions::Vector{DirectPoint{D}}, op::SymOperation{D})

Store the symmetry representation matrix `ρ` of a symmetry operation `op` induced by the site
symmetry group of a band representation, along with the positions of the orbitals used for
building the representation. The matrix `ρ` is a square matrix of size `length(positions)`.
"""
struct SiteInducedSGRepElement{D}
    ρ::Matrix{ComplexF64}
    positions::Vector{DirectPoint{D}}
    op::SymOperation{D}
    function SiteInducedSGRepElement{D}(
        ρ::AbstractMatrix,
        positions::Vector{DirectPoint{D}},
        op::SymOperation{D},
    ) where D
        @boundscheck N = LinearAlgebra.checksquare(ρ)
        length(positions) == N || error("length of positions must match the size of ρ")
        new{D}(Matrix{ComplexF64}(ρ), positions, op)
    end
end

# functor behavior for `SiteInducedSGRepElement`
function (sgrep::SiteInducedSGRepElement{D})(k::AbstractVector{<:Real}) where {D}
    g = sgrep.op
    gk = compose(g, ReciprocalPoint{D}(k))
    ρ′ = similar(sgrep.ρ)
    for (β, qᵦ) in enumerate(sgrep.positions)
        for (α, qₐ) in enumerate(sgrep.positions)
            tᵦₐ = compose(g, qₐ) - qᵦ
            e = cispi(-2dot(gk, tᵦₐ)) # e^{-i(gk)·tᵦₐ}
            ρ′[β, α] = sgrep.ρ[β, α] * e # ρᵦₐ(g) e^{-i(gk)·tᵦₐ}
        end
    end

    return ρ′
end

"""
    sgrep_induced_by_siteir(br::NewBandRep, op::SymOperation, positions::Vector{DirectPoint{D}}) -> SiteInducedSGRepElement{D}
    sgrep_induced_by_siteir(cbr::CompositeBandRep, op::SymOperation, positions::Vector{DirectPoint{D}}) -> SiteInducedSGRepElement{D}
    sgrep_induced_by_siteir(tbm::TightBindingModel, op::SymOperation) -> SiteInducedSGRepElement{D}
    sgrep_induced_by_siteir(tbm::ParameterizedTightBindingModel, op::SymOperation) -> SiteInducedSGRepElement{D}

Computes the symmetry representation matrix of a symmetry operation `op` induced by the site
symmetry group of a band representation `br` or `cbr`, including the phase factors that depend
on momentum 𝗸. The matrix `ρ` is a square matrix of size `length(positions)`.
"""
function sgrep_induced_by_siteir(
    br::Union{NewBandRep{D}, CompositeBandRep{D}},
    op::SymOperation{D},
    positions::Vector{DirectPoint{D}} = orbital_positions(br),
) where D
    ρ = sgrep_induced_by_siteir_excl_phase(br, op)
    size(ρ, 1) == length(positions) || error("incompatible dimensions of `ρ` & `positions`")

    return SiteInducedSGRepElement{D}(ρ, positions, op)
end
function sgrep_induced_by_siteir(
    tbm::Union{TightBindingModel{D}, ParameterizedTightBindingModel{D}},
    op::SymOperation{D},
) where D
    return sgrep_induced_by_siteir(CompositeBandRep(tbm), op, orbital_positions(tbm))
end
