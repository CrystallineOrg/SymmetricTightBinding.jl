
"""
    site_induced_sgrep_excl_phase(
        br::BandRep, op::AbstractOperation, [ordering::OrbitalOrdering]
    )
    site_induced_sgrep_excl_phase(cbr::CompositeBandRep, op::AbstractOperation)
                                                                 --> Matrix{ComplexF64}

Return the representation matrix of a symmetry operation `op` induced by the site
symmetry group of a band representation `br` or composite band representation `cbr`,
excluding the global momentum-dependent phase factor.

The rows and columns of the matrix are ordered according to `ordering` (by default,
`OrbitalOrdering(br)`); for `cbr`, according to [`OrbitalOrdering`](@ref) of each band
representation.

For a spinful (double-valued) band representation, i.e., when `isspinful(br)` is true, `op`
must be a `DSymOperation`, since the representation depends on the SU(2) element in addition
to the spatial operation. For spinless (single-valued) band representations, `op` must be a
`SymOperation`.

# Note
This function assumes Convention 1 for the Fourier transform, so the momentum dependence is
introduced as a global phase factor. This is not true if Convention 2 is used. See 
`/docs/src/theory.md` for more details.
"""
function site_induced_sgrep_excl_phase(
    br::BandRep{D},
    op::AbstractOperation{D},
    ordering::OrbitalOrdering{D} = OrbitalOrdering(br),
) where {D}
    _check_operation_spin(br, op)
    # NB: `bandreps` in Crystalline already applies `physical_realify` if
    #     `timereversal` is true, so we don't need to manually redo it for `siteir` below
    siteir = br.siteir
    siteg = primitivize(group(siteir))
    wps = orbit(siteg)
    g = op

    # orbital (row/column) index of the `j`th partner function at the `α`th site of the
    # orbit is `idxs[α, j]`, as given by `ordering` (whose `site_idx`s index the same orbit
    # as `wps`)
    idxs = Matrix{Int}(undef, length(wps), irdim(siteir))
    for (n, o) in enumerate(ordering)
        idxs[o.site_idx, o.partner_idx] = n
    end

    ρ = zeros(ComplexF64, length(ordering), length(ordering))
    for (α, (gₐ, qₐ)) in enumerate(zip(cosets(siteg), wps))
        check = false
        for (β, (gᵦ, qᵦ)) in enumerate(zip(cosets(siteg), wps))
            tᵦₐ = constant(g * parent(qₐ) - parent(qᵦ)) # ignore free parts of the WP
            # compute h = gᵦ⁻¹ tᵦₐ⁻¹ g gₐ
            h = compose(
                compose(compose(inv(gᵦ), typeof(g)(-tᵦₐ), false), g, false),
                gₐ,
                false,
            )
            idx_h = findfirst(h′ -> isapprox(h, h′, nothing, false), siteg)
            if !isnothing(idx_h) # h ∈ siteg and qₐ and qᵦ are connected by `g`
                # assign the "block" (β, α) (rows/columns of orbitals at qᵦ/qₐ); we go
                # "through" `idxs` because the partner-function-major ordering doesn't
                # actually produce (β, α)-site blocks (but (i, j)-partner-function blocks)
                @views ρ[idxs[β, :], idxs[α, :]] .= siteir.matrices[idx_h]

                # We build the representation acting as the transpose, i.e.,
                # gΦ(k) = ρᵀ(g)Φ(Rk), where Φ(k) is the site-symmetry function of the
                # bandrep. This yields ρⱼᵦᵢₐ(g) = e(-i(gk)·v) Γⱼᵢ(g) δ(gqₐ, qᵦ), where Γ is
                # the representation of the site-symmetry group, and δ(gqₐ, qᵦ) is the
                # Kronecker delta mod τ ∈ T. We are building it as the transpose because we
                # want to ensure the proper composition order: ρ(g₁g₂) = ρ(g₁)ρ(g₂).
                # NB: We do not include the (usually redundant) exponential (k-dependent)
                #     phases. Note that these phases are NOT REDUNDANT if we mean to
                #     use the sgrep as the group action on eigenstates, e.g., for
                #     determining the irreps of a tight-binding Hamiltonian; for this, use
                #     `site_induced_sgrep` instead.
                check = true
                break
            end
        end
        check || error(lazy"failed to find any nonzero block (br=$br, siteg=$siteg, op=$op)")
    end

    return ρ
end

function site_induced_sgrep_excl_phase(
    cbr::CompositeBandRep{D},
    op::AbstractOperation{D},
) where {D}
    f = Base.Fix2(site_induced_sgrep_excl_phase, op) # = br -> sgrep_…(br, op)
    return _apply_across_matrix_blocks_of_composite_bandrep(f, cbr, ComplexF64)
end

function _apply_across_matrix_blocks_of_composite_bandrep(
    f::F,
    cbr::CompositeBandRep,
    ::Type{T} # eltype
) where {F, T}
    N = occupation(cbr)
    ρ = zeros(T, N, N) # TODO: `spzeros` would be nice here (but then requires SparseArrays)
    j = 0
    for (cᵢ, brᵢ) in zip(cbr.coefs, cbr.brs)
        iszero(cᵢ) && continue
        Nᵢ = occupation(brᵢ)
        ρᵢ = f(brᵢ)
        for _ in 1:Int(cᵢ) # slot in the block `ρᵢ` into `ρ` for each copy of the BR
            ρ[j+1:j+Nᵢ, j+1:j+Nᵢ] .= ρᵢ
            j += Nᵢ
        end
    end
    return ρ
end

# check that `x` (e.g., a `BandRep`, `TightBindingModel`, `AbstractIrrep`, or `AbstractGroup`)
# and the operation(s) `op` agree on `isspinful`; error otherwise
@inline function _check_operation_spin(x, op::AbstractOperation)
    isspinful(x) == isspinful(op) && return nothing
    if isspinful(x)
        error(lazy"a spinful setting requires a double group operation (`DSymOperation`): got a `$(typeof(op))`")
    else
        error(lazy"a spinless setting requires a spinless operation (`SymOperation`): got a `$(typeof(op))`")
    end
end
_check_operation_spin(x, ops::AbstractVector{<:AbstractOperation}) = _check_operation_spin(x, first(ops))

# ---------------------------------------------------------------------------------------- #
# Site-induced symmetry representation matrix _with_ phase factors

"""
    SiteInducedSGRepElement{D}(
        ρ::AbstractMatrix,
        positions::Vector{DirectPoint{D}},
        op::AbstractOperation{D}
    )

Represents a matrix-valued element of a site-induced representation of a space group,
including a global momentum-dependent phase factor.

This structure behaves like a functor: calling it with a momentum `k :: AbstractVector` 
returns the matrix representation at `k`.

## Fields (internal)
- `ρ :: Matrix{ComplexF64}` : The momentum-independent matrix part of the representation.
- `positions :: Vector{DirectPoint{D}}`: Real-space positions corresponding to the orbitals
  in the orbit of the associated site-symmetry group.
"""
struct SiteInducedSGRepElement{D, O<:AbstractOperation{D}}
    ρ::Matrix{ComplexF64}
    positions::Vector{DirectPoint{D}}
    op::O
    function SiteInducedSGRepElement{D}(
        ρ::AbstractMatrix,
        positions::Vector{DirectPoint{D}},
        op::O,
    ) where {D, O<:AbstractOperation{D}}
        @boundscheck N = LinearAlgebra.checksquare(ρ)
        length(positions) == N || error("length of positions must match the size of ρ")
        new{D, O}(Matrix{ComplexF64}(ρ), positions, op)
    end
end

# functor behavior for `SiteInducedSGRepElement`
function (sgrep::SiteInducedSGRepElement{D})(k::AbstractVector{<:Real}) where {D}
    g = sgrep.op
    v = translation(g)
    gk = compose(g, ReciprocalPoint{D}(k))
    Dₖ = cispi(-2dot(gk, v)) * sgrep.ρ # e^{-i(gk)·v} ρ(g)

    return Dₖ
end

"""
    site_induced_sgrep(
        br::Union{BandRep, CompositeBandRep},
        op::AbstractOperation, [positions::Vector{<:DirectPoint}]
    )
    site_induced_sgrep(
        tbm::Union{TightBindingModel,ParameterizedTightBindingModel}, op::AbstractOperation
    )
        --> SiteInducedSGRepElement

Computes the representation matrix of a symmetry operation `op` induced by the site symmetry
group associated with an elementary or composite band representation `br` , _including_ the global
momentum-dependent phase factor, returning a `SiteInducedSGRepElement`, which is a functor
over momentum inputs.

A (possibly parameterized) tight-binding model `tbm` can be specified instead of a band representation,
in which case the latter is inferred from the former.

For spinful band representations, `op` must be a `DSymOperation`, i.e., carry its SU(2)
element.
"""
function site_induced_sgrep(
    br::Union{BandRep{D}, CompositeBandRep{D}},
    op::AbstractOperation{D},
    positions::Vector{DirectPoint{D}} = orbital_positions(br),
) where D
    ρ = site_induced_sgrep_excl_phase(br, op)
    size(ρ, 1) == length(positions) || error("incompatible dimensions of `ρ` & `positions`")

    return SiteInducedSGRepElement{D}(ρ, positions, op)
end
function site_induced_sgrep(
    tbm::Union{TightBindingModel{D}, ParameterizedTightBindingModel{D}},
    op::AbstractOperation{D},
) where D
    return site_induced_sgrep(CompositeBandRep(tbm), op, orbital_positions(tbm))
end
