"""
Obtain a bandrep decomposition for the symmetry vector of the bands provided `m` with a minimal
number of auxiliary bands in the interval `[μᴸ_min,μᴸ_max]`.

If the photonic bands are connected to zero frequency corrections to the singularity at Γ are
made. This parameter is settle by default to `true`
"""
function find_bandrep_decompositions(
    m::AbstractSymmetryVector{D},
    brs::Collection{NewBandRep{D}};
    μᴸ_min::Integer=0,
    μᴸ_max::Integer=μᴸ_min + 2 * occupation(m),
    connected_to_zero_frequency::Bool=true
) where {D}

    if D < 3 || !connected_to_zero_frequency
        μ = occupation(m)
        μs_brs = occupation.(brs)
        idxs_k = eachindex(first(brs))

        # we don't need any longitudinal modes so we can directly find the expansions using
        # `m` as a constraint
        idxs_sol = find_all_admissible_expansions(brs, μs_brs, μ, Vector(m), idxs_k)

        longitudinal = Crystalline.CompositeBandRep_from_indices(Int[], brs)
        apolar = Crystalline.CompositeBandRep_from_indices.(idxs_sol, Ref(brs))

        !isempty(apolar) || error("Check the symmetry vector and space group used")

        return TightBindingCandidateSet(longitudinal, apolar, [Float64[]])

    else

        μᴸ = μᴸ_min - 1
        while μᴸ < μᴸ_max
            μᴸ += 1
            idxsᴸs = find_auxiliary_modes(μᴸ, brs)
            (isempty(idxsᴸs) && μᴸ ≠ 0) && continue
            # compute all possible decomposition of m into valid combinations of nᴸ and nᵀ⁺ᴸ
            candidatesv = find_apolar_modes(m, idxsᴸs, brs)
            isempty(candidatesv) || return candidatesv
        end
        error("""failed to find possible auxiliary-apolar decompositions for provided \
                symmetry vector in search range for auxiliary modes; increasing kwarg \
                `μᴸ_max` may help, if a decomposition exists""")
    end
end

function sgrep_induced_by_siteir(br::NewBandRep{D}, op::SymOperation{D}) where {D}
    # NB: `calc_bandreps` in Crystalline already applies `physical_realify` if
    #     `timereversal` is true, so we don't need to manually redo it for `siteir` below
    siteir = br.siteir
    siteir_dim = irdim(siteir)
    siteg = group(siteir)
    wps = orbit(siteg)
    mult = length(wps)
    g = op

    ρ = BlockArray{ComplexF64}(zeros(ComplexF64, siteir_dim * mult, siteir_dim * mult),
        fill(siteir_dim, mult), fill(siteir_dim, mult))

    for (α, (gₐ, qₐ)) in enumerate(zip(cosets(siteg), wps))
        check = false
        for (β, (gᵦ, qᵦ)) in enumerate(zip(cosets(siteg), wps))
            tᵦₐ = constant(g * parent(qₐ) - parent(qᵦ)) # ignore free parts of the WP
            # compute h = gᵦ⁻¹ tᵦₐ⁻¹ g gₐ
            h = compose(compose(compose(inv(gᵦ), SymOperation(-tᵦₐ), false), g, false), gₐ, false)
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

                # NB: we do not include the (usually redundant) exponential (k-dependent) phases.
                # Note that these phases are NOT REDUNDANT if we are intending to use the sgrep as 
                # the group action on eigenstates, e.g., for determining the irreps of a tight-binding
                # Hamiltonian
                check = true
                break
            end
        end
        check || error("failed to find any nonzero block")
    end

    return ρ
end

function sgrep_induced_by_siteir(cbr::CompositeBandRep{D}, op::SymOperation{D}) where {D}
    ρ = Matrix{Complex}(undef, 0, 0)
    for (idxc, c) in enumerate(cbr.coefs)
        iszero(c) && continue
        ρ_idxc = sgrep_induced_by_siteir(cbr.brs[idxc], op)
        for _ in 1:Int(c)
            ρ = ρ ⊕ ρ_idxc
        end
    end
    return ρ
end

# NB: we do not include the (usually redundant) exponential (k-dependent) phases below.
#     Note that these phases are NOT REDUNDANT if we are intending to use the sgrep as 
#     the group action on eigenstates, e.g., for determining the irreps of a tight-binding
#     Hamiltonian
"""
    sgrep_induced_by_siteir_generators(
        br::NewBandRep{D},
        gens::AbstractVector{SymOperation{D}})
    --> Tuple{Vector{SymOperation{D}}, Vector{Matrix{ComplexF64}}}

Induce a representation for the generators of the SG from a representation of the site-symmetry 
group of a particular maximal WP.
"""
function sgrep_induced_by_siteir_generators(
    br::NewBandRep{D},
    gens::AbstractVector{SymOperation{D}}=generators(num(br), SpaceGroup{D}),
) where {D}
    return gens, sgrep_induced_by_siteir.(Ref(br), gens)
end

"""
Direct sum of matrices
"""
function ⊕(As::AbstractMatrix...)
    return cat(As..., dims=Val((1, 2)))
end


function sgrep_induced_by_siteir_generators(
    cbr::CompositeBandRep{D},
    gens::AbstractVector{SymOperation{D}}=generators(num(cbr), SpaceGroup{D}),
) where {D}
    # TODO: primitivize the generators of the space group
    return gens, sgrep_induced_by_siteir.(Ref(cbr), gens)
end