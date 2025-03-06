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

# we do not include the (usually redundant) exponential (k-dependent) phases below
"""
Induce a representation for the generators of the SG from a representation of the site-symmetry 
group of a particular maximal WP.
"""
function sgrep_induced_by_siteir_generators(
    br::NewBandRep{D},
    gens::AbstractVector{SymOperation{D}}=generators(num(br), SpaceGroup{D})
) where {D}
    siteir = br.siteir
    siteir_dim = irdim(siteir)
    siteg = group(siteir)
    wps = orbit(siteg)
    mult = length(wps)


    ρs = [BlockArray{ComplexF64}(
        zeros(ComplexF64, siteir_dim * mult, siteir_dim * mult),
        fill(siteir_dim, mult), fill(siteir_dim, mult)) for _ in eachindex(gens)]
    for (n, g) in enumerate(gens)
        ρ = ρs[n]
        for (α, (gₐ, qₐ)) in enumerate(zip(cosets(siteg), wps))
            check = false
            for (β, (gᵦ, qᵦ)) in enumerate(zip(cosets(siteg), wps))
                tᵦₐ = constant(g * parent(qₐ) - parent(qᵦ)) # ignore free parts of the WP
                # compute h = gᵦ⁻¹ tᵦₐ⁻¹ g gₐ
                h = compose(compose(compose(inv(gᵦ), SymOperation(-tᵦₐ), false), g, false), gₐ, false)
                idx_h = findfirst(==(h), siteg)
                if !isnothing(idx_h) # h ∈ siteg and qₐ and qᵦ are connected by g
                    ρ[Block(β, α)] .= siteir.matrices[idx_h]
                    check = true
                    break
                end
            end
            check || error("failed to find any nonzero block")
        end
    end

    return gens, ρs
end

"""
Direct sum of matrices
"""
function ⊕(As::AbstractMatrix...)
    return cat(As..., dims=Val((1, 2)))
end


function sgrep_induced_by_siteir_generators(brs::CompositeBandRep{D}) where {D}
    # TODO: primitivize the generators of the space group
    gens = generators(num(brs), SpaceGroup{D})
    ρs = [Matrix{Complex}(undef, 0, 0) for _ in eachindex(gens)]

    for (idxc, c) in enumerate(brs.coefs)
        sgrep = sgrep_induced_by_siteir_generators(brs.brs[idxc], gens)
        for idxg in eachindex(gens)
            iszero(c) && continue
            ρ_idxg = sgrep[2][idxg]
            for _ in 1:Int(c)
                ρs[idxg] = ρs[idxg] ⊕ ρ_idxg
            end
        end
    end

    return gens, ρs
end