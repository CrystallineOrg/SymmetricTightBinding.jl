# Just a file with the functions that we are building to solve all the things

function prune_klab_irreps_brs!(brs::BandRepSet, klab::String="Γ")
    prune_kidx = findfirst(==(klab), brs.klabs)
    isnothing(prune_kidx) && error(lazy"could not find $klab among included k-points")
    deleteat!(brs.klabs, prune_kidx)

    prune_iridxs = findall(irlab -> klabel(irlab) == klab, brs.irlabs)
    isempty(prune_iridxs) && error(lazy"could not find $klab among included irreps")
    deleteat!(brs.irlabs, prune_iridxs)

    foreach(brs.bandreps) do br
        deleteat!(br.irvec, prune_iridxs)
    end

    return brs
end

function prune_klab_irreps_brs(brs::BandRepSet, klab::String="Γ")
    irlabs′ = copy(brs.irlabs)
    brs′ = BandRepSet(
        brs.sgnum,
        map(brs.bandreps) do br
            BandRep(
                br.wyckpos,
                br.sitesym,
                br.label,
                br.dim,
                br.decomposable,
                br.spinful,
                copy(br.irvec),
                irlabs′
            )
        end,
        copy(brs.kvs),
        copy(brs.klabs),
        irlabs′,
        brs.allpaths,
        brs.spinful,
        brs.timereversal
    )
    return prune_klab_irreps_brs!(brs′, klab)
end

function pick_klab_irreps_brs!(brs::BandRepSet, klab::String="Γ")
    prune_kidx = findall(!=(klab), brs.klabs)
    isnothing(prune_kidx) && error(lazy"could not find $klab among included k-points")
    deleteat!(brs.klabs, prune_kidx)

    prune_iridxs = findall(irlab -> klabel(irlab) != klab, brs.irlabs)
    isempty(prune_iridxs) && error(lazy"could not find $klab among included irreps")
    deleteat!(brs.irlabs, prune_iridxs)

    foreach(brs.bandreps) do br
        deleteat!(br.irvec, prune_iridxs)
    end

    return brs
end

function pick_klab_irreps_brs(brs::BandRepSet, klab::String="Γ")
    irlabs′ = copy(brs.irlabs)
    brs′ = BandRepSet(
        brs.sgnum,
        map(brs.bandreps) do br
            BandRep(
                br.wyckpos,
                br.sitesym,
                br.label,
                br.dim,
                br.decomposable,
                br.spinful,
                copy(br.irvec),
                irlabs′
            )
        end,
        copy(brs.kvs),
        copy(brs.klabs),
        irlabs′,
        brs.allpaths,
        brs.spinful,
        brs.timereversal
    )
    return pick_klab_irreps_brs!(brs′, klab)
end


function prune_klab_irreps_vecs!(v::BandSummary, klab::String="Γ")
    prune_iridxs = findall(irlab -> klabel(irlab) == klab, v.brs.irlabs)
    isempty(prune_iridxs) && error(lazy"could not find $klab among included irreps")
    deleteat!(v.n, prune_iridxs)

    prune_klab_irreps_brs!(v.brs, klab)

    return v
end

function prune_klab_irreps_vecs(v::BandSummary, klab::String="Γ")
    brs = v.brs
    irlabs´ = copy(brs.irlabs)
    v´ = BandSummary(
        v.topology,
        v.bands,
        copy(v.n),
        BandRepSet(
            brs.sgnum,
            map(brs.bandreps) do br
                BandRep(
                    br.wyckpos,
                    br.sitesym,
                    br.label,
                    br.dim,
                    br.decomposable,
                    br.spinful,
                    copy(br.irvec),
                    irlabs´
                )
            end,
            copy(brs.kvs),
            copy(brs.klabs),
            irlabs´,
            brs.allpaths,
            brs.spinful,
            brs.timereversal
        ),
        v.indicators,
        v.indicator_group
    )
    return prune_klab_irreps_vecs!(v´, klab)
end

function pick_klab_irreps_vecs!(v::BandSummary, klab::String="Γ")
    prune_iridxs = findall(irlab -> klabel(irlab) != klab, v.brs.irlabs)
    isempty(prune_iridxs) && error(lazy"could not find $klab among included irreps")
    deleteat!(v.n, prune_iridxs)

    pick_klab_irreps_brs!(v.brs, klab)

    return v
end

function pick_klab_irreps_vecs(v::BandSummary, klab::String="Γ")
    brs = v.brs
    irlabs´ = copy(brs.irlabs)
    v´ = BandSummary(
        v.topology,
        v.bands,
        copy(v.n),
        BandRepSet(
            brs.sgnum,
            map(brs.bandreps) do br
                BandRep(
                    br.wyckpos,
                    br.sitesym,
                    br.label,
                    br.dim,
                    br.decomposable,
                    br.spinful,
                    copy(br.irvec),
                    irlabs´
                )
            end,
            copy(brs.kvs),
            copy(brs.klabs),
            irlabs´,
            brs.allpaths,
            brs.spinful,
            brs.timereversal
        ),
        v.indicators,
        v.indicator_group
    )
    return pick_klab_irreps_vecs!(v´, klab)
end

function obtain_symmetry_vectors(ms::PyObject, sg_num::Int)
    brs = bandreps(sg_num) # elementary band representations
    lgs = littlegroups(sg_num) # little groups
    filter!(((klab, _),) -> klab ∈ klabels(brs), lgs) # restrict to k-points in `brs`
    map!(lg -> primitivize(lg, false), values(lgs)) # convert to primitive setting (without reducing translations)
    lgirsd = pick_lgirreps(lgs; timereversal=true) # small irreps associated with `lgs`

    symeigsd = Dict{String,Vector{Vector{ComplexF64}}}()
    for (klab, lg) in lgs
        kv = mp.Vector3(position(lg)()...)
        ms.solve_kpoint(kv)

        symeigsd[klab] = [Vector{ComplexF64}(undef, length(lg)) for n in 1:ms.num_bands]
        for (i, gᵢ) in enumerate(lg)
            W = mp.Matrix(eachcol(rotation(gᵢ))...) # decompose gᵢ = {W|w}
            w = mp.Vector3(translation(gᵢ)...)
            symeigs = ms.compute_symmetries(W, w) # compute ⟨Eₙₖ|gᵢDₙₖ⟩ for all bands
            setindex!.(symeigsd[klab], symeigs, i) # update container of symmetry eigenvalues
        end
    end

    # --- fix singular photonic symmetry content at Γ, ω=0 --- # TODO: Maybe we can skip this step or try to change the condition on Γ
    fixup_gamma_symmetry!(symeigsd, lgs)

    # --- analyze connectivity and topology of symmetry data ---
    summaries = analyze_symmetry_data(symeigsd, lgirsd, brs)

    return summaries
end

function find_auxiliary_modes(t::Int, d::Vector{Int64}, brs::BandRepSet)
    long_cand = PBC.filling_symmetry_constrained_expansions(t, Int[], d, brs, Int[])

    return long_cand
end

function physical(vᵀ::BandSummary, nᵀ⁺ᴸ, nᴸ)

    return all(>=(0), vᵀ.n - (nᵀ⁺ᴸ - nᴸ))
end

function find_all_band_representations(vᵀ::BandSummary, long_modes::Vector{Vector{Int64}}, d::Vector{Int64}, brs::BandRepSet)
    brs´ = prune_klab_irreps_brs(brs, "Γ")
    vᵀ´ = prune_klab_irreps_vecs(vᵀ, "Γ")
    idxs = collect(1:size(matrix(brs´), 1))

    output = Tuple{Vector{Vector{Int64}},Vector{Int64},Vector{Bool}}[]
    for i in 1:length(long_modes)
        nᴸ = long_modes[i]
        vᴸ´ = sum(brs´[nᴸ])
        vᵀ⁺ᴸ´ = vᵀ´.n + vᴸ´
        μᵀ⁺ᴸ = vᵀ⁺ᴸ´[end]

        nᵀ⁺ᴸ = PBC.filling_symmetry_constrained_expansions(μᵀ⁺ᴸ, vᵀ⁺ᴸ´, d, brs´, idxs)

        if nᵀ⁺ᴸ != []
            phys = [physical(vᵀ, sum(brs[j]), sum(brs[nᴸ])) for j in nᵀ⁺ᴸ]
            push!(output, (nᵀ⁺ᴸ, nᴸ, phys))
        end
    end
    return output
end

# Computes the generalized inverse `Xᵍ` of `X`, computed from the Smith normal form.
function generalized_inv(X::AbstractMatrix{<:Integer})
    F = smith(X)
    Λ = MPBUtils.diagm(F)
    Λg = zeros(Float64, size(Λ)[2], size(Λ)[1])
    for (n, λₙ) in enumerate(F.SNF)
        Λg[n, n] = iszero(λₙ) ? λₙ : inv(λₙ)
    end
    Xᵍ = F.Tinv * Λg * F.Sinv # generalized inverse

    return Xᵍ
end