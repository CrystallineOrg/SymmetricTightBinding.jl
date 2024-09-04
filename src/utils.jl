# 
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
    long_cand = find_all_admissible_expansions(
        brs, d, t, #= occupation =#
        Int[], Int[]) #= idxs =#

    return long_cand
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

function physical(vᵀᵧ::BandSummary, nᵀ⁺ᴸᵧ, nᴸᵧ, sg_num::Int)
    lgirs = realify(lgirreps(sg_num)["Γ"])
    _, Q = physical_zero_frequency_gamma_irreps(
        lgirs;
        supergroup_constraints=true,
        force_fixed=true,
        lattice_reduce=true)

    # sort Q in the same order as the SG irreps order
    irs = string.(lgirs.irs) # change them from 
    irlabs = vᵀᵧ.brs.irlabs
    perm = sortperm(irs)[invperm(sortperm(irlabs))]
    Q_ordered = Q[perm, :]

    Q⁻¹ = generalized_inv(Q_ordered)
    nᵀᵧ = nᵀ⁺ᴸᵧ - nᴸᵧ
    y = Q⁻¹ * (nᵀᵧ-vᵀᵧ.n)[1:end-1]

    return (all([(y[i] - round(y[i])) == 0 for i in eachindex(y)]), y)
end

function find_all_band_representations(vᵀ::BandSummary, long_modes::Vector{Vector{Int64}}, d::Vector{Int64}, brs::BandRepSet, sg_num::Int)
    brs´ = prune_klab_irreps_brs(brs, "Γ")
    vᵀ´ = prune_klab_irreps_vecs(vᵀ, "Γ")
    idxs = collect(1:size(matrix(brs´), 1))

    brsᵧ = pick_klab_irreps_brs(brs, "Γ")
    vᵀᵧ = pick_klab_irreps_vecs(vᵀ, "Γ")

    output = Tuple{Vector{Vector{Int64}},Vector{Int64},Any}[]
    for i in 1:length(long_modes)
        nᴸ = long_modes[i]
        vᴸ´ = sum(brs´[nᴸ])
        vᵀ⁺ᴸ´ = vᵀ´.n + vᴸ´
        μᵀ⁺ᴸ = vᵀ⁺ᴸ´[end]

        nᵀ⁺ᴸ = find_all_admissible_expansions(brs´, d, μᵀ⁺ᴸ, vᵀ⁺ᴸ´, idxs)

        if !isempty(nᵀ⁺ᴸ)
            phys = [physical(vᵀᵧ, sum(brsᵧ[j]), sum(brsᵧ[nᴸ]), sg_num) for j in nᵀ⁺ᴸ]
            push!(output, (nᵀ⁺ᴸ, nᴸ, phys))
        end
    end
    return output
end

function find_physical_band_representations(vᵀ::BandSummary, long_modes::Vector{Vector{Int64}}, d::Vector{Int64}, brs::BandRepSet, sg_num::Int)
    brs´ = prune_klab_irreps_brs(brs, "Γ")
    vᵀ´ = prune_klab_irreps_vecs(vᵀ, "Γ")
    idxs = collect(1:size(matrix(brs´), 1))

    brsᵧ = pick_klab_irreps_brs(brs, "Γ")
    vᵀᵧ = pick_klab_irreps_vecs(vᵀ, "Γ")

    output = Tuple{Vector{Int64},Vector{Int64},Any}[]
    for i in 1:length(long_modes)
        nᴸ = long_modes[i]
        vᴸ´ = sum(brs´[nᴸ])
        vᵀ⁺ᴸ´ = vᵀ´.n + vᴸ´
        μᵀ⁺ᴸ = vᵀ⁺ᴸ´[end]

        nᵀ⁺ᴸ = find_all_admissible_expansions(brs´, d, μᵀ⁺ᴸ, vᵀ⁺ᴸ´, idxs)

        if nᵀ⁺ᴸ != []
            for j in eachindex(nᵀ⁺ᴸ)
                phys = physical(vᵀᵧ, sum(brsᵧ[nᵀ⁺ᴸ[j]]), sum(brsᵧ[nᴸ]), sg_num)
                if phys[1]
                    push!(output, (nᵀ⁺ᴸ[j], nᴸ, phys[2]))
                end
            end
        end
    end
    return output
end