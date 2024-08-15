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