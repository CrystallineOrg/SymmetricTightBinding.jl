"""
Erase the content of a certain HSP from the input object. This HSP is settle by default to Γ.
"""
function prune_klab_irreps!(brs::Collection{<:NewBandRep}, klab::String="Γ")
    prune_kidx = findfirst(==(klab), klabels(brs))
    isnothing(prune_kidx) && error(lazy"could not find $klab among included k-points")

    lgirsv = irreps(first(brs))
    foreach(brs) do br
        deleteat!(multiplicities(br.n), prune_kidx)
        @assert irreps(br.n) === lgirsv
    end
    deleteat!(lgirsv, prune_kidx)

    return brs
end

"""
Copies the symmetry vectors of the bandreps inside the collection so they are not overwritten
"""
function _symmetry_vector_shallow_copy(brs::Collection{<:NewBandRep})
    lgirsv′ = copy(irreps(first(brs)))
    brs′ = Collection(map(brs) do br
        NewBandRep(
            br.siteir,
            SymmetryVector(lgirsv′, copy(multiplicities(br)), occupation(br)),
            br.spinful,
            br.timereversal)
    end)
    return brs′
end
"""
Erase the content of a certain HSP from the input object. This HSP is settle by default to Γ.
"""
function prune_klab_irreps(brs::Collection{<:NewBandRep}, klab::String="Γ")
    return prune_klab_irreps!(_symmetry_vector_shallow_copy(brs), klab)
end

function prune_klab_irreps!(v::SymmetryVector, klab::String="Γ")
    prune_iridxs = findall(lgirs -> klabel(first(lgirs)) == klab, irreps(v))
    isempty(prune_iridxs) && error(lazy"could not find $klab among included irreps")
    deleteat!(multiplicities(v), prune_iridxs)
    deleteat!(irreps(v), prune_iridxs)
    return v
end
function prune_klab_irreps!(v::AbstractSymmetryVector, klab::String="Γ")
    prune_klab_irreps!(SymmetryVector(v), klab)
end

function prune_klab_irreps(v::SymmetryVector, klab::String="Γ")
    lgirsv´ = copy(irreps(v))
    multsv´ = copy(multiplicities(v))
    v´ = SymmetryVector(lgirsv´, multsv´, occupation(v))
    return prune_klab_irreps!(v´, klab)
end
function prune_klab_irreps(v::AbstractSymmetryVector, klab::String="Γ")
    prune_klab_irreps(SymmetryVector(v), klab)
end

"""
Picks the content of a certain HSP from the input object, erasing all the content for other 
HSPs. This HSP is settle by default to Γ.
"""
function pick_klab_irreps!(brs::Collection{<:NewBandRep}, klab::String="Γ")
    prune_kidx = findall(!=(klab), klabels(brs))
    isnothing(prune_kidx) && error(lazy"could not find $klab among included k-points")

    lgirsv = irreps(first(brs))
    foreach(brs) do br
        deleteat!(multiplicities(br.n), prune_kidx)
        @assert irreps(br.n) === lgirsv
    end
    deleteat!(lgirsv, prune_kidx)

    return brs
end

function pick_klab_irreps(brs::Collection{<:NewBandRep}, klab::String="Γ")
    return pick_klab_irreps!(_symmetry_vector_shallow_copy(brs), klab)
end


function pick_klab_irreps!(v::SymmetryVector, klab::String="Γ")
    prune_iridxs = findall(lgirs -> klabel(first(lgirs)) != klab, irreps(v))
    isempty(prune_iridxs) && error(lazy"could not find $klab among included irreps")
    deleteat!(multiplicities(v), prune_iridxs)
    deleteat!(irreps(v), prune_iridxs)
    return v
end
function pick_klab_irreps!(v::AbstractSymmetryVector, klab::String="Γ")
    pick_klab_irreps!(SymmetryVector(v), klab)
end


function pick_klab_irreps(v::SymmetryVector, klab::String="Γ")
    lgirsv´ = copy(irreps(v))
    multsv´ = copy(multiplicities(v))
    v´ = SymmetryVector(lgirsv´, multsv´, occupation(v))
    return pick_klab_irreps!(v´, klab)
end
function pick_klab_irreps(v::AbstractSymmetryVector, klab::String="Γ")
    pick_klab_irreps(SymmetryVector(v), klab)
end