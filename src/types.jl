const TBC_ELTYPE = @NamedTuple{longitudinal::Vector{NewBandRep{3}},
    apolarv::Vector{Vector{NewBandRep{3}}},
    physicalv::Vector{Bool},
    pv::Vector{Vector{Float64}}}

struct TightBindingCandidates <: AbstractVector{TBC_ELTYPE}
    idxsᵀ⁺ᴸss::Vector{Vector{Vector{Int64}}} # all posible solutions nᵀ⁺ᴸ for each idxsᴸs
    # nᴸ 
    idxsᴸs::Vector{Vector{Int64}} # all posible idxsᴸs nᴸ with solutions
    p::Vector{Vector{Vector{Float64}}} # vector of integers neccesary to compute the surrogate
    # representation
    brs::Collection{NewBandRep{3}}
end

Base.size(tbc::TightBindingCandidates) = (length(tbc.idxsᴸs),)
function Base.getindex(tbc::TightBindingCandidates, i::Int)
    longitudinal = tbc.brs[tbc.idxsᴸs[i]]
    apolarv = [tbc.brs[brs_idx] for brs_idx in tbc.idxsᵀ⁺ᴸss[i]]
    pv = tbc.p[i]
    return (; longitudinal, apolarv, pv)
end
Base.setindex!(::TightBindingCandidates, v, i::Int) = error("setindex! is not supported")