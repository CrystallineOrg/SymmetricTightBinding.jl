const TBC_ELTYPE = @NamedTuple{longitudinal::Vector{NewBandRep{3}},
    apolarv::Vector{Vector{NewBandRep{3}}},
    physicalv::Vector{Bool},
    pv::Vector{Vector{Float64}}}

struct TightBindingCandidates <: AbstractVector{TBC_ELTYPE}
    solutions::Vector{Vector{Vector{Int64}}} # all posible solutions nᵀ⁺ᴸ for each long_modes
    # nᴸ 
    long_modes::Vector{Vector{Int64}} # all posible long_modes nᴸ with solutions
    phys::Vector{Vector{Bool}} # physicallity of each solution nᵀ⁺ᴸ
    p::Vector{Vector{Vector{Float64}}} # vector of integers neccesary to compute the surrogate
    # representation
    brs::Collection{NewBandRep{3}}
end
Base.size(tbc::TightBindingCandidates) = (length(tbc.long_modes),)
function Base.getindex(tbc::TightBindingCandidates, i::Int)
    longitudinal = tbc.brs[tbc.long_modes[i]]
    apolarv = [tbc.brs[brs_idx] for brs_idx in tbc.solutions[i]]
    physicalv = tbc.phys[i]
    pv = tbc.p[i]
    return (; longitudinal, apolarv, physicalv, pv)
end
Base.setindex!(::TightBindingCandidates, v, i::Int) = error("setindex! is not supported")

struct PhysicalTightBindingCandidates
    solutions::Vector{Vector{Vector{Int64}}} # all posible solutions nᵀ⁺ᴸ for each long_modes
    # nᴸ 
    long_modes::Vector{Vector{Int64}} # all posible long_modes nᴸ with solutions
    p::Vector{Vector{Vector{Float64}}} # vector of integers neccesary to compute the surrogate
    # representation
    brs::Collection{NewBandRep{3}}
end
function Base.getindex(tbc::PhysicalTightBindingCandidates, i::Int)
    longitudinal = tbc.brs[tbc.long_modes[i]]
    apolarv = [tbc.brs[brs_idx] for brs_idx in tbc.solutions[i]]
    pv = tbc.p[i]
    return (; longitudinal, apolarv, pv)
end