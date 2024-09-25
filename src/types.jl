struct TightBindingCandidateSet <: AbstractVector{Vector{NewBandRep{3}}}
    idxsᴸ::Vector{Int}            # `nᴸ = sum(brs[idxsᴸ])`
    idxsᵀ⁺ᴸs::Vector{Vector{Int}} # `nᵀ⁺ᴸ[i] = sum(brs[idxsᵀ⁺ᴸs[i]))` (specific to `nᴸ`)
    ps::Vector{Vector{Float64}}   # `ps[i]` associates to `nᵀ⁺ᴸ[i]`; w/ nᵀ(ω=0) = n_fixed + Q*p
    brs::Collection{NewBandRep{3}}
end

Base.size(tbc::TightBindingCandidateSet) = (length(tbc.idxsᵀ⁺ᴸs),)
auxiliary(tbc::TightBindingCandidateSet) = tbc.brs[tbc.idxsᴸ]
function Base.getindex(tbc::TightBindingCandidateSet, i::Int)
    apolar_decomposition = [tbc.brs[brs_idx] for brs_idx in tbc.idxsᵀ⁺ᴸs[i]]
    return apolar_decomposition
end
Base.setindex!(::TightBindingCandidateSet, v, i::Int) = error("setindex! is not supported")