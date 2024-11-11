struct TightBindingCandidateSet <: AbstractVector{CompositeBandRep{3}}
    longitudinal::CompositeBandRep{3}    # `nᴸ`
    apolarv::Vector{CompositeBandRep{3}} # `nᵀ⁺ᴸ[i]` over `i ∈ eachindex(apolarv)`    
    ps::Vector{Vector{Float64}}          # `ps[i]` associates to `nᵀ⁺ᴸ[i]`; w/ nᵀ(ω=0) = n_fixed + Q*p
end

Base.size(tbc::TightBindingCandidateSet) = (length(tbc.apolarv),)
longitudinal(tbc::TightBindingCandidateSet) = tbc.longitudinal
Base.getindex(tbc::TightBindingCandidateSet, i::Int) = tbc.apolarv[i]
Base.setindex!(::TightBindingCandidateSet, v, i::Int) = error("setindex! is not supported")

struct HoppingClass <: AbstractVector{RVec{D}} where {D}
    distance::RVec{D} # δ = wi + R - qj 
    terms::Tuple{WyckoffPosition{D},WyckoffPosition{D}} # WPs concerning the hopping (with R included)
end