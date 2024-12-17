struct TightBindingCandidateSet <: AbstractVector{CompositeBandRep{3}}
    longitudinal::CompositeBandRep{3}    # `nᴸ`
    apolarv::Vector{CompositeBandRep{3}} # `nᵀ⁺ᴸ[i]` over `i ∈ eachindex(apolarv)`    
    ps::Vector{Vector{Float64}}          # `ps[i]` associates to `nᵀ⁺ᴸ[i]`; w/ nᵀ(ω=0) = n_fixed + Q*p
end
Base.size(tbc::TightBindingCandidateSet) = (length(tbc.apolarv),)
longitudinal(tbc::TightBindingCandidateSet) = tbc.longitudinal
Base.getindex(tbc::TightBindingCandidateSet, i::Int) = tbc.apolarv[i]
Base.setindex!(::TightBindingCandidateSet, v, i::Int) = error("setindex! is not supported")

struct SymmetricHopping
    representatives::RVec # the representative of the hopping distance
    orbit::Vector{RVec} # the orbit of the hopping distance
    hop_terms::Vector{Vector{NTuple{3,RVec}}} # the hopping terms `(a,b,R)` associated to 
    #                                           each hopping distance
end
Base.size(s::SymmetricHopping) = (length(s.orbit),)
representatives(s::SymmetricHopping) = s.representatives
orbit(s::SymmetricHopping) = s.orbit
hop(s::SymmetricHopping) = s.hop_terms
Base.getindex(s::SymmetricHopping, i::Int) = s.hop_terms[i]
Base.setindex(s::SymmetricHopping, v, i::Int) = error("setindex! is not supported")