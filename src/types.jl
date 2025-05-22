"""
    TightBindingCandidateSet

A structure for storing information about a set of tight-binding candidates for the photonic
crystal.

## Fields
- `longitudinal :: CompositeBandRep{3}`: the longitudinal band representation `nᴸ`
- `apolarv :: Vector{CompositeBandRep{3}}`: a vector of apolar band representations
  `nᵀ⁺ᴸ[i]` over `i ∈ eachindex(apolarv)`
- `ps :: Vector{Vector{Float64}}`: a vector of vectors of coefficients `ps[i]` associated to
    `nᵀ⁺ᴸ[i]`. The coefficients are used to trivialize the transverse band representations
    at ω=0, i.e., `nᵀ(ω=0) = n_fixed + Q*p`, where `n_fixed` is the fixed part shown in Tables
    (S6-8) in [this article](https://link.aps.org/doi/10.1103/PhysRevX.12.021066).
"""
struct TightBindingCandidateSet <: AbstractVector{CompositeBandRep{3}}
    longitudinal::CompositeBandRep{3}    # `nᴸ`
    apolarv::Vector{CompositeBandRep{3}} # `nᵀ⁺ᴸ[i]` over `i ∈ eachindex(apolarv)`
    ps::Vector{Vector{Float64}}          # `ps[i]` associates to `nᵀ⁺ᴸ[i]`; w/ nᵀ(ω=0) = n_fixed + Q*p
end
Base.size(tbc::TightBindingCandidateSet) = (length(tbc.apolarv),)
longitudinal(tbc::TightBindingCandidateSet) = tbc.longitudinal
Base.getindex(tbc::TightBindingCandidateSet, i::Int) = tbc.apolarv[i]
Base.setindex!(::TightBindingCandidateSet, v, i::Int) = error("setindex! is not supported")