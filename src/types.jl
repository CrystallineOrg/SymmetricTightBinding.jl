struct TightBindingCandidates
    solutions::Vector{Vector{Vector{Int64}}} # all posible solutions nᵀ⁺ᴸ for each long_modes
    # nᴸ 
    long_modes::Vector{Vector{Int64}} # all posible long_modes nᴸ with solutions
    phys::Vector{Vector{Bool}} # physicallity of each solution nᵀ⁺ᴸ
    p::Vector{Vector{Vector{Float64}}} # vector of integers neccesary to compute the surrogate
    # representation
    bandreps::BandRepSet
end

struct PhysicalTightBindingCandidates
    solutions::Vector{Vector{Vector{Int64}}} # all posible solutions nᵀ⁺ᴸ for each long_modes
    # nᴸ 
    long_modes::Vector{Vector{Int64}} # all posible long_modes nᴸ with solutions
    p::Vector{Vector{Vector{Float64}}} # vector of integers neccesary to compute the surrogate
    # representation
end
