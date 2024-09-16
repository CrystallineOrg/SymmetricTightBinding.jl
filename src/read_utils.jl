struct TightBindingCandidates
    solutions::Vector{Vector{Vector{Int64}}} # all posible solutions nᵀ⁺ᴸ for each long_modes
    # nᴸ 
    long_modes::Vector{Vector{Int64}} # all posible long_modes nᴸ with solutions
    phys::Vector{Vector{Bool}} # physicallity of each solution nᵀ⁺ᴸ
    p::Vector{Vector{Vector{Float64}}} # vector of integers neccesary to compute the surrogate
    # representation
    bandreps::BandRepSet # for printing pretty
end

function Base.show(io::IO, ::MIME"text/plain", candidates::TightBindingCandidates)
    for i in 1:length(candidates.long_modes)
        nᴸ = [candidates.bandreps[j...] for j in candidates.long_modes[i]]
        print(io, "Solutions using the auxiliary mode: ")

        for j in nᴸ[1:end-1]
            print(io, j.label, " at ", j.wyckpos, " ⊕ ")
        end

        print(io, nᴸ[end].label, " at ", nᴸ[end].wyckpos, "\n")

        nᵀ⁺ᴸ = [[candidates.bandreps[m...] for k in j for m in k] for j in candidates.solutions[i]]
        count = 1
        for j in nᵀ⁺ᴸ
            print(io, "   ↪Solution #$count: ")
            for k in j[1:end-1]
                print(io, k.label, " at ", k.wyckpos, " ⊕ ")
            end
            print(io, j[end].label, " at ", j[end].wyckpos, "\n")
            count += 1
        end
    end
end

struct PhysicalTightBindingCandidates
    solutions::Vector{Vector{Vector{Int64}}} # all posible solutions nᵀ⁺ᴸ for each long_modes
    # nᴸ 
    long_modes::Vector{Vector{Int64}} # all posible long_modes nᴸ with solutions
    p::Vector{Vector{Vector{Float64}}} # vector of integers neccesary to compute the surrogate
    # representation
end