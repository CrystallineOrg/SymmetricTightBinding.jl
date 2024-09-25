function Base.show(io::IO, candidates::TightBindingCandidateSet)
    printstyled(io, "nᴸ"; bold=true, color=:light_black)
    printstyled(io, " = ", join(auxiliary(candidates), "+"), ": "; color=:light_black)

    printstyled(io, "nᵀ⁺ᴸ"; bold=true)
    print(io, " ∈ [")
    for (j, idxsᵀ⁺ᴸ) in enumerate(candidates.idxsᵀ⁺ᴸs)
        join(io, candidates.brs[idxsᵀ⁺ᴸ], "+")
        j ≠ length(candidates) && print(io, ", ")
    end
    print(io, "]")
end
