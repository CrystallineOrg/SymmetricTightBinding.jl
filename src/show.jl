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

function Base.show(io::IO, ::MIME"text/plain", candidates::TightBindingCandidateSet)
    summary(io, candidates)
    println(io, ":")

    printstyled(io, "nᴸ"; bold=true)
    print(io, " = ")
    join(io, auxiliary(candidates), " + ")
    println(io)

    for (j, idxsᵀ⁺ᴸ) in enumerate(candidates.idxsᵀ⁺ᴸs)
        printstyled(io, "⁽", Crystalline.supscriptify(string(j)), "⁾ ";
            color=:light_black)
        printstyled(io, "nᵀ⁺ᴸ"; bold=true)
        print(io, " = ")
        join(io, candidates.brs[idxsᵀ⁺ᴸ], " + ")

        printstyled(io, " (𝐩 = "; color=:light_black)
        if isapprox(candidates.ps[j], round.(candidates.ps[j]), atol=1e-10)
            printstyled(io, round.(Int, candidates.ps[j]); color=:light_black)
        else
            printstyled(io, candidates.ps[j]; color=:light_red)
        end
        printstyled(io, ")"; color=:light_black)
        j ≠ length(candidates) && println(io)
    end
end