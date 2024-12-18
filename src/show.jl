function Base.show(io::IO, candidates::TightBindingCandidateSet)
    printstyled(io, "nᴸ"; bold=true, color=:light_black)
    printstyled(io, " = ", longitudinal(candidates), ": "; color=:light_black)

    printstyled(io, "nᵀ⁺ᴸ"; bold=true)
    print(io, " ∈ [")
    for (j, nᵀ⁺ᴸ) in enumerate(candidates)
        print(io, nᵀ⁺ᴸ)
        j ≠ length(candidates) && print(io, ", ")
    end
    print(io, "]")
end

function Base.show(io::IO, ::MIME"text/plain", candidates::TightBindingCandidateSet)
    summary(io, candidates)
    println(io, ":")

    printstyled(io, "nᴸ"; bold=true)
    print(io, " = ")
    print(io, longitudinal(candidates))
    println(io)

    for (j, nᵀ⁺ᴸ) in enumerate(candidates)
        printstyled(io, "⁽", Crystalline.supscriptify(string(j)), "⁾ ";
            color=:light_black)
        printstyled(io, "nᵀ⁺ᴸ"; bold=true)
        print(io, " = ")
        print(io, nᵀ⁺ᴸ)

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

function Base.show(io::IO, ::MIME"text/plain", ho::HoppingOrbit)
    println(io, typeof(ho), ":")
    for (i, (δᵢ, abRs)) in enumerate(zip(ho.orbit, ho.hoppings))
        printstyled(io, " δ", Crystalline.subscriptify(string(i)), " = ", δᵢ; bold=i == 1)
        print(io, ": [")
        for (j, (a, b, R)) in enumerate(abRs)
            print(io, "$a → $b + $R")
            j ≠ length(abRs) && print(io, ", ")
        end
        println(io, "]")
    end
end