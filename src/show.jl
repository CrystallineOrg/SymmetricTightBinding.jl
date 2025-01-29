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
    # before getting started, determine maximum length of δᵢ entries, to align:
    aligns = map(enumerate(ho.orbit)) do (i, δᵢ)
        s = sprint(print, Crystalline.subscriptify(string(i)), δᵢ)
        textwidth(s)
    end
    max_align = maximum(aligns)
    # now print info about each orbit element and its hopping terms
    print(io, typeof(ho), " (")
    printstyled(io, "a"; color=:green)
    print(io, " + δ = ")
    printstyled(io, "b"; color=:red)
    print(io, " + ")
    printstyled(io, "R"; color=:blue)
    println(io, "):")
    for (i, (δᵢ, abRs)) in enumerate(zip(ho.orbit, ho.hoppings))
        print(io, " ")
        printstyled(io, "δ", Crystalline.subscriptify(string(i)), " = ", δᵢ; underline=i==1)
        print(io, )
        print(io, ": ", " "^(max_align-aligns[i]), "[")
        for (j, (a, b, R)) in enumerate(abRs)
            printstyled(io, "("; color=:light_black)
            printstyled(io, a; color=:green)
            print(io, " → ")
            printstyled(io, b; color=:red)
            print(io, " + ")
            printstyled(io, R; color=:blue)
            printstyled(io, ")"; color=:light_black)
            j ≠ length(abRs) && print(io, ", ")
        end
        print(io, "]")
        i ≠ length(ho.orbit) && println(io)
    end
end