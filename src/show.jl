function Base.show(io::IO, ::MIME"text/plain", candidates::TightBindingCandidates)
    summary(io, candidates)
    println(io, ":")
    for is in 1:length(candidates.long_modes)
        nᴸs = [candidates.brs[j] for j in candidates.long_modes[is]]
        print(io, " nᴸ = ")

        join(io, nᴸs, " + ")
        println(io)

        nᵀ⁺ᴸss = [[candidates.brs[m] for k in js for m in k] for js in candidates.solutions[is]]
        for (count, nᵀ⁺ᴸs) in enumerate(nᵀ⁺ᴸss)
            printstyled(io, "   ⁽", Crystalline.supscriptify(string(count)), "⁾ "; color=:light_black)
            print(io, "nᵀ⁺ᴸ = ")
            join(io, nᵀ⁺ᴸs, " + ")
            println(io)
        end
    end
end