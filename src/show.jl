function Base.show(io::IO, ::MIME"text/plain", candidates::TightBindingCandidates)
    summary(io, candidates)
    println(io, ":")
    for is in 1:length(candidates.idxsᴸs)
        nᴸs = [candidates.brs[j] for j in candidates.idxsᴸs[is]]
        print(io, " nᴸ = ")

        join(io, nᴸs, " + ")
        println(io)

        nᵀ⁺ᴸss = [[candidates.brs[m] for k in js for m in k] for js in candidates.idxsᵀ⁺ᴸss[is]]
        for js in eachindex(nᵀ⁺ᴸss)
            printstyled(io, "   ⁽", Crystalline.supscriptify(string(js)), "⁾ ";
                color=:light_black)
            print(io, "nᵀ⁺ᴸ = ")
            join(io, nᵀ⁺ᴸss[js], " + ")
            printstyled(io, " (𝐩 = ", candidates.p[is][js], ")"; color=:light_black)
            println(io)
        end
    end
end
