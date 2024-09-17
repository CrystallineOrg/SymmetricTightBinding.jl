function Base.show(io::IO, ::MIME"text/plain", candidates::TightBindingCandidates)
    summary(io, candidates)
    println(io, ":")
    for is in 1:length(candidates.long_modes)
        nᴸs = [candidates.brs[j] for j in candidates.long_modes[is]]
        print(io, " nᴸ = ")

        join(io, nᴸs, " + ")
        println(io)

        nᵀ⁺ᴸss = [[candidates.brs[m] for k in js for m in k] for js in candidates.solutions[is]]
        for js in eachindex(nᵀ⁺ᴸss)
            printstyled(io, "   ⁽", Crystalline.supscriptify(string(js)), "⁾ ";
                color=:light_black)
            print(io, "nᵀ⁺ᴸ = ")
            join(io, nᵀ⁺ᴸss[js], " + ")
            print(io, "; Physical? = ", candidates.phys[js][1])
            print(io, "; 𝐩 = ", candidates.p[js][1])
            println(io)
        end
    end
end

function Base.show(io::IO, ::MIME"text/plain", candidates::PhysicalTightBindingCandidates)
    summary(io, candidates)
    println(io, ":")
    for is in 1:length(candidates.long_modes)
        nᴸs = [candidates.brs[j] for j in candidates.long_modes[is]]
        print(io, " nᴸ = ")

        join(io, nᴸs, " + ")
        println(io)

        nᵀ⁺ᴸss = [[candidates.brs[m] for k in js for m in k] for js in candidates.solutions[is]]
        for js in eachindex(nᵀ⁺ᴸss)
            printstyled(io, "   ⁽", Crystalline.supscriptify(string(js)), "⁾ ";
                color=:light_black)
            print(io, "nᵀ⁺ᴸ = ")
            join(io, nᵀ⁺ᴸss[js], " + ")
            print(io, "; 𝐩 = ", candidates.p[js][1])
            println(io)
        end
    end
end