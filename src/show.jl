function Base.show(io::IO, ::MIME"text/plain", candidates::TightBindingCandidates)
    summary(io, candidates)
    println(io, ":")
    for is in 1:length(candidates.long_modes)
        nᴸs = [candidates.bandreps[j] for j in candidates.long_modes[is]]
        print(io, " nᴸ = ")

        labs = String[]
        for j in nᴸs
            lab = "(" * j.wyckpos * "|" * replace(j.label, "↑G"=>"") * ")"
            push!(labs, lab)
        end
        join(io, labs, " + ")
        println(io)

        nᵀ⁺ᴸs = [[candidates.bandreps[m] for k in js for m in k] for js in candidates.solutions[is]]
        for (count, j) in enumerate(nᵀ⁺ᴸs)
            printstyled(io, "   ⁽", Crystalline.supscriptify(string(count)), "⁾ "; color=:light_black)
            print(io, "nᵀ⁺ᴸ = ")
            labsᵀ⁺ᴸ = String[]
            for k in j
                lab = "(" * k.wyckpos * "|" * replace(k.label, "↑G"=>"") * ")"
                push!(labsᵀ⁺ᴸ, lab)
            end
            join(io, labsᵀ⁺ᴸ, " + ")
            println(io)
        end
    end
end