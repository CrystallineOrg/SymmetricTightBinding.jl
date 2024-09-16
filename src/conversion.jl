# Convert a `BandSummary` to a `SymmetryVector`
function bandsum2symvec(bs::BandSummary, lgirsv::Vector{Collection{LGIrrep{D}}}) where D
    multsv = [zeros(Int, length(lgirs)) for lgirs in lgirsv]
    for (i, lgirs) in enumerate(lgirsv)
        for (j, lgir) in enumerate(lgirs)
            irlab = label(lgir)
            idx = findfirst(==(irlab), bs.brs.irlabs)
            isnothing(idx) && error(lazy"could not find irrep $irlab in BandSummary")
            multsv[i][j] = bs.n[idx]
        end
    end
    μ = length(bs.bands)

    return SymmetryVector{D}(lgirsv, multsv, #= occupation =# μ)
end

function bandsum2symvec(bs::BandSummary, brs :: Collection{<:NewBandRep})
    return bandsum2symvec(bs, irreps(first(brs)))
end

""" 
TODO

"""
function string2symvec(s::String, lgirsv::Vector{Collection{LGIrrep{D}}}) where D
    s = replace(s, " "=>"", "*"=>"")
    multsv = [zeros(Int, length(lgirs)) for lgirs in lgirsv]
    for (i, lgirs) in enumerate(lgirsv)
        for (j, lgir) in enumerate(lgirs)
            irlab = label(lgir)
            idxs = findfirst(irlab, s)
            isnothing(idxs) && continue
            pos₂ = first(idxs)
            m = Crystalline.searchpriornumerals(s, pos₂, Int)
            multsv[i][j] = m
        end
    end
    μ = sum(eachindex(first(lgirsv))) do j
        first(multsv)[j]*Crystalline.irdim(first(lgirsv)[j])
    end
    SymmetryVector(lgirsv, multsv, μ)
end