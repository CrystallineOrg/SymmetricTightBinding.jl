"""
    bandsum2symvec(bs::BandSummary, lgirsv::Vector{Collection{LGIrrep{D}}}) where D
    bandsum2symvec(bs::BandSummary, brs::Collection{<:NewBandRep}) --> Vector{SymmetryVector{D}}

Converts a `BandSummary` into a `SymmetryVector` using the irreps of the bands in `lgirsv` or `brs`.
"""
function bandsum2symvec(bs::BandSummary, lgirsv::Vector{Collection{LGIrrep{D}}}) where {D}
    multsv = [zeros(Int, length(lgirs)) for lgirs in lgirsv]
    for (i, lgirs) in enumerate(lgirsv)
        for (j, lgir) in enumerate(lgirs)
            ir_lab = label(lgir)
            idx = findfirst(==(ir_lab), bs.brs.irlabs)
            isnothing(idx) && error(lazy"could not find irrep $ir_lab in BandSummary")
            multsv[i][j] = bs.n[idx]
        end
    end
    μ = length(bs.bands)

    return SymmetryVector{D}(lgirsv, multsv, μ) #= occupation =#
end
function bandsum2symvec(bs::BandSummary, brs::Collection{<:NewBandRep})
    return bandsum2symvec(bs, irreps(first(brs)))
end
