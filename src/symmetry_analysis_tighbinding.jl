"""
    symmetry_analysis(ptbm::ParameterizedTightBindingModel{D}; multiplicities_kws...)

Determine a decomposition of the bands associated with `ptbm` into a set of
`SymmetryVector`s, with each symmetry vector corresponding to a set of
compatibility-respecting (i.e., energy separable along high-symmetry **k**-lines) bands.

## Keyword arguments
- `multiplicities_kws...`: keyword arguments passed to `MPBUtils.analyze_symmetry_data`
  used in determining the multiplicities of irreps across high-symmetry **k**-points.

## Example
```jl-repl
julia> using Crystalline, TETB

julia> brs = calc_bandreps(221); coefs = zeros(Int, length(brs)); coefs[[1, 2]] .= 1;

julia> cbr = CompositeBandRep(coefs, brs)
40-irrep CompositeBandRep{3}:
 (3d|A₁g) + (3d|A₁ᵤ) (6 bands)

julia> tbm = tb_hamiltonian(cbr, [zeros(Int, dim(cbr))]); # a 4-term tight-binding model

julia> ptbm = tbm([1.0, 0.2, -1.0, 0.3]); # fix free coefficients

julia> symmetry_analysis(ptbm)
2-element Vector{SymmetryVector{3}}:
 [M₅⁺+M₁⁻, X₃⁺+X₁⁻+X₂⁻, Γ₁⁻+Γ₃⁻, R₄⁺] (3 bands)
 [M₁⁺+M₅⁻, X₁⁺+X₂⁺+X₃⁻, Γ₁⁺+Γ₃⁺, R₄⁻] (3 bands)
```
In the above example, the bands separate into two symmetry vectors, one for each of the
original EBRs in `cbr`.
"""
function symmetry_analysis(
    ptbm::ParameterizedTightBindingModel{D},
    multiplicities_kws...,
) where D
    tbm = ptbm.tbm
    isempty(tbm.terms) && error("`ptbm` is an empty tight-binding model")
    cbr = tbm.cbr # CompositeBandRep

    lgirsv = irreps(cbr) # get irreps associated to the EBRs
    lgs = [group(first(lgirs)) for lgirs in lgirsv]
    ops = unique(Iterators.flatten(lgs))

    # determine the induced space group rep associated with `cbr` across all `ops`
    sgrep_d = Dict(op => sgrep_induced_by_siteir(cbr, op) for op in ops)

    N = ptbm.tbm.N # number of bands in `ptbm`
    symeigsd = Dict{String, Vector{Vector{ComplexF64}}}()
    positions = _positions_in_bandrep_induced_sgrep(cbr) # for Bloch phases (NB: we already have this in elements of `sgrep_d`)

    # determine symmetry eigenvalues for each band in each little group
    for lg in lgs
        k = constant(position(lg))
        phases = cispi.(2 .* dot.(Ref(k), positions)) # e^{ik·r}

        H = Hermitian(ptbm(k))
        es, vs = eigen(H)
        vs = Diagonal(phases) * vs # add Bloch phases
        symeigs = [Vector{ComplexF64}(undef, length(lg)) for _ in 1:N]
        for (j, op) in enumerate(lg)
            ρ = sgrep_d[op](k)
            for (n, v) in enumerate(eachcol(vs))
                symeigs[n][j] = dot(v, transpose(ρ) * v)
            end
        end
        symeigsd[klabel(lg)] = symeigs
    end

    # TODO: transfer these utilities from MBPUtils to Crystalline, and also
    # to drop their dependence on taking a BandRepSet rather than a Collection{NewBandRep}
    _brs = convert(BandRepSet, cbr.brs)
    lgirsd = Dict(klabel(first(lgirs)) => lgirs for lgirs in irreps(cbr))
    summaries =
        MPBUtils.analyze_symmetry_data(symeigsd, lgirsd, _brs; multiplicities_kws...)
    # TODO: ↓ Ideally, wouldn't be needed once we've transferred from MPBUtils to Crystalline
    ns = bandsum2symvec.(summaries, Ref(lgirsv))

    return ns
end
