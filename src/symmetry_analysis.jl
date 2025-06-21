"""
    collect_compatible(ptbm::ParameterizedTightBindingModel{D}; multiplicities_kws...)

Determine a decomposition of the bands associated with `ptbm` into a set of
`SymmetryVector`s, with each symmetry vector corresponding to a set of
compatibility-respecting (i.e., energy separable along high-symmetry **k**-lines) bands.

## Keyword arguments
- `multiplicities_kws...`: keyword arguments passed to `Crystalline.collect_compatible`
  used in determining the multiplicities of irreps across high-symmetry **k**-points.

## Example
```julia-repl
julia> using Crystalline, SymmetricTightBinding

julia> brs = calc_bandreps(221);

julia> cbr = @composite brs[1] + brs[2]
40-irrep CompositeBandRep{3}:
 (3d|A₁g) + (3d|A₁ᵤ) (6 bands)

julia> tbm = tb_hamiltonian(cbr); # a 4-term, 6-band model

julia> ptbm = tbm([1.0, 0.1, -1.0, 0.1]); # fix free coefficients

julia> collect_compatible(ptbm)
2-element Vector{SymmetryVector{3}}:
 [M₅⁺+M₁⁻, X₃⁺+X₁⁻+X₂⁻, Γ₁⁻+Γ₃⁻, R₄⁺] (3 bands)
 [M₁⁺+M₅⁻, X₁⁺+X₂⁺+X₃⁻, Γ₁⁺+Γ₃⁺, R₄⁻] (3 bands)
```
In the above example, the bands separate into two symmetry vectors, one for each of the
original EBRs in `cbr`.
"""
function Crystalline.collect_compatible(
    ptbm::ParameterizedTightBindingModel{D},
    multiplicities_kws...,
) where D
    tbm = ptbm.tbm
    isempty(tbm.terms) && error("`ptbm` is an empty tight-binding model")
    cbr = CompositeBandRep(tbm)

    lgirsv = irreps(cbr) # get irreps associated to the EBRs
    lgs = [group(first(lgirs)) for lgirs in lgirsv]
    ops = unique(Iterators.flatten(lgs))

    # determine the induced space group rep associated with `cbr` across all `ops`
    sgrep_d = Dict(op => sgrep_induced_by_siteir(ptbm, op) for op in ops)

    symeigsv = Vector{Vector{Vector{ComplexF64}}}(undef, length(lgs))

    # determine symmetry eigenvalues for each band in each little group
    for (kidx, lg) in enumerate(lgs)
        sgreps = [sgrep_d[op] for op in lg]
        symeigs = symmetry_eigenvalues(ptbm, lg, sgreps)
        symeigsv[kidx] = collect(eachcol(symeigs))
    end

    ns = collect_compatible(symeigsv, cbr.brs; multiplicities_kws...)
    return ns
end

"""
    symmetry_eigenvalues(
        ptbm::ParameterizedTightBindingModel{D},
        ops::AbstractVector{SymOperation{D}},
        k::ReciprocalPointLike{D},
        [sgreps::AbstractVector{SiteInducedSGRepElement{D}}]
    )
    symmetry_eigenvalues(
        ptbm::ParameterizedTightBindingModel{D},
        lg::LittleGroup{D},
        [sgreps::AbstractVector{SiteInducedSGRepElement{D}}]
    )
        --> Matrix{ComplexF64}
    
Compute the symmetry eigenvalues of a coefficient-parameterized tight-binding model `ptbm`
at the **k**-point `k` for the symmetry operations `ops`. A `LittleGroup` can also be
provided instead of `ops` and `k`.

Representations of the symmetry operations `ops` as acting on the orbitals of the
tight-binding setting can optionally be provided in `sgreps` (see `sgrep_induced_by_siteir`)
and are otherwise initialized by the function.

The symmetry eigenvalues are returned as a matrix, with rows running over the elements of
`ops` and columns running over the bands of `ptbm`.
"""
function symmetry_eigenvalues(
    ptbm::ParameterizedTightBindingModel{D},
    ops::AbstractVector{SymOperation{D}},
    k::ReciprocalPointLike{D},
    sgreps::AbstractVector{SiteInducedSGRepElement{D}} = sgrep_induced_by_siteir.(
        Ref(ptbm.tbm.cbr),
        ops,
    ),
) where D
    length(k) == D || error("dimension mismatch")
    length(sgreps) == length(ops) || error("length of `sgreps` must match length of `ops`")

    # NB: Currently, the site-symmetry induced reps assume the "Convention 2" Fourier
    #     transform, which does depend on "in-unit-cell" coordinates; so we must add
    #     such phases here. Longer term, we might want to not do that (change to 
    #     "Convention 1") though since it could simplify the induced rep phases to an 
    #     overall phase instead of the tᵦₐ-business
    # NOTE: maybe the last comment is not the best since the bloch periodic functions are 
    #       periodic in real space but not in reciprocal space. This means that we will
    #       need to correct such phases manually. Changing to "Convention 2", as we do now,
    #       might be the best option for simplicity.
    _, vs = solve(ptbm, k; bloch_phase = Val(true))
    symeigs = Matrix{ComplexF64}(undef, length(ops), ptbm.tbm.N)
    for (j, sgrep) in enumerate(sgreps)
        ρ = sgrep(k)
        for (n, v) in enumerate(eachcol(vs))
            symeigs[j, n] = dot(v, transpose(ρ), v)
        end
    end
    return symeigs
end

function symmetry_eigenvalues(
    ptbm::ParameterizedTightBindingModel{D},
    lg::LittleGroup{D},
    sgreps::AbstractVector{SiteInducedSGRepElement{D}} = sgrep_induced_by_siteir.(
        Ref(ptbm.tbm.cbr),
        lg,
    ),
) where D
    kv = position(lg)
    isspecial(kv) || error("input k-point has free parameters")
    k = constant(kv)
    return symmetry_eigenvalues(ptbm, operations(lg), k, sgreps)
end

"""
    collect_irrep_annotations(ptbm::ParameterizedTightBindingModel; kws...)

Collect the irrep labels across the high-symmetry **k**-points referenced by the underlying
composite band representation of `ptbm`, across the bands of the model.

Useful for annotating irrep labels in band structure plots (via the Makie extension call
`plot(ks, energies; annotations=collect_irrep_annotations(ptbm))`)
"""
function Crystalline.collect_irrep_annotations(ptbm::ParameterizedTightBindingModel; kws...)
    lgirsv = irreps(ptbm.tbm.cbr) # get irreps associated to the EBRs
    symeigsv = [eachcol(symmetry_eigenvalues(ptbm, group(lgirs))) for lgirs in lgirsv]
    return collect_irrep_annotations(symeigsv, lgirsv; kws...)
end