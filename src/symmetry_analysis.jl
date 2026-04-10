# Note [⚠️ phase]: `symmetry_eigenvalues` returns the complex conjugate of the Convention 1
#   character to match Crystalline.jl's `calc_bandreps` convention.
#   See `docs/src/devdocs/symmetry_eigenvalue_conventions.md`.

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

    clgirsv = irreps(cbr) # irreps associated to the EBRs (conventional setting operations)
    lgirsv = primitivize.(clgirsv) # must be `modw=false` (default for Collection dispatch)
    lgs = group.(lgirsv)  # little groups associated to the EBRs (primitive setting)
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

!!! note
    The inputs `ops`, `k`, and `lg` must be provided in a primitive setting. See
    Crystalline.jl's `primitivize`.

!!! warning "⚠️ character phase convention"
    The symmetry eigenvalues returned by this function are the complex conjugate of the
    Convention 1 result (see `docs/src/devdocs/symmetry_eigenvalue_conventions.md`)
    in order to match the convention used by Crystalline.jl's `calc_bandreps` and `lgirreps`
    functions. See the above documentation for more details and 
    https://github.com/thchr/Crystalline.jl/issues/12 for the relevant issue in
    Crystalline.jl.
"""
function symmetry_eigenvalues(
    ptbm::ParameterizedTightBindingModel{D},
    ops::AbstractVector{SymOperation{D}},
    k::ReciprocalPointLike{D},
    sgreps::AbstractVector{SiteInducedSGRepElement{D}} = begin
        sgrep_induced_by_siteir.(Ref(ptbm.tbm.cbr), ops)
    end,
) where D
    length(k) == D || error("dimension mismatch")
    length(sgreps) == length(ops) || error("length of `sgreps` must match length of `ops`")

    # NB: `solve` with `bloch_phase=Val(false)` returns eigenvectors `vs` of H(k) in the
    #     Convention 1 coefficient basis (without Bloch position phases). In Convention 1,
    #     the symmetry eigenvalues are then `χ[n] = (Θ_G vs[n])† D_k vs[n]` where Θ_G & D_k
    #     defined in `docs/src/theory.md` and `docs/src/devdocs/` (and methods below).
    #
    # [⚠️ phase]: Crystalline.jl's `calc_bandreps` and `lgirreps` computes characters in a
    #     convention that is the complex conjugate of the Convention 1 result (see
    #     thchr/Crystalline.jl/#12).
    #     To be able to interface with Crystalline.jl, and until thchr/Crystalline.jl/#12 is
    #     resolved, we thus actually return `χ_Crystalline = conj(χ_Convention1)`.
    _, vs = solve(ptbm, k; bloch_phase = Val(false))
    symeigs = Matrix{ComplexF64}(undef, length(ops), ptbm.tbm.N)
    v_kpG = similar(vs, size(vs, 1)) # preallocate for Θᴳ * v
    for (j, sgrep) in enumerate(sgreps)
        g = sgrep.op
        gk = compose(g, ReciprocalPoint{D}(k)) # NB: for k ∈ Gₖ, there exist G st g∘k = k+G
        G = gk - k # the possible reciprocal vector-difference G between k & g∘k; for Θᴳ
        Θᴳ = reciprocal_translation_phase(orbital_positions(ptbm), G)
        D_k = sgrep(k) # = D_k(g) = e^{-2πi(gk)·t} ρ(h) (Convention 1)
        for (n, v) in enumerate(eachcol(vs))
            v_kpG = mul!(v_kpG, Θᴳ, v) # = Θᴳ * v (without re-allocating `v_kpG`)
            χ = dot(v_kpG, D_k, v)  # Convention 1: (Θ_G w)† D_k w
            χ_Crystalline = conj(χ) # [⚠️ phase]: convert to Crystalline.jl's convention
            symeigs[j, n] = χ_Crystalline
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
    clgirsv = irreps(ptbm.tbm.cbr) # irreps associated to the EBRs (conventional setting)
    lgirsv = primitivize.(clgirsv) # convert associated groups & irreps to primitive setting
    symeigsv = [eachcol(symmetry_eigenvalues(ptbm, group(lgirs))) for lgirs in lgirsv]
    return collect_irrep_annotations(symeigsv, lgirsv; kws...)
end