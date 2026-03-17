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
    lgirsv = primitivize.(clgirsv)
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

    # NB: Currently, the site-symmetry induced reps assume the "Convention 1" Fourier
    #     transform. This Fourier transform does depend on "in-unit-cell" coordinates;
    #     so we must correct such phases here, as indicated in `/docs/usr/theory.md`.
    #
    # NOTE: since we picked "Convention 1" for the Fourier transform, we need to correct an
    #       extra phase factor to correct the non-periodicity of the Bloch functions under
    #       this convention.
    _, vs = solve(ptbm, k; bloch_phase = Val(false))
    symeigs = Matrix{ComplexF64}(undef, length(ops), ptbm.tbm.N)
    for (j, sgrep) in enumerate(sgreps)
        g = sgrep.op
        gk = compose(g, ReciprocalPoint{D}(k)) # NB: for k ∈ Gₖ, there exist G st g∘k = k+G
        G = gk - k # the possible reciprocal vector-difference G between k & g∘k; for Θᴳ
        # NB: we use -G (i.e., `conj(Θᴳ)`) rather than G because the symmetry eigenvalue
        #     formula `⟨ψ|ĝ|ψ⟩ = w† Θ_G† D_k w` uses the physical Convention 1 result with
        #     Θ_G†; but `calc_bandreps` in Crystalline.jl (following Cano et al.) computes
        #     characters as `Tr(Θ_G D_k)` (not Θ_G†). To match, we compute `w† Θ_G D_k w`,
        #     achieved by placing `conj(Θ_G)` in the conjugated slot of the dot product.
        Θᴳ_conj = reciprocal_translation_phase(orbital_positions(ptbm), -G) # = conj(Θᴳ)
        # NB: the `sgrep` functor computes `D_k(g) = e^{-2πi(gk)·v} ρ(h)` (physical Conv 1),
        #     but `calc_bandreps` in Crystalline.jl uses the conjugated global phase
        #     `e^{+2πi(gk)·v}` (cf. Crystalline.jl issue #12). To match, we conjugate the
        #     global phase by multiplying by `e^{+4πi(gk)·v}` (flipping the sign of the
        #     exponent).
        v_g = translation(g)
        phase_correction = cispi(4dot(gk, v_g))
        ρ = phase_correction * sgrep(k)
        for (n, v) in enumerate(eachcol(vs))
            v_kpG = Θᴳ_conj * v
            symeigs[j, n] = dot(v_kpG, ρ, v) # = v† Θᴳ conj(D_k) v
            # TODO: preallocate and `mul!` the `Θᴳ_conj * v` term to avoid allocations
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