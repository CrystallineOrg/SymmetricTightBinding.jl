# Note [⚠️ phase]: Crystalline.jl's irreps assume Bloch states `e^{-ik·r}u_k`, opposite to
#   our convention; `flip_bloch_phase` converts them before comparing with the symmetry
#   eigenvalues from `symmetry_eigenvalues`. See issue #137 and
#   `docs/src/devdocs/symmetry_eigenvalue_conventions.md`.

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

julia> brs = bandreps(221);

julia> cbr = @composite brs[1] + brs[2]
40-irrep CompositeBandRep{3} (spinless):
 (3d|A₁g) + (3d|A₁ᵤ) (6 bands)

julia> tbm = tb_hamiltonian(cbr); # a 4-term, 6-band model

julia> ptbm = tbm([1.0, 0.1, -1.0, 0.1]); # fix free coefficients

julia> collect_compatible(ptbm)
2-element Vector{SymmetryVector{3, LGIrrep{3}}}:
 [Γ₁⁻+Γ₃⁻, R₄⁺, M₅⁺+M₁⁻, X₃⁺+X₁⁻+X₂⁻] (3 bands)
 [Γ₁⁺+Γ₃⁺, R₄⁻, M₁⁺+M₅⁻, X₁⁺+X₂⁺+X₃⁻] (3 bands)
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
    # NB: `timereversal = false`, so that we evaluate at -k where the unconverted tables of
    #     `cbr.brs` (used below) apply; this works with and without TR
    clgirsv = flip_bloch_phase(clgirsv; timereversal = false)
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

!!! warning "⚠️ Bloch phase convention"
    The symmetry eigenvalues are computed in this package's Bloch convention,
    ``e^{+i𝐤·𝐫}u_𝐤(𝐫)``. Crystalline.jl's irreps, as returned by e.g. `lgirreps`, assume
    the opposite convention, and must be converted before they are compared (see
    `docs/src/devdocs/symmetry_eigenvalue_conventions.md`).
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
            symeigs[j, n] = dot(v_kpG, D_k, v)  # Convention 1: (Θ_G w)† D_k w
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

!!! warning
    Without time-reversal symmetry, the irrep labels at a **k**-point that is not
    time-reversal invariant may be named after a different **k**-point (e.g., irreps
    `KA₁, KA₂, …` at the **k**-point `K` in plane group *p*3, or `H₁, H₂, …` at `K` in space
    group 143). The labels still correctly describe the bands at the annotated **k**-point.
    See https://github.com/CrystallineOrg/SymmetricTightBinding.jl/issues/137.
"""
function Crystalline.collect_irrep_annotations(ptbm::ParameterizedTightBindingModel; kws...)
    cbr = ptbm.tbm.cbr
    clgirsv = irreps(cbr) # irreps associated to the EBRs (conventional setting)
    clgirsv = flip_bloch_phase(clgirsv; timereversal = first(cbr.brs).timereversal)
    lgirsv = primitivize.(clgirsv) # convert associated groups & irreps to primitive setting
    # NB: use the k-label of the group, not of the irreps: after `flip_bloch_phase`, irreps
    #     named e.g. K₁ may sit at the point KA, and the labels must go where the bands are
    return Dict(map(lgirsv) do lgirs
        symeigs = eachcol(symmetry_eigenvalues(ptbm, group(lgirs)))
        klabel(group(lgirs)) => collect_irrep_annotations(symeigs, lgirs; kws...)
    end)
end

"""
    flip_bloch_phase(
        lgirsv::AbstractVector{<:Collection{<:AbstractLGIrrep{D}}};
        timereversal::Bool
    ) --> Vector{<:Collection{<:AbstractLGIrrep{D}}}

Convert the little group irreps `lgirsv`, over one or more special **k**-points, as
tabulated by Crystalline.jl, from its Bloch phase convention, ``e^{-i𝐤·𝐫}u_𝐤(𝐫)``, to
that of this package, ``e^{+i𝐤·𝐫}u_𝐤(𝐫)``, keeping their labels. `timereversal` indicates
whether time-reversal symmetry is present.

The irreps `lgirs` that Crystalline.jl tabulates at `k` describe our Bloch states at `-k`.
The converted irreps `lgirs′` are placed at `q = position(group(lgirs′))`:
1. If `g∘k ≡ -k` for some operation `g` of the space group (including `g = 1`, if `k ≡ -k`):
   `q = k`, and the irreps are `D′(g⁻¹hg) = D(h)`.
2. Otherwise, with time-reversal: `q = k`, and the irreps are `D′(h) = D(h)*`.
3. Otherwise, without time-reversal: `q = -k`, and the irreps are unchanged. The **k**-label
   is then that of the tabulated **k**-point whose star contains `-k` (e.g., `KA` for `K` in
   plane group *p*3), taken from among those of `lgirsv` if possible.

In all cases, the `i`th operation of the converted group corresponds to the `i`th operation
of `group(lgirs)`.

!!! warning
    The converted irreps must not be passed to Crystalline.jl functions that involve the
    translation phase (e.g., `israyrep`, `realify`, `remap_to_kstar`, or the functor
    `lgir(αβγ)` at non-special **k**-points); character-based functions like
    `find_multiplicities` are safe.
"""
function flip_bloch_phase(
    lgirsv::AbstractVector{<:Collection{IR}};
    timereversal::Bool
) where {D, IR<:AbstractLGIrrep{D}}
    sgnum = num(group(first(lgirsv)))
    all(lgirs -> num(group(lgirs)) == sgnum, lgirsv) ||
        error("all irreps must belong to the same space group")
    cntr = centering(sgnum, D)
    sgops = reduce_ops(spacegroup(sgnum, Val(D); spinful = Val(isspinful(IR))), cntr)
    instar(kv, kv′) = any(g -> isapprox(g * kv, kv′, cntr, #=modw=# true), sgops)
    lgs = nothing # all of Crystalline.jl's little groups of `sgnum`; loaded only if needed

    return map(lgirsv) do lgirs
        lg = group(lgirs)
        kv = position(lg)
        isspecial(kv) || error("only special k-points are supported")
        idx = findfirst(g -> isapprox(g * kv, -kv, cntr, #=modw=# true), sgops)
        isconj = false
        if !isnothing(idx) # case 1: our states at `k`
            g = sgops[idx]
            q, klab = kv, klabel(lg)
            ops′ = [compose(inv(g), compose(h, g, false), false) for h in lg] # h′ = g⁻¹hg
        elseif timereversal # case 2: our states at `k` are the TR partners of those at `-k`
            q, klab, ops′ = kv, klabel(lg), operations(lg)
            isconj = true
        else                # case 3: our states at `-k`
            q, ops′ = -kv, operations(lg)
            i = findfirst(lgirs′ -> instar(position(group(lgirs′)), q), lgirsv)
            klab = if !isnothing(i)
                klabel(group(lgirsv[i]))
            else
                isnothing(lgs) && (lgs = littlegroups(sgnum, Val(D)))
                @something(findfirst(lg′ -> isspecial(position(lg′)) &&
                                            instar(position(lg′), q), lgs),
                           error(lazy"no tabulated k-point has $q in its star"))
            end
        end
        lg′ = typeof(lg)(sgnum, q, klab, ops′)
        Collection([IR(lgir.cdml, lg′, isconj ? conj.(lgir()) : lgir(), nothing,
                       lgir.reality, lgir.iscorep) for lgir in lgirs])
    end
end
