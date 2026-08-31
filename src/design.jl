# Objective-driven model *design*, as opposed to model *fitting*.
#
# `fit` (Optim.jl extension) answers "which hopping amplitudes reproduce this reference
# spectrum?". Often, though, there is no reference spectrum — only a property that the band
# structure ought to have; e.g., "make the Γ₄⁺ triplet an energy-isolated three-fold
# degeneracy". Such questions are optimization problems too: they merely swap the
# least-squares loss for a loss that measures the desired property.
#
# This file provides the first such design objective, `IrrepIsolationObjective`, together
# with the ingredients that any similar objective will want: a uniform k-mesh
# (`uniform_kmesh`), symmetry-based location of the band multiplet carrying a given irrep
# (`IrrepTarget`, `locate_multiplet`), and a verification routine (`isolation_report`).
# Objectives here follow the `fgh!(F, G, H, cs)` calling convention of Optim.jl's
# `only_fgh!`, so that `make_fit_objective`/`multistart_fit` — the moment-seeded,
# basin-hopping search engine already used by `fit` — can drive them unchanged; the driver
# `isolate_irrep` does exactly that (implemented in the Optim.jl extension).

# ---------------------------------------------------------------------------------------- #
# k-mesh

"""
    uniform_kmesh(::Val{D}, Nk::Integer)  -->  Vector{SVector{D, Float64}}

A uniform ``N_k^D`` **k**-mesh over the reduced Brillouin zone ``[-½, ½)^D``, in the basis
of the primitive reciprocal lattice vectors.

Unlike the cell-centered mesh used for density-of-states integration (cf.
`densityofstates`), the mesh includes the zone center Γ: band extrema, and the
symmetry-enforced degeneracies that can make or break a design objective, overwhelmingly sit
at high-symmetry points, and a mesh that misses them can report a gap that is symmetry-
forbidden.

!!! tip "Choose `Nk` divisible by 6"
    Which high-symmetry points land on the mesh depends on `Nk`: an even `Nk` picks up the
    ½-type points (X, M, R, …) and an `Nk` divisible by 3 the ⅓-type ones (K, …). An `Nk`
    divisible by 6 catches both — `Nk = 8`, by contrast, misses K entirely, and a hexagonal
    model sampled on it can look isolated when it is not.
"""
function uniform_kmesh(::Val{D}, Nk::Integer) where D
    Nk ≥ 1 || error("`Nk` must be positive")
    axis = [-0.5 + (j - 1) / Nk for j in 1:Nk]
    return [SVector{D, Float64}(pt) for pt in Iterators.product(ntuple(_ -> axis, Val(D))...)] |> vec
end

# ---------------------------------------------------------------------------------------- #
# irrep targets

"""
    IrrepTarget{D}

A band multiplet, specified by the irrep it transforms as, that a design objective aims at.

Constructed as `IrrepTarget(tbm, irlab)` from a [`TightBindingModel`](@ref) `tbm` and an
irrep label `irlab` (e.g., `"Γ₄⁺"`), which must be among the irreps of the band
representation underlying `tbm` (cf. `irreps(CompositeBandRep(tbm))`).

## Fields
- `k`: the (primitive-setting) **k**-point of the irrep
- `χ`: the irrep's characters, over the operations of its little group
- `d`: the irrep's dimension, i.e., the number of bands in the multiplet
- `Ms`: per-operation matrices ``M^g`` converting eigenvectors of ``H(\\mathbf{k})`` into
  symmetry eigenvalues via ``\\chi_n(g) = (v_n^\\dagger M^g v_n)^*`` (cf.
  `symmetry_operator_matrices`); cached, since they depend only on `k` and on the model's
  band representation, not on its hopping amplitudes
"""
struct IrrepTarget{D}
    irlab :: String
    k     :: SVector{D, Float64}
    χ     :: Vector{ComplexF64}
    d     :: Int
    Ms    :: Vector{Matrix{ComplexF64}}
end

function IrrepTarget(tbm::TightBindingModel{D}, irlab::AbstractString) where D
    cbr = CompositeBandRep(tbm)

    # the irrep must both exist at one of the `cbr`'s k-points & actually feature in its band
    # content; otherwise, no multiplet could ever be located
    sv = SymmetryVector(cbr)
    labs = irreplabels(sv)
    j = findfirst(==(irlab), labs)
    if isnothing(j)
        error("irrep \"$irlab\" is not among the irreps of `tbm` (candidates: \
               $(join(labs, ", ")))")
    elseif iszero(collect(Iterators.flatten(multiplicities(sv)))[j])
        error("irrep \"$irlab\" does not feature in the band content of `tbm` ($sv)")
    end

    lgirsv = primitivize.(irreps(cbr)) # NB: `symmetry_eigenvalues` needs a primitive setting
    i = something(findfirst(lgirs -> any(lgir -> label(lgir) == irlab, lgirs), lgirsv))
    lgirs = lgirsv[i]
    lgir = lgirs[something(findfirst(lgir -> label(lgir) == irlab, lgirs))]

    kv = position(group(lgir))
    isspecial(kv) || error("the k-point of irrep \"$irlab\" has free parameters")
    k = SVector{D, Float64}(constant(kv))

    ops = operations(group(lgir))
    positions = orbital_positions(tbm)
    sgreps = [sgrep_induced_by_siteir(cbr, op, positions) for op in ops]
    Ms = symmetry_operator_matrices(sgreps, positions, k)

    return IrrepTarget{D}(String(irlab), k, characters(lgir), irdim(lgir), Ms)
end

"""
    locate_multiplet(target::IrrepTarget, Es, vs; degen_rtol, degen_atol, atol)
                                                       -->  Union{Nothing, UnitRange{Int}}

Return the band indices of the multiplet transforming as `target`'s irrep, given the
eigenvalues `Es` and eigenvectors `vs` of ``H(\\mathbf{k})`` at `target.k` (in the
Convention 1 basis, i.e., as returned by `solve(…; bloch_phase = Val(false))` or by
diagonalizing a [`TightBindingCache`](@ref) assembly).

Bands are grouped by degeneracy; a group of the right size whose summed symmetry eigenvalues
match the irrep's characters (to `atol`) is returned; if the irrep features more than once
in the model's band content, the lowest-energy such multiplet is returned. `nothing` is
returned if no group matches — which, since the degeneracy is symmetry-enforced, happens only if the multiplet is
*accidentally* degenerate with further bands, so that its energy is shared with, and cannot
be assigned away from, other irreps.
"""
function locate_multiplet(
    target::IrrepTarget,
    Es::AbstractVector{<:Real},
    vs::AbstractMatrix{<:Complex};
    degen_rtol::Real = 1e-8,
    degen_atol::Real = 1e-8,
    atol::Real = 1e-4,
)
    for bands in _degenerate_band_groups(Es, degen_rtol, degen_atol)
        length(bands) == target.d || continue
        matches = true
        for (j, M) in enumerate(target.Ms)
            χ = zero(ComplexF64)
            for n in bands
                v = @view vs[:, n]
                χ += conj(dot(v, M, v)) # [⚠️ phase]: cf. `symmetry_eigenvalues`
            end
            if abs(χ - target.χ[j]) > atol
                matches = false
                break
            end
        end
        matches && return bands
    end
    return nothing
end

# ---------------------------------------------------------------------------------------- #
# isolation objective

"""
    IrrepIsolationObjective(tbm::TightBindingModel{D}, irlab::AbstractString; kws...)

A design objective whose minimization drives the band multiplet transforming as the irrep
`irlab` (e.g., `"Γ₄⁺"`) toward being *energy-isolated*: the constant-energy surface at the
multiplet's energy ``E_t`` must intersect the bands nowhere but at the degeneracy itself.

Used as a functor following Optim.jl's `fgh!` convention, `obj(F, G, H, cs)`: the loss is
returned if `F` is non-`nothing` and its gradient written into `G` if `G` is non-`nothing`
(no Hessian is provided, so `H` must be `nothing`; use a first-order optimizer). See
`make_fit_objective` and `multistart_fit` for how to hand it to Optim.jl, or `isolate_irrep`
for the assembled driver.

## Isolation, as a constraint on sorted bands

Let the multiplet occupy the (sorted) band indices `m:m+d-1` at the irrep's **k**-point
``\\mathbf{k}_0``, with energy ``E_t``. Its own `d` bands are the ones that usually spoil
isolation: they emanate from ``E_t`` at ``\\mathbf{k}_0`` and, left unconstrained, wander
back across it elsewhere in the zone. Isolation therefore demands that they *split* — `p` of
them dropping below ``E_t`` and `d-p` rising above it — and then stay put. Equivalently,
there must exist a filling ``\\nu = m+p-1`` such that

``E_\\nu(\\mathbf{k}) < E_t < E_{\\nu+1}(\\mathbf{k})``    for every
``\\mathbf{k} \\neq \\mathbf{k}_0``.

A split must be proper, `1 ≤ p ≤ d-1`: `p = 0` or `p = d` would put ``E_t`` at a band edge
of the multiplet, where no constant-energy surface separates it from the multiplet's own
bands.

Since the multiplet's bands approach ``E_t`` as ``\\mathbf{k} \\to \\mathbf{k}_0``, that
requirement provably fails near the degeneracy, and **k**-points within `kexclude` of
``\\mathbf{k}_0`` are exempted from it. They are *not* dropped altogether, though: the bands
*outside* the multiplet stay a finite distance from ``E_t`` right up to the degeneracy, and
pushing them away there is exactly what keeps the node clean. Every sampled **k**-point
therefore also carries

``E_{m-1}(\\mathbf{k}) < E_t < E_{m+d}(\\mathbf{k})``,

i.e., "no more than the multiplet's own `d` bands may reach ``E_t``". Sorting takes care of
the remaining bands, so each **k**-point contributes at most four constraints. `kexclude` is
thus a genuine modelling choice: it sets the scale below which the *node* may be approached,
and hence what the attained gap means. (One may equivalently sample `ks` along
high-symmetry lines that keep away from ``\\mathbf{k}_0``, e.g., a decimated `irrfbz_path`
from Brillouin.jl — but note that the outside-band constraints then go unsampled near
``\\mathbf{k}_0``.)

## Loss

Writing ``\\delta_j`` for the signed gaps of the constraints above (positive when
satisfied), the loss maximizes the smallest of them, measured relative to the spectrum's
overall energy scale ``\\sigma`` (its root-mean-square spread over `ks`, so that the loss is
invariant under an overall rescaling of the hopping amplitudes),

``L = -\\text{softmin}_\\beta(\\delta_j/\\sigma)
    = \\beta^{-1}\\log \\sum_j e^{-\\beta\\delta_j/\\sigma}``,

a smooth (log-sum-exp) surrogate for ``-\\min_j \\delta_j/\\sigma`` that concentrates on the
near-touching bands, as `β` sets how sharply. The loss is minimized over the admissible
splits `p` as well as over the amplitudes, unless `split` fixes `p`.

Two small quadratic penalties, ``\\lambda_\\sigma\\log^2(\\sigma/\\sigma_0)`` and
``\\lambda_\\text{shift}\\bar{E}^2`` (with ``\\bar{E}`` the mean band energy), fix the two
gauge directions — overall scale and overall energy offset — that the loss is otherwise
blind to; neither affects which models are optimal, only which representative of a
scale/offset family is returned.

A `lasso` weight ``\\lambda_1`` adds a further ``\\ell_1`` penalty,
``\\lambda_1 N_c^{-1}\\sigma^{-1}\\sum_i|c_i|``, encouraging *sparse* models — ones that
achieve isolation with as few hopping terms as possible. It is normalized by the number of
terms ``N_c`` and by ``\\sigma``, so that it is scale-invariant like the rest of the loss
(an unnormalized ``\\ell_1`` term would be minimized by simply shrinking every amplitude)
and so that a given ``\\lambda_1`` means the same thing across models of different size.
Expect it to buy sparsity at the cost of gap: verify the outcome with
[`isolation_report`](@ref).

If the multiplet cannot be located (cf. [`locate_multiplet`](@ref)), the loss takes the
large sentinel value `FAILED_ISOLATION_LOSS` with vanishing gradient, so that the
surrounding multi-start search simply discards the trial.

## Keyword arguments
- `ks` (default, `uniform_kmesh(Val(D), Nk)`): the **k**-points over which isolation is
  imposed. The irrep's own **k**-point is always prepended, whether or not it is in `ks`.
- `Nk` (default, `6`): linear mesh density of the default **k**-mesh.
- `kexclude` (default, `0.1`): radius of the neighbourhood of ``\\mathbf{k}_0`` — in
  reduced coordinates, minimal image — within which the multiplet is not required to have
  split, as described above. Must be positive.
- `split` (default, `nothing`): if set to an integer `p ∈ 1:d-1`, fixes the number of
  multiplet bands required to fall *below* ``E_t``; otherwise every admissible `p` is tried
  at each evaluation and the best is used. (The two ends, `p` and `d-p`, are related by
  ``H \\to -H`` whenever the model admits it, so fixing `split` mostly saves work.)
- `β` (default, `20.0`): sharpness of the log-sum-exp soft-minimum; larger values track the
  true minimum gap more closely, at the cost of a less smooth loss landscape.
- `λ_scale` (default, `0.1`) and `σ₀` (default, `1.0`): weight and target of the scale-gauge
  penalty.
- `λ_shift` (default, `0.1`): weight of the energy-offset-gauge penalty.
- `lasso` (default, `nothing`): if set to a positive number, the weight ``\\lambda_1`` of
  the sparsity-encouraging ``\\ell_1`` penalty above; `nothing` disables it.

## Fields (diagnostics)

After a call, `obj.bands` holds the located multiplet's band indices (empty if it could not
be located), `obj.split` the `p` it settled on, and `obj.relgap` the *hard* (non-smoothed)
minimum gap ``\\min_j\\delta_j/\\sigma`` at the evaluated amplitudes — positive exactly when
the multiplet is isolated over the constrained **k**-points.
"""
mutable struct IrrepIsolationObjective{D, K}
    const cache   :: TightBindingCache{D, HERMITIAN, Float64, K}
    const target  :: IrrepTarget{D}
    const split_κ :: Vector{Bool}           # k-points outside the exclusion zone (cf. below)
    const splits  :: UnitRange{Int}         # admissible #(multiplet bands below Et)
    const β       :: Float64
    const λ_scale :: Float64
    const σ₀      :: Float64
    const λ_shift :: Float64
    const lasso   :: Float64                # ℓ₁ weight; 0 ⇔ disabled
    # cached spectral-moment data: `a[i] = ∑ₖ tr hᵢ(k)` and `Q̄[i,j] = ∑ₖ tr[hᵢ(k)hⱼ(k)]`,
    # giving the mean band energy Ē = aᵀc/Nᵏᵇ and mean square ⟨E²⟩ = cᵀQ̄c/Nᵏᵇ exactly and
    # without diagonalization (cf. `SpectralMoments` in the Optim.jl extension)
    const a       :: Vector{Float64}
    const Q̄       :: Matrix{Float64}
    const Nᵏᵇ     :: Int                    # #k-points × #bands
    # work arrays: eigenvalues & eigenvectors across the mesh, retained between the loss and
    # gradient passes of a single evaluation, plus the constraint list (cf. `_constraints!`)
    const Esm     :: Matrix{Float64}        # Esm[n, κ]    = Eₙ(ks[κ])
    const vsm     :: Array{ComplexF64, 3}   # vsm[:, n, κ] = vₙ(ks[κ])
    const cs_buf  :: Vector{Tuple{Int, Int, Float64}}
    # diagnostics from the most recent evaluation
    bands         :: UnitRange{Int}
    split         :: Int
    relgap        :: Float64
end

"Sentinel loss returned when the target multiplet cannot be located; cf. `locate_multiplet`."
const FAILED_ISOLATION_LOSS = 1e3

function IrrepIsolationObjective(
    tbm::TightBindingModel{D, HERMITIAN},
    irlab::AbstractString;
    Nk::Integer = 6,
    ks::AbstractVector{<:ReciprocalPointLike{D}} = uniform_kmesh(Val(D), Nk),
    kexclude::Real = 0.1,
    split::Union{Nothing, Integer} = nothing,
    β::Real = 20.0,
    λ_scale::Real = 0.1,
    σ₀::Real = 1.0,
    λ_shift::Real = 0.1,
    lasso::Union{Nothing, Real} = nothing,
) where D
    target = IrrepTarget(tbm, irlab)
    d, N = target.d, tbm.N
    d < N || error(
        "the multiplet of \"$irlab\" spans all $N bands of `tbm`: there are no other bands \
         for it to be isolated from")
    d > 1 || error(
        "\"$irlab\" is one-dimensional: a single band cannot be split across $(irlab)'s \
         energy, so isolation is not a meaningful goal here")
    kexclude > 0 || error("`kexclude` must be positive (it must at least exclude `k₀`)")
    splits = if isnothing(split)
        1:d-1
    else
        1 ≤ split ≤ d-1 || error("`split` must lie in `1:$(d-1)` for the $d-dimensional \
                                  irrep \"$irlab\"")
        split:split
    end

    # the target k-point must be sampled; prepend it, so that it is at cache index 1
    ks′ = vcat([SVector{D, Float64}(target.k)], [SVector{D, Float64}(k) for k in ks])
    cache = TightBindingCache(tbm, ks′)
    split_κ = [_reduced_distance(k, target.k) ≥ kexclude for k in ks′]
    any(split_κ) || error("every k-point of `ks` lies within `kexclude` of the irrep's \
                           k-point: the multiplet's own splitting would go unconstrained")

    Nᶜ, Nᵏ = length(tbm), length(ks′)
    hs = cache.hs
    a = [sum(hsₖ -> real(tr(hsₖ[i])), hs) for i in 1:Nᶜ]
    Q̄ = [sum(hsₖ -> real(dot(hsₖ[i], hsₖ[j])), hs) for i in 1:Nᶜ, j in 1:Nᶜ]

    return IrrepIsolationObjective{D, typeof(ks′)}(
        cache, target, split_κ, splits, β, λ_scale, σ₀, λ_shift, something(lasso, 0.0),
        a, Q̄, Nᵏ * N,
        Matrix{Float64}(undef, N, Nᵏ), Array{ComplexF64, 3}(undef, N, N, Nᵏ),
        Vector{Tuple{Int, Int, Float64}}(undef, 4Nᵏ),
        1:0, 0, NaN)
end

# distance between two k-points in reduced coordinates, minimized over reciprocal lattice
# translations. A proxy for the true (metric) distance — the reciprocal basis is not part of
# a `TightBindingModel` — and hence anisotropic for non-cubic lattices; adequate for merely
# fencing off a neighbourhood of the degeneracy.
_reduced_distance(k, k₀) = sqrt(sum(δ -> (δ - round(δ))^2, k .- k₀))

# The constraint list realizing isolation at a given multiplet `bands` and split `p`, written
# into `o.cs_buf` as `(κ, n, s)` triples denoting a signed gap `δ = s(Eₙ(ks[κ]) - Et)`, and
# ordered by `κ` so that a gradient pass needs one Feynman–Hellmann call per k-point.
# At every sampled k-point, including inside the excluded neighbourhood of the degeneracy:
#   - band m-1 must stay below Et and band m+d above it, i.e., no more than the `d` bands of
#     the multiplet may reach Et. Sorting makes the remaining outside bands automatic — and
#     these constraints are perfectly meaningful next to the degeneracy, where the bands
#     outside the multiplet are separated from Et by a finite amount.
# At every k-point *outside* that neighbourhood, additionally:
#   - the valence band ν = m+p-1 must stay below Et and the conduction band ν+1 above it,
#     i.e., the multiplet's own bands must have split. (These are the constraints that the
#     degeneracy necessarily violates as k → k₀, and hence the only ones exempted there.)
function _constraints!(o::IrrepIsolationObjective, bands::UnitRange{Int}, p::Integer)
    m, d, N = first(bands), length(bands), o.cache.tbm.N
    ν = m + p - 1
    buf, j = o.cs_buf, 0
    for κ in eachindex(o.cache.ks)
        m > 1     && (buf[j+=1] = (κ, m-1, -1.0))
        m + d ≤ N && (buf[j+=1] = (κ, m+d, +1.0))
        if o.split_κ[κ]
            buf[j+=1] = (κ, ν,   -1.0)
            buf[j+=1] = (κ, ν+1, +1.0)
        end
    end
    return j
end

# the softmin reduction over the `nc` first constraints of `o.cs_buf`: returns the hard
# minimum normalized gap `δ̂min` and the (stabilized) sum `S = ∑ⱼexp[-β(δ̂ⱼ-δ̂min)]`, from
# which the loss is `-δ̂min + β⁻¹log S`
function _softmin(o::IrrepIsolationObjective, nc::Integer, Et::Real, σ::Real)
    δ̂min = Inf
    for j in 1:nc
        (κ, n, s) = o.cs_buf[j]
        δ̂min = min(δ̂min, s * (o.Esm[n, κ] - Et) / σ)
    end
    S = 0.0
    for j in 1:nc
        (κ, n, s) = o.cs_buf[j]
        S += exp(-o.β * (s * (o.Esm[n, κ] - Et) / σ - δ̂min))
    end
    return δ̂min, S
end

function (o::IrrepIsolationObjective)(F, G, H, cs::AbstractVector{<:Real})
    isnothing(H) || error("`IrrepIsolationObjective` provides no Hessian: use a first-order \
                           optimizer")
    isnothing(G) || fill!(G, zero(eltype(G)))
    Nᵏ = length(o.cache.ks)

    # diagonalize across the mesh, retaining eigenvectors for the gradient pass below
    for κ in 1:Nᵏ
        Es, vs = eigen!(o.cache(cs, κ))
        o.Esm[:, κ] .= Es
        o.vsm[:, :, κ] .= vs
    end

    # locate the target multiplet & its energy Et
    bands = locate_multiplet(o.target, (@view o.Esm[:, 1]), (@view o.vsm[:, :, 1]))
    if isnothing(bands)
        o.bands, o.split, o.relgap = 1:0, 0, NaN
        return isnothing(F) ? nothing : FAILED_ISOLATION_LOSS
    end
    o.bands = bands
    Et = sum(@view o.Esm[bands, 1]) / o.target.d

    # energy scale σ & offset Ē, from the cached spectral moments (no diagonalization)
    Ē  = dot(o.a, cs) / o.Nᵏᵇ
    σ² = max(dot(cs, o.Q̄, cs) / o.Nᵏᵇ - Ē^2, eps())
    σ  = sqrt(σ²)

    # pick the split `p` that isolates best, then leave its constraints in `o.cs_buf`
    p, nc, δ̂min, S = 0, 0, -Inf, NaN
    for p′ in o.splits
        nc′ = _constraints!(o, bands, p′)
        δ̂min′, S′ = _softmin(o, nc′, Et, σ)
        # lower loss ⇔ larger -δ̂min + β⁻¹log S; compare on that combination directly
        if p == 0 || -δ̂min′ + log(S′)/o.β < -δ̂min + log(S)/o.β
            p, nc, δ̂min, S = p′, nc′, δ̂min′, S′
        end
    end
    p == last(o.splits) || _constraints!(o, bands, p) # restore the winning constraint list
    o.split, o.relgap = p, δ̂min

    if !isnothing(G)
        # gauge-penalty gradients & the moment-derived ∇σ (∇σ² = 2[Q̄c - Ēa]/Nᵏᵇ)
        ∇σ = (o.Q̄ * cs .- Ē .* o.a) ./ (o.Nᵏᵇ * σ)
        @. G += o.λ_scale * 2log(σ/o.σ₀) * ∇σ / σ + o.λ_shift * 2Ē * o.a / o.Nᵏᵇ
        if !iszero(o.lasso) # ∂/∂c of λ₁‖c‖₁/(Nᶜσ)
            λ̂ = o.lasso / (length(cs) * σ)
            ℓ₁ = sum(abs, cs) # NB: hoisted out of the `@.` below, which would otherwise
            @. G += λ̂ * (sign(cs) - ℓ₁ * ∇σ / σ) # broadcast the reduction element-wise
        end

        # Feynman–Hellmann gradients; ∇Et is the multiplet-averaged gradient (equal across
        # the multiplet, by symmetry, up to numerical noise)
        ∇Es = energy_gradient_wrt_hopping(o.cache, 1, ((@view o.Esm[:,1]), (@view o.vsm[:,:,1])))
        ∇Et = sum(n -> ∇Es[n], bands) ./ o.target.d
        κ_prev = 0
        for j in 1:nc
            (κ, n, s) = o.cs_buf[j]
            if κ ≠ κ_prev # `cs_buf` is κ-ordered: one gradient call per k-point suffices
                ∇Es = energy_gradient_wrt_hopping(
                    o.cache, κ, ((@view o.Esm[:, κ]), (@view o.vsm[:, :, κ])))
                κ_prev = κ
            end
            δ = s * (o.Esm[n, κ] - Et)
            w = exp(-o.β * (δ/σ - δ̂min)) / S   # softmin weight of this constraint
            ∇Eₙ = ∇Es[n]
            # ∂L/∂c = -∑ⱼ wⱼ ∂(δⱼ/σ)/∂c, with ∂(δ/σ) = ∂δ/σ - δ∂σ/σ²
            @inbounds for i in eachindex(G)
                ∂δ = s * (∇Eₙ[i] - ∇Et[i])
                G[i] -= w * (∂δ / σ - δ * ∇σ[i] / σ²)
            end
        end
    end

    isnothing(F) && return nothing
    return (-δ̂min + log(S)/o.β + o.λ_scale * log(σ/o.σ₀)^2 + o.λ_shift * Ē^2 +
            o.lasso * sum(abs, cs) / (length(cs) * σ))
end

# ---------------------------------------------------------------------------------------- #
# verification

"""
    isolation_report(ptbm::ParameterizedTightBindingModel{D}, irlab::AbstractString;
                     Nk = 24, ks = uniform_kmesh(Val(D), Nk), kexclude = 0.1, split = nothing)

Check whether the band multiplet transforming as the irrep `irlab` is energy-isolated in
`ptbm`: whether, at every sampled **k**-point, no band outside the multiplet reaches the
multiplet's energy ``E_t``, and whether, outside the excluded neighbourhood of the irrep's
own **k**-point, the multiplet's `d` own bands have split cleanly into `split` below ``E_t``
and `d-split` above it. See [`IrrepIsolationObjective`](@ref) for the full statement and for
the roles of `kexclude` and `split`; the keywords mean the same thing here, and should be
given the same values that the search used (`kexclude` in particular, since the attained gap
is measured at its boundary). Intended as an independent verification of the outcome of
[`isolate_irrep`](@ref), on a **k**-sampling finer than the one it optimized over.

!!! note "Sample the high-symmetry points"
    A symmetry-enforced degeneracy at a high-symmetry **k**-point can rule out a filling
    outright, and a **k**-sampling that misses that point will happily report isolation that
    does not exist. Prefer an `Nk` divisible by 6 (so that both the ½- and ⅓-type
    high-symmetry points are sampled), or supply `ks` that visit them explicitly.

Returns a `NamedTuple` with fields:
- `isolated`: whether the multiplet is isolated over the sampled **k**-points
- `bands`: the multiplet's band indices (empty if it could not be located, in which case
  `isolated` is `false` and the remaining fields are `NaN`)
- `split`: the number of multiplet bands falling below ``E_t``; the best such split, if
  `split` was not fixed
- `Et`: the multiplet's energy
- `below`, `above`: the closest that any band required to lie below (resp. above) ``E_t``
  comes to it; `Inf` if there are no such bands
- `gap`: the smaller of `below` and `above`; positive exactly when `isolated`
- `σ`: the spectrum's root-mean-square energy spread over the sampled **k**-points
- `relgap`: `gap/σ`, the scale-invariant figure of merit that `IrrepIsolationObjective`
  maximizes
"""
function isolation_report(
    ptbm::ParameterizedTightBindingModel{D, HERMITIAN},
    irlab::AbstractString;
    Nk::Integer = 24,
    ks::AbstractVector{<:ReciprocalPointLike{D}} = uniform_kmesh(Val(D), Nk),
    kexclude::Real = 0.1,
    split::Union{Nothing, Integer} = nothing,
) where D
    target = IrrepTarget(ptbm.tbm, irlab)
    Es₀, vs₀ = solve(ptbm, target.k; bloch_phase = Val(false))
    bands = locate_multiplet(target, Es₀, vs₀)
    if isnothing(bands)
        return (; isolated = false, bands = 1:0, split = 0, Et = NaN, below = NaN,
                  above = NaN, gap = NaN, σ = NaN, relgap = NaN)
    end
    m, d, N = first(bands), target.d, ptbm.tbm.N
    Et = sum(@view Es₀[bands]) / d

    # tabulate the spectrum, flagging which k-points must show the multiplet already split
    Esm, split_κ = [Es₀], [false] # the irrep's own k-point is never split-constrained
    ΣE, ΣE², Nᵏᵇ = sum(Es₀), sum(abs2, Es₀), N
    for k in ks
        Es = spectrum_single_k(ptbm, k)
        ΣE += sum(Es); ΣE² += sum(abs2, Es); Nᵏᵇ += N
        push!(Esm, Es)
        push!(split_κ, _reduced_distance(k, target.k) ≥ kexclude)
    end
    any(split_κ) || error("every k-point of `ks` lies within `kexclude` of the irrep's \
                           k-point: the multiplet's own splitting would go unchecked")
    σ = sqrt(max(ΣE²/Nᵏᵇ - (ΣE/Nᵏᵇ)^2, 0.0))

    # everywhere: no band outside the multiplet may reach Et (bands m-1 and m+d suffice, by
    # sorting); outside the exclusion zone: the multiplet itself must have split at ν = m+p-1
    best = (; split = 0, below = NaN, above = NaN, gap = -Inf)
    for p in (isnothing(split) ? (1:d-1) : (split:split))
        1 ≤ p ≤ d-1 || error("`split` must lie in `1:$(d-1)` for \"$irlab\"")
        ν = m + p - 1
        below, above = Inf, Inf
        for (Es, isplit) in zip(Esm, split_κ)
            m > 1     && (below = min(below, Et - Es[m-1]))
            m + d ≤ N && (above = min(above, Es[m+d] - Et))
            if isplit
                below = min(below, Et - Es[ν])
                above = min(above, Es[ν+1] - Et)
            end
        end
        gap = min(below, above)
        gap > best.gap && (best = (; split = p, below, above, gap))
    end

    return (; isolated = best.gap > 0, bands, best.split, Et, best.below, best.above,
              best.gap, σ, relgap = best.gap/σ)
end
