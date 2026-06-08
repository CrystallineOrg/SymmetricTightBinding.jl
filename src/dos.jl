# Density of states via the (generalized) Gilat–Raubenheimer method.
#
# We mesh the reduced Brillouin zone [-½, ½]^D into a uniform grid of `Nk^D` cubic cells.
# In each cell we linearize every band around the cell center `k₀`,
#     Eₙ(k) ≈ E₀ + vₙ·(k − k₀),   vₙ = ∇ₖEₙ(k₀)   (the group velocity),
# whose exact contribution to the DOS is the per-cell integral
#     C(δE) = ∫_cell δ(δE − v·Δk) d^D(Δk),   δE = E − E₀,
# i.e. the (D−1)-dimensional measure of the iso-energy surface {v·Δk = δE} inside the cell
# divided by |v|. `C` has a closed-form, piecewise-quadratic expression (`_GGRContribution`)
# with breakpoints set by the sorted velocity components (`_breakpoints`).
#
# Because `energy_gradient_wrt_momentum` already returns velocities in *reduced* coordinates
# (∂E/∂kᵢ with kᵢ relative to the primitive reciprocal vector bᵢ), and our mesh is a cube in
# those same reduced coordinates, the velocity frame and the mesh frame coincide: we feed
# reduced velocities directly into the *cubic* GR formula with half-side `h = 1/(2Nk)`. No 
# affine transform and no lattice basis are needed (the affine step in the original
# "generalized" GR is only there to reconcile Cartesian velocities with a reduced-coordinate
# mesh).
#
# The result is the physical *states-per-unit-cell* DOS: ∫g(E)dE equals the number of bands.
# (Working in reduced coordinates, the (2π)^D and unit-cell-volume factors of the physical
# DOS cancel exactly, since V_BZ = (2π)^D / V_uc.)
#
# --- Why "bins"? (this is the standard GGR accumulation, not a new idea) ------------------
# The reference GGR implementation builds an energy *histogram*: it adds each cell's
# contribution into the energy bins that the cell's energy window overlaps. We do exactly
# the same, only with bin edges placed at the midpoints between the requested `energies`
# (so that arbitrary, non-uniform energy grids are supported).
#
# Picture the energy axis. The requested points φⱼ are bin *centers*; the bin edges eⱼ sit
# at the midpoints between neighboring φⱼ. A single band in a single cell contributes the
# function C(δE): a little "tent" (or box, in 1D) of total area (2h)^D, centered on the band
# energy E₀ and supported on [E₀-wmax, E₀+wmax] with wmax = h·Σ|vᵢ|.
#
#       φ₁         φ₂         φ₃         φ₄         φ₅       ← requested energies (centers)
#    ┌─────┬──────────┬──────────┬──────────┬──────────┐
#    │ b₁  │    b₂    │    b₃    │    b₄    │    b₅     │    ← bins (edges eⱼ @ φ-midpoints)
#    e₁    e₂         e₃         e₄         e₅         e₆
#                            C(δE)
#                             ╱▔╲                           ← (cell, band) contribution;
#                            ╱    ╲                            area = (2h)^D, centered at E₀
#    ───────────────────────╱──────╲───────────────────► E
#                       E₀-wmax  E₀  E₀+wmax
#                            └──┬──┘
#         for each overlapped bin j:   mass[j] += ∫_{bin j ∩ support} C(δE) dδE
#
# Summing over all cells and bands and finally taking g[j] = mass[j] / width(bin j) gives
# the DOS as a (bin-averaged) density at each φⱼ.
#
# The key property is *mass conservation*: because the bins tile the axis, a cell's entire
# area (2h)^D is always deposited *somewhere*, even when its support is far narrower than a
# bin (a nearly-flat band) and falls between two φⱼ. A naive "evaluate C(φⱼ-E₀) pointwise"
# scheme instead drops such cells, and its sum rule ∫g dE then collapses toward 0 as `Nk`
# grows (we measured 2.0 → 0.18 in 1D). Binning avoids this: ∫g dE = N_bands at any `Nk`.

"""
    densityofstates(
        ptbm::ParameterizedTightBindingModel{D}, 
        energies; 
        kws...
    ) --> Vector{Float64}

Compute the density of states (DOS) of the Hermitian tight-binding model `ptbm`, evaluated
at each element of `energies`, using the (generalized) Gilat–Raubenheimer method [^1].

The DOS is returned as a `Vector{Float64}` of the same length as `energies`. It is
normalized *per unit cell*, i.e., ``\\int g(E)\\,\\mathrm{d}E = N`` where `N` is the number
of bands. This is the physical states-per-unit-cell DOS and is independent of the lattice
basis. See the *Algorithm* section below for details.

## Arguments
- `ptbm :: ParameterizedTightBindingModel{D, HERMITIAN}`. Currently, only Hermitian models
  are supported.
- `energies`: a real-valued, indexable iterable of energies at which to evaluate the DOS. 
  Must be
  sorted in ascending order and free of duplicates (both are checked); spacing may be
  non-uniform. For two or more energies, the DOS at `energies[j]` is the average density
  over the bin centered on `energies[j]` (bin edges at midpoints between neighboring
  energies); for a single energy, the raw Gilat–Raubenheimer density at that energy is
  returned.

## Keyword arguments
- `Nk::Integer` (default, 50): number of mesh cells per reduced dimension.
- `offset` (default, 0.0): a real scalar or `D`-tuple added to every cell-center coordinate
  (in reduced coordinates). The default places the cell-centers a half-cell off Γ. Manual
  adjustments of `offset` can be used to further avoid high-symmetry points, to mitigate
  band-degeneracy artifacts or points of zero group velocity.
- `transform` (default, `nothing` ~ identity): an optional function mapping band eigenvalues
  to a transformed spectral variable (e.g. `sqrt` to go from `ω²`-like eigenvalues to
  frequencies `ω`). When not `nothing`, both the eigenvalues and band velocities used by the
  GGR method are transformed.

## Algorithm

The DOS is evaluated with the generalized Gilat–Raubenheimer method [^1]. The reduced
Brillouin zone ``[-½, ½]^D`` is meshed into a uniform grid of `Nk^D` cubic cells. Within
each cell, bands are linearized about the cell center ``\\mathbf{k}_0`` using its group
velocity ``\\mathbf{v} = \\nabla_{\\mathbf{k}} E_n(\\mathbf{k}_0)`` (obtained via the
Feynman–Hellmann theorem from [`energy_gradient_wrt_momentum`](@ref)).
The exact contribution of the linearized band to the DOS is the per-cell integral

```math
C(\\delta E) = \\int_{\\text{cell}} \\delta(\\delta E - \\mathbf{v}\\cdot\\Delta\\mathbf{k})
                  \\,\\mathrm{d}^D(\\Delta\\mathbf{k}), \\qquad \\delta E = E - E_0,
```

i.e. the ``(D-1)``-dimensional measure of the iso-energy surface inside the cell divided by
``|\\mathbf{v}|``; it has a closed-form, piecewise-quadratic expression. Summing these
contributions over all cells and bands gives the DOS, converging as `𝒪(1/Nk^2)` and is
considerably smoother than histogram or Gaussian-broadening schemes at equal mesh density.

## Example
```julia
brs = calc_bandreps(17, Val(2));                 # build nearest-neighbor graphene model
cbr = @composite brs[5];                         #   → (2b|A₁) EBR (2 bands)
ptbm = tb_hamiltonian(cbr, [[0,0]])([0.0, 1.0]); #   → unit nearest-neighbor hopping

energies = range(-6.5, 6.5, 501);                # eigenvalues span [-6, 6]

dos = densityofstates(ptbm, energies);           # evaluate DOS across all `energies`

sum(dos) * step(energies)                        # DOS integral ≈ number of bands
2.0000 …
```

[^1]: Liu, Lu, et al., *Generalized Gilat–Raubenheimer method for density-of-states
    calculation in photonic crystals*,
    [J. Opt. **20**, 044005 (2018)](https://doi.org/10.1088/2040-8986/aaae52). The original
    method is due to Gilat & Raubenheimer, Phys. Rev. **144**, 390 (1966).
"""
function densityofstates(
    ptbm::ParameterizedTightBindingModel{D, S},
    energies;
    Nk::Integer = 50,
    offset::Union{Real, Tuple{Vararg{Real, D}}} = ntuple(_ -> 0.0, Val(D)),
    transform::F = nothing,
) where {D, S, F}
    if S !== HERMITIAN
        error("`densityofstates` is only implemented for HERMITIAN models (got $S): a \
               non-Hermitian spectrum is complex and the real-energy DOS is ill-defined")
    end
    if D ∉ (1, 2, 3)
        error(lazy"`densityofstates` is only implemented for D ∈ {1,2,3}; got D = $D")
    end
    Nk ≥ 1 || error(lazy"`Nk` must be a positive integer; got $Nk")

    h = 1 / (2Nk)         # half cell side length (reduced coordinates)
    cellmass = (2h)^D     # ∫C dδE over a single cell = cell volume (deposited per band/cell)

    M = length(energies)
    M == 0 && return Float64[]
    issorted(energies) || error("`energies` must be sorted in ascending order")
    allunique(energies) || error("`energies` must not contain duplicate values")
    E₁ = first(energies)                              # only used on the M == 1 (single E) path
    edges = M == 1 ? Float64[] : _bin_edges(energies) # length M+1; bin `j` is [edges[j],edges[j+1]]
    accum = zeros(Float64, M)                         # raw GR density (M==1) or per-bin mass (M≥2)

    # sweep the GR mesh, depositing every (cell, band)'s contribution into `accum`
    Nᵇ = ptbm.tbm.N
    ∇Hs = ntuple(_ -> Matrix{ComplexF64}(undef, Nᵇ, Nᵇ), Val(D)) # reused scratch across the mesh
    for (k, weight) in _dos_kmesh(Val(D), Nk, offset)
        Es, us = solve(ptbm, k; bloch_phase = Val(false))
        vs = energy_gradient_wrt_momentum(ptbm, k, (Es, us); ∇Hs) # group velocities
        for n in 1:Nᵇ
            # remap the eigenvalues & velocity to the transformed variable φₙ = transform(Eₙ),
            # rescaling the velocity by the chain rule, vφ = transform'(Eₙ)·v
            Eₙ = Es[n]
            if transform === nothing
                φₙ = Eₙ
                vφ = Tuple(vs[n])
            else
                φₙ = Float64(transform(Eₙ)::Real)         # transform may return any `Real`
                dfac = _fd_derivative(transform, Eₙ)       # transform'(Eₙ)
                vφ = ntuple(i -> Float64(dfac * vs[n][i]), Val(D))
            end

            if M == 1
                # single energy: no bins to integrate over — accumulate the raw GR density
                # (the sum of per-cell contributions `C(E - φₙ)`, the linearized DOS at `E`)
                accum[1] += weight * _GGRContribution(vφ, h)(E₁ - φₙ)
            else
                _deposit!(accum, edges, weight, φₙ, vφ, h, cellmass)
            end
        end
    end

    M == 1 && return accum # single energy: raw GR density at that energy
    return @inbounds [accum[j] / (edges[j+1] - edges[j]) for j in 1:M]
end

# Deposit a single (cell, band) contribution `C(δE)` of total mass `cellmass` into the
# energy-bin masses `mass` (bins delimited by `edges`), weighted by `weight`. The band sits
# at abscissa `φₙ` with (reduced-coordinate) velocity `vφ`; `h` is the cell half-side.
function _deposit!(
    mass::Vector{Float64}, edges::Vector{Float64}, weight::Float64, φₙ::Float64,
    vφ::NTuple{D, Float64}, h::Float64, cellmass::Float64,
) where {D}
    M = length(mass)
    C = _GGRContribution(vφ, h) # precompute once (sorted speeds+breakpoints); reused below
    wmax = last(C.w)            # support half-width (the largest breakpoint)
    if iszero(wmax) # perfectly flat cell: deposit its whole mass into φₙ's bin
        j = clamp(searchsortedlast(edges, φₙ), 1, M)
        mass[j] += weight * cellmass
    else
        lo, hi = φₙ - wmax, φₙ + wmax
        jlo = clamp(searchsortedlast(edges, lo), 1, M)
        jhi = clamp(searchsortedlast(edges, hi), 1, M)
        for j in jlo:jhi
            a = max(edges[j], lo) - φₙ    # integrate C over (bin ∩ support), shifted
            b = min(edges[j+1], hi) - φₙ  # to the δE = E − φₙ variable
            b > a && (mass[j] += weight * _C_integral(a, b, C))
        end
    end
    return mass
end

# central finite-difference derivative of a scalar function (for the `transform` chain rule)
function _fd_derivative(f::F, x::Float64) where {F}
    ε = cbrt(eps(Float64)) * max(one(x), abs(x))
    return (f(x + ε) - f(x - ε)) / (2ε)
end

# ---------------------------------------------------------------------------------------- #
# k-mesh: cell centers of a uniform `Nk^D` grid over the reduced BZ [-½,½]^D, with a uniform
# weight of 1.0 each. The `(k, weight)` interface is the seam for a future symmetry-reduced
# mesh (which would emit fewer representatives with integer/fractional multiplicities).

function _dos_kmesh(Dⱽ::Val{D}, Nk::Integer, offset) where {D}
    off = offset isa Real ? ntuple(_ -> Float64(offset), Dⱽ) :
                            NTuple{D, Float64}(offset)
    # per-dimension cell centers: cell `j` (j = 1…Nk) is [cⱼ-h, cⱼ+h] with cⱼ = -½+(j-½)/Nk,
    # so the cells exactly tile [-½, ½]; the +(j-½)/Nk centering already offsets Γ by h
    centers = ntuple(d -> [(-0.5 + (j - 0.5) / Nk + off[d]) for j in 1:Nk], Dⱽ)
    grid = Iterators.product(centers...)
    return ((SVector{D, Float64}(pt), 1.0) for pt in grid)
end

# bin edges placed at midpoints between the (sorted) energies; the 2 end bins are reflected
# to half the adjacent spacing so the bins tile the axis and stay centered on each energy
function _bin_edges(energies)
    M = length(energies)
    edges = Vector{Float64}(undef, M + 1)
    @inbounds for j in 2:M
        edges[j] = (energies[j-1] + energies[j]) * 0.5
    end
    @inbounds edges[1] = 2 * energies[1] - edges[2]
    @inbounds edges[M+1] = 2 * energies[M] - edges[M]
    return edges
end

# ---------------------------------------------------------------------------------------- #
# Per-cell GR contribution C(δE): the analytic value of ∫_cell δ(δE − v·Δk) d^D(Δk) for a
# cubic cell of half-side `h`, i.e. the (D-1)-volume of the iso-energy slice divided by |v|.
# `C` is symmetric in δE → −δE. A perfectly flat cell yields `Inf` at δE = 0 (the integrable
# 1/|v| singularity) and 0 otherwise - callers handle flat cells before evaluating `C`,
# except on the single-energy density path.
#
# The descending-abs-sorted velocity magnitudes, the positive breakpoints `w` (the energy
# half-widths at which the analytic form of `C` changes), and `|v|²` are all invariants of
# a (cell, band). They are precomputed once into a `_GGRContribution` kernel, which is then
# evaluated as a functor `C(δE)` at each Simpson node and reused across every energy bin,
# avoiding repeating `abs`/`sort`/breakpoint work for efficiency.
#
# The 1D and 2D forms are derived (and mass-checked: ∫C dδE = (2h)^D) directly; the 3D form
# follows the generalized GR reference implementation (Liu et al.; `DOS_GGR.m`).

# velocity magnitudes |vᵢ| in descending order (largest first); the contributions and
# breakpoints below all assume this ordering
_descending_abs_sort(v::NTuple) = sort(abs.(v); rev=true)

# positive breakpoints of `C(δE)`, for velocity magnitudes `v` already sorted descending
# (v1 ≥ v2 ≥ … ≥ 0)
_breakpoints(v::NTuple{1, Float64}, h) = (v[1] * h,)
_breakpoints(v::NTuple{2, Float64}, h) = ((v[1] - v[2]) * h, (v[1] + v[2]) * h)
function _breakpoints(v::NTuple{3, Float64}, h)
    v1, v2, v3 = v
    return (h * abs(v1 - v2 - v3), h * (v1 - v2 + v3), h * (v1 + v2 - v3), h * (v1 + v2 + v3))
end

# Prebuilt per-cell contribution kernel: sorted speeds `v`, breakpoints `w` (Nw = 1,2,4 for
# D = 1,2,3), and `|v|²` (3D only), are computed once from a raw velocity. Evaluate `C(δE)`
# via the functor methods below. The breakpoints are ascending, so `last(w)` is the support
# half-width (`C(δE) = 0` for `|δE| > last(w)`).
struct _GGRContribution{D, Nw}
    v  :: NTuple{D, Float64}   # |vᵢ| sorted descending: v1 ≥ v2 ≥ … ≥ 0
    w  :: NTuple{Nw, Float64}  # positive breakpoints, ascending
    V2 :: Float64              # |v|² (only used for D = 3)
    h  :: Float64
end
function _GGRContribution(v_raw::NTuple{D, Float64}, h::Float64) where {D}
    v = _descending_abs_sort(v_raw) # v1 ≥ v2 ≥ … ≥ 0; formulas below assume this ordering
    w = _breakpoints(v, h)
    V2 = D == 3 ? dot(v, v) : 0.0 # |v|² = Σ vᵢ²
    return _GGRContribution(v, w, V2, h)
end

# 1D: point-like iso-surface; C = 1/|v| inside the cell's energy span, else 0.
function (C::_GGRContribution{1})(δE::Float64)
    a = C.v[1]
    return abs(δE) ≤ C.w[1] ? inv(a) : 0.0 # a == 0 ⇒ only δE == 0 reaches here ⇒ Inf
end

# 2D: line-segment length of {v·Δk = δE} ∩ square, divided by |v|. With a ≥ b ≥ 0 and
# breakpoints w₁ = (a−b)h, w₂ = (a+b)h:
#   δ ≤ w₁      : C = 2h/a            (segment spans both edges; plateau)
#   w₁ < δ ≤ w₂ : C = (w₂ − δ)/(ab)   (segment cuts a corner)
function (C::_GGRContribution{2})(δE::Float64)
    a, b = C.v
    w1, w2 = C.w
    h = C.h
    δ = abs(δE)
    if δ ≤ w1
        return 2h / a            # a == 0 ⇒ only δ == 0 reaches here ⇒ Inf
    elseif b > 0 && δ ≤ w2
        return (w2 - δ) / (a * b)
    else
        return 0.0
    end
end

# 3D: cross-sectional area of {v·Δk = δE} ∩ cube, divided by |v|. With v1 ≥ v2 ≥ v3 ≥ 0,
# breakpoints w₁ ≤ w₂ ≤ w₃ ≤ w₄, and `V2 = |v|²` (entering via the (h|v|)² terms).
function (C::_GGRContribution{3})(δE::Float64)
    v1, v2, v3 = C.v
    w1, w2, w3, w4 = C.w
    V2, h = C.V2, C.h
    δ = abs(δE)
    if δ ≤ w1
        if v1 ≥ v2 + v3 # iso-surface fully spans the cube (planar slab); no v2,v3 division
            return 4h^2 / v1 # v1 == 0 ⇒ only δ == 0 reaches here ⇒ Inf
        else
            return (2h^2 * (v1*v2 + v2*v3 + v3*v1) - (δ^2 + h^2*V2)) / (v1*v2*v3)
        end
    elseif δ ≤ w2
        return ((h^2*(v1*v2 + 3v2*v3 + v3*v1) - h*δ*(-v1 + v2 + v3) - (δ^2 + h^2*V2)/2) /
                (v1*v2*v3))
    elseif δ ≤ w3
        return 2 * (h^2*(v1 + v2) - h*δ) / (v1*v2)
    elseif δ ≤ w4
        return (h*(v1 + v2 + v3) - δ)^2 / (2*v1*v2*v3)
    else
        return 0.0
    end
end

# ---------------------------------------------------------------------------------------- #
# Exact integral ∫_a^b C(δE) dδE of a prebuilt contribution `C`. `C` is piecewise-quadratic
# in δE with breakpoints at 0 and ±C.w, so a 3-point Simpson rule on each panel between
# consecutive breakpoints (clipped to [a,b]) is exact.
function _C_integral(a::Float64, b::Float64, C::_GGRContribution)
    s = 0.0
    p = a
    Cp = C(p)
    while p < b
        q = _next_node(p, b, C.w) # next breakpoint in (p, b], else b
        m = 0.5 * (p + q)
        Cq = C(q)
        s += (q - p) / 6 * (Cp + 4 * C(m) + Cq)
        p, Cp = q, Cq
    end
    return s
end

# smallest panel boundary strictly greater than `p` and ≤`b`: a breakpoint (±w or 0) if one
# lies in (p, b), otherwise `b` itself
@inline function _next_node(p::Float64, b::Float64, bps::NTuple{N, Float64}) where {N}
    q = b
    (p < 0.0 < q) && (q = 0.0)
    for w in bps
        (p < w < q) && (q = w)
        (p < -w < q) && (q = -w)
    end
    return q
end
