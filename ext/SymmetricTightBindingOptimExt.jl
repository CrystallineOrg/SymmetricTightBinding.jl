module SymmetricTightBindingOptimExt

using SymmetricTightBinding
using SymmetricTightBinding: solve
using LinearAlgebra: eigen!, eigvals!, Hermitian, diag, dot, pinv, tr
using Optim
using Optim: NLSolversBase
import SymmetricTightBinding: fit, multistart_fit, make_fit_objective, spectralmoments
import SymmetricTightBinding: isolate_irrep
using SymmetricTightBinding: IrrepIsolationObjective, uniform_kmesh

# Basin-hopping exploration constants (cf. `multistart_fit` & `hop_start`).
# A "hop" perturbs the incumbent best by `step .* randn .* (abs.(best_cs) .+ FLOOR*ρ)`,
# with `step = BASE*(1+n)^(1/4)` after `n` stagnated trials and ρ the moment-derived
# coefficient scale (cf. `SpectralMoments`). The two hop constants play different,
# dimensionally distinct roles:
# - `BASIN_HOP_BASE` is a multiplicative gain: hop sizes relative to each coefficient's own
#   magnitude (the initial hop perturbs each coefficient by ~30% of itself), growing slowly
#   under stagnation and reset upon meaningful improvement;
# - `BASIN_HOP_FLOOR` is an additive offset (in units of ρ) inside the per-coefficient
#   amplitude: it keeps coefficients with `|csᵢ| ≈ 0` explorable — a purely relative
#   perturbation could never move an exactly-zero coefficient.
const BASIN_RESTART_EVERY = 4     # default of the `restart_every` keyword argument
const BASIN_HOP_BASE      = 0.3   # relative hop amplitude (multiplicative gain)
const BASIN_HOP_FLOOR     = 0.25  # hop-amplitude floor, in units of the moment scale ρ
const BASIN_IMPROVE_RTOL  = 0.005 # relative improvement resetting the stagnation counter

# ---------------------------------------------------------------------------------------- #
# Define loss as sum of absolute squared error (MSE, up to scaling)

function fgh!(
    F, G, H, cs, cache::TightBindingCache{D, S}, Em_r;
    lasso::Union{Nothing,Real} = nothing
) where {D, S}
    S ≠ HERMITIAN && error("loss function can currently only handle HERMITIAN models")
    isnothing(G) || fill!(G, zero(eltype(G)))
    isnothing(H) || fill!(H, zero(eltype(H)))

    for (κ, Es_r) in enumerate(eachrow(Em_r))
        Hₖ = cache(cs, κ) # assembled into `cache`'s work array (mutated by eigen below)
        if isnothing(G) && isnothing(H)
            # fast-path, avoiding eigenvector eval. if this is a residuals-only calculation
            Es = eigvals!(Hₖ)
            us = nothing
        else
            Es, us = eigen!(Hₖ) # no Bloch phases, deliberately
        end

        # MSE loss
        if !isnothing(F)
            F += sum(abs2∘splat(-), zip(Es_r, Es))
        end

        # loss gradient & Gauss–Newton approx. of Hessian (⋆)
        if !isnothing(G) || !isnothing(H)
            ∇Es = energy_gradient_wrt_hopping(cache, κ, (Es, us))
            for (E_r, E, ∇E) in zip(Es_r, Es, ∇Es)
                isnothing(G) || (G .+= (-2 * (E_r - E)) .* ∇E)
                isnothing(H) || (H .+= 2 .* ∇E .* ∇E')
            end
        end
    end

    # lasso penalty term & gradient
    if !isnothing(lasso)
        lasso *= size(Em_r, 1) # rescaling `lasso` weight to ensure relative contributions
                               # of LSE vs LASSO are invariant to number of k-points
        !isnothing(F) && (F += lasso * sum(abs, cs))
        !isnothing(G) && (G .+= lasso .* sign.(cs))
    end

    return F
end

"""
    make_fit_objective(_fgh!)
    make_fit_objective(cache::TightBindingCache, Em_r; lasso=nothing)
    make_fit_objective(tbm::TightBindingModel, Em_r, ks; lasso=nothing)

Bundle a loss into a single Optim.jl objective, suitable for [`multistart_fit`](@ref).

The generic form takes any loss following the `fgh!` calling convention
`_fgh!(F, G, H, cs)`: `G`/`H` are filled in-place when non-`nothing` (`_fgh!` can, and for
performance usually should, be an in-place mutating function) and the loss value is returned
when `F` is non-`nothing`. Optim extracts whichever value/gradient/Hessian combination each
optimizer actually needs, passing `nothing` for the rest, so the same objective serves
zeroth-, first-, and second-order optimizers alike; for the latter, a Gauss–Newton `H` makes
e.g. `Newton()` act as Gauss–Newton & `NewtonTrustRegion()` as Levenberg–Marquardt.

The two remaining forms build the sorted-eigenvalue least-squares loss that `fit` itself
uses, comparing against the reference spectrum `Em_r` over the k-points of `cache` (or over
`ks`, tabulating a [`TightBindingCache`](@ref) internally). The `lasso` keyword adds an
``\\ell_1`` penalty on the hopping amplitudes, as in `fit`.
"""
make_fit_objective(_fgh!) = NLSolversBase.only_fgh!(_fgh!)
function make_fit_objective(cache::TightBindingCache, Em_r;
                            lasso::Union{Nothing,Real} = nothing)
    return make_fit_objective((F, G, H, cs) -> fgh!(F, G, H, cs, cache, Em_r; lasso))
end
function make_fit_objective(tbm::TightBindingModel, Em_r, ks;
                            lasso::Union{Nothing,Real} = nothing)
    return make_fit_objective(TightBindingCache(tbm, ks), Em_r; lasso)
end

# (⋆): The Gauss–Newton approximation of the Hessian `H` of a least-squares error (LSE) is 
# `∇²F ≈ 2 ∑ₖ ∑ₙ ∇Eₙ∇Eₙᵀ` (F denoting LSE loss value), which corresponds to dropping the
# residual-curvature term `2∑ₖ∑ₙ (Eₙ-Eₙʳ)∇²Eₙ`. If we use this approximation with Optim's
# second-order optimizers (`Newton`, `NewtonTrustRegion`), we effectively end up with
# Gauss–Newton (line-search) or Levenberg–Marquardt-like (trust-region) least-squares
# solvers. 
# There is a detail around the LASSO penalty though: in principle, this is not an LSE-like
# term, but it also doesn't contribute any curvature (|·| has zero second derivative), so 
# `H` only carries the Gauss–Newton term.
# The mapping from Newton method to Gauss-Newton is e.g. described e.g., at
# https://en.wikipedia.org/wiki/Gauss–Newton_algorithm#Derivation_from_Newton%27s_method

# ---------------------------------------------------------------------------------------- #
# Spectral-moment machinery, shared by `fit`'s initialization & exploration draws

"""
    SpectralMoments(tbm::TightBindingModel, Em_r, ks)

Moment-derived quantities of a band-fitting problem, exploiting linearity of
H(k) = ∑ᵢ cᵢhᵢ(k): the first moment ∑ₙ Eₙ(k) = tr H(k) = ∑ᵢ cᵢ tr hᵢ(k) is *linear* in the
coefficients `cs`, and the aggregate second moment ∑ₖₙ Eₙ²(k) = ∑ₖ tr H²(k) = cᵀQ̄c is
quadratic, with Gram matrix Q̄ᵢⱼ = ∑ₖ tr[hᵢ(k)hⱼ(k)]. Both moments are sorting-free (basis-
independent) functions of the coefficients, and hence cheap, exact handles on the fitting
problem's geometry.

## Fields
- `c₀`: least-squares solution of the (linear) first-moment (trace) equations
- `Q̄`: the Gram matrix above
- `M₂`: the reference spectrum's aggregate second moment ∑ₖₙ [Eₙʳ(k)]²
- `ρ`: scalar coefficient scale √(M₂/tr Q̄)
- `ρs`: per-term coefficient scales √(M₂/(Nᶜ·Q̄ᵢᵢ)), i.e., the size of `csᵢ` if term `i`
  were to carry an equal (1/Nᶜ) share of `M₂`; useful as per-term bounds for bounded
  global-search methods
"""
struct SpectralMoments
    c₀ :: Vector{Float64}
    Q̄  :: Matrix{Float64}
    M₂ :: Float64
    ρ  :: Float64
    ρs :: Vector{Float64}
end

function SpectralMoments(cache::TightBindingCache, Em_r::AbstractMatrix{<:Real})
    hs = cache.hs # hᵢ(k), as hs[κ][i]; Hermitian & coefficient-independent (cached)
    Nᶜ = length(cache.tbm)
    A  = [real(tr(hs[κ][i])) for κ in eachindex(hs), i in 1:Nᶜ] # Nᵏ×Nᶜ: A[κ,i] = tr hᵢ(kᴋ)
    Q̄  = [sum(hsₖ -> real(dot(hsₖ[i], hsₖ[j])), hs) for i in 1:Nᶜ, j in 1:Nᶜ]
    m₁ = vec(sum(Em_r; dims = 2))  # ∑ₙ Eₙʳ(k), for each k
    M₂ = sum(abs2, Em_r)           # ∑ₖₙ [Eₙʳ(k)]²
    # min-norm least-squares seed via `pinv` (SVD): robust to rank-deficient `A`. A purely
    # off-diagonal term has tr hᵢ(k) ≡ 0 → a zero column in `A`; the first moment carries no
    # information about it, so `pinv` correctly seeds its coefficient at 0 (its exploration
    # scale still comes from `ρs`). Unlike `A \ m₁`, this never throws on a zero/collinear
    # column, whatever the shape of `A` (a plain solve throws `SingularException` once `A`
    # is square, i.e. once #k-points ≤ #terms).
    c₀ = pinv(A) * m₁
    ρ  = sqrt(M₂ / max(tr(Q̄), eps()))
    ρs = sqrt.(M₂ ./ (Nᶜ .* max.(diag(Q̄), eps())))
    return SpectralMoments(c₀, Q̄, M₂, ρ, ρs)
end
function SpectralMoments(tbm::TightBindingModel, Em_r::AbstractMatrix{<:Real}, ks)
    return SpectralMoments(TightBindingCache(tbm, ks), Em_r)
end

"""
    SpectralMoments(cache::TightBindingCache; scale::Real = 1.0)

Moment-derived quantities for a *design* problem, i.e., one posed by an objective rather
than by a reference spectrum (cf. `isolate_irrep`), and hence with no `Em_r` to take moments
of.

The Gram matrix `Q̄` is as in the reference-based constructor, but the reference moments are
replaced by a caller-supplied energy `scale`: the second moment is set to that of a spectrum
with root-mean-square band energy `scale`, and the first-moment seed `c₀` — which the
reference alone could fix — is set to zero. The exploration scales `ρ`/`ρs` then just say
what size of hopping amplitudes produces a spectrum of the requested scale, which is all
that `multistart_fit`'s random restarts need. Since `c₀ = 0` is a degenerate starting point
(an identically vanishing Hamiltonian), `multistart_fit` replaces it by a random
moment-scaled draw when handed reference-free moments.
"""
function SpectralMoments(cache::TightBindingCache; scale::Real = 1.0)
    hs = cache.hs
    Nᶜ = length(cache.tbm)
    Q̄  = [sum(hsₖ -> real(dot(hsₖ[i], hsₖ[j])), hs) for i in 1:Nᶜ, j in 1:Nᶜ]
    M₂ = length(hs) * cache.tbm.N * float(scale)^2 # ∑ₖₙ E² at an RMS band energy of `scale`
    c₀ = zeros(Float64, Nᶜ)
    ρ  = sqrt(M₂ / max(tr(Q̄), eps()))
    ρs = sqrt.(M₂ ./ (Nᶜ .* max.(diag(Q̄), eps())))
    return SpectralMoments(c₀, Q̄, M₂, ρ, ρs)
end
function spectralmoments(cache::TightBindingCache; kws...)
    return SpectralMoments(cache; kws...)
end

# implementation of the `SymmetricTightBinding.spectralmoments` stub: the struct itself is
# deliberately extension-local (small inter-package interface); the stub lets dependent
# packages (e.g., PhotonicTightBinding.jl) construct it for use with `multistart_fit`
function spectralmoments(cache::TightBindingCache, Em_r::AbstractMatrix{<:Real})
    return SpectralMoments(cache, Em_r)
end
function spectralmoments(tbm::TightBindingModel, Em_r::AbstractMatrix{<:Real}, ks)
    return SpectralMoments(tbm, Em_r, ks)
end

# fresh moment-scaled start: displace `c₀` along a random direction `v` with magnitude `α`
# solved from (c₀+αv)ᵀQ̄(c₀+αv) = M₂, i.e., every start reproduces the reference's second
# moment and thereby has the right overall spectral scale. If no real root exists (the
# trace fit alone overshoots M₂), the displacement is scaled to M₂ instead.
# NB: a negative root is fine; it merely flips the (sign-symmetric) direction `v`
function moment_start(m::SpectralMoments)
    v = randn(length(m.c₀))
    q₂ = max(dot(v, m.Q̄, v), eps())
    q₁ = 2 * dot(m.c₀, m.Q̄, v)
    disc = q₁^2 - 4 * q₂ * (dot(m.c₀, m.Q̄, m.c₀) - m.M₂)
    α = disc ≥ 0 ? (-q₁ + sqrt(disc)) / (2q₂) : sqrt(m.M₂ / q₂)
    return m.c₀ .+ α .* v
end

# basin hop: perturb the incumbent best with relative amplitudes, floored by the moment
# scale ρ (so vanishing coefficients still get explored) & growing under stagnation
function hop_start(m::SpectralMoments, best_cs::AbstractVector, since_improve::Integer)
    step = BASIN_HOP_BASE * (1 + since_improve)^(1/4)
    return best_cs .+ step .* randn(length(best_cs)) .* (abs.(best_cs) .+ BASIN_HOP_FLOOR * m.ρ)
end

"""
    fit(tbm::TightBindingModel{D},
        Em_r::AbstractMatrix{<:Real},
        ks::AbstractVector{<:ReciprocalPointLike{D}},
        kws...)                                  --> ParameterizedTightBindingModel{D}

Fit the hopping amplitudes of a tight-binding model `tbm` to the reference energies `Em_r`,
assumed sampled over **k**-points `ks`. `Em_r[i,n]` denotes the band energy at `ks[i]` in
band `n` (and bands are assumed energetically sorted).

Fitting is performed using a local optimizer (configurable via `optimizer` from Optim.jl)
with mean-squared error loss. The local optimizer is used as a basis for a "multi-start"
global optimization, combining moment-seeded starts with basin-hopping exploration: the
first trial starts from the least-squares solution of the (linear) first-moment (trace)
equations; most subsequent trials perturb the incumbent best fit ("hops"), with amplitudes
that grow under stagnation, interspersed with fresh random restarts whose magnitudes are
chosen to reproduce the reference spectrum's second moment.
The global search returns early if the mean fit error, per band and per **k**-point, is less
than `atol`.

The function is defined as an Optim.jl extension to SymmetricTightBinding.jl: i.e., Optim.jl
must be explicitly loaded to use this function.

## Keyword arguments
- `optimizer` (default, `Optim.NewtonTrustRegion()`): a local optimizer from Optim.jl.
  First-order optimizers exploit the analytic (Feynman–Hellmann) gradient of the loss;
  second-order optimizers (e.g., `Newton()`, `NewtonTrustRegion()`) additionally exploit a
  Gauss–Newton approximation of its Hessian, `∇²F ≈ 2∑ₖ∑ₙ ∇Eₙ∇Eₙᵀ`, and thereby act as
  Gauss–Newton (line-search) or Levenberg–Marquardt-like (trust-region) least-squares
  solvers.
- `max_multistarts` (default, `max(100, 25length(tbm))`): maximum number of multi-start
  trials; scales with the number of hopping terms, since higher-dimensional models require
  more exploration.
- `restart_every` (default, `4`): every `restart_every`th multi-start trial is a fresh
  moment-scaled restart rather than a basin-hopping perturbation of the incumbent best
  (see above). Set to e.g. `typemax(Int)` to disable fresh restarts entirely (pure basin
  hopping).
- `atol` (default, `1e-3`): threshold for early return, specifying the mean energetic error
  (averaged over bands and **k**-points) below which the search stops.
- `verbose` (default, `false`): whether to print information on optimization progress.
- `options` (default, `Optim.Options(g_abstol=1e-2, f_reltol=1e-5)`): an `Optim.Options(…)`
  structure of optimization options, used during the local optimization of the multi-start
  search. The default is deliberately loose, suiting the low-precision demands of the
  multi-start search.
- `polish` (default, `true`): whether to polish off the multi-start optimization with a
  final local optimization step, at the tighter `polish_options` tolerances. This is useful
  to ensure that the best candidate from the multi-start search is fully converged.
- `polish_options` (default, `Optim.Options()`, i.e., Optim.jl's own defaults): options for
  that final step. Worth loosening — `Optim.Options(g_abstol=1e-6, f_reltol=1e-9)`, say —
  whenever the loss is not smooth enough for the gradient norm to reach Optim's default
  `g_abstol` of `1e-8`: the polish then simply runs to the 1000-iteration cap, which can
  cost more time than the entire multi-start search that preceded it.
- `lasso` (default, `nothing`): if set to a positive number, applies a LASSO penalty to the
  hopping amplitudes, encouraging model sparsity (i.e., encouraging small hopping amplitudes
  to vanish). Setting to `nothing` disables the LASSO penalty.
- `init` (default, `nothing`): if provided, a vector of hopping amplitudes used as the
  deterministic first multi-start trial, replacing the first-moment (trace) fit.

## Example

As a synthetic example, we might use `fit` to recover the coefficients of a randomly
parameterized tight-binding model, using its spectrum sampled over 10 **k**-points:

```jldoctest
julia> using Crystalline, SymmetricTightBinding, Brillouin, Optim, Random

julia> sgnum = 221;

julia> brs = calc_bandreps(sgnum);

julia> cbr = @composite brs[1] + brs[7];

julia> tbm = tb_hamiltonian(cbr);

julia> Random.seed!(123);

julia> ptbm_r = tbm(randn(length(tbm)))
4-term 6×6 ParameterizedTightBindingModel{3} (hermitian) over (3d|A₁g)⊕(3d|B₂g) with amplitudes:
 [-0.64573, -1.4633, -1.6236, -0.21767]

julia> kp = irrfbz_path(sgnum, directbasis(sgnum, Val(3)));

julia> ks = interpolate(kp, 10);

julia> Em_r = spectrum(ptbm_r, ks);

julia> ptbm_fit = fit(tbm, Em_r, ks)
4-term 6×6 ParameterizedTightBindingModel{3} (hermitian) over (3d|A₁g)⊕(3d|B₂g) with amplitudes:
 [-0.64573, -1.4633, -1.6236, -0.21767]

julia> ptbm_fit.cs ≈ ptbm_r.cs
true
```
"""
function fit(
    tbm::TightBindingModel{D},
    Em_r::AbstractMatrix{<:Real},
    ks::AbstractVector{<:SymmetricTightBinding.ReciprocalPointLike{D}};
    optimizer::Optim.AbstractOptimizer = NewtonTrustRegion(),
    max_multistarts::Integer = max(100, 25length(tbm)),
    restart_every::Integer = BASIN_RESTART_EVERY,
    atol::Real = 1e-3, # minimum threshold error, per k-point & per band, averaged over both
    verbose::Bool = false,
    options::Optim.Options = Optim.Options(;
        g_abstol = 1e-2,
        f_reltol = 1e-5,
    ),
    polish::Bool = true,
    polish_options::Optim.Options = Optim.Options(),
    lasso::Union{Nothing,Real} = nothing,
    init::Union{Nothing, AbstractVector{<:Real}} = nothing,
) where D

    cache = TightBindingCache(tbm, ks) # hᵢ(k) tabulated once, shared by objective & moments
    obj = make_fit_objective(cache, Em_r; lasso)
    moments = spectralmoments(cache, Em_r)
    tol = length(ks) * tbm.N * atol^2 # sum of absolute squares tolerance
    best_cs, _ = multistart_fit(obj, moments;
                                optimizer, max_multistarts, restart_every, tol, options,
                                polish, polish_options, verbose, init)
    return tbm(best_cs)
end


# ---------------------------------------------------------------------------------------- #
# multi-start driver

"""
    multistart_fit(obj, moments::SpectralMoments; kws...)  -->  (best_cs, best_loss)

Moment-seeded, basin-hopping multi-start minimization of an Optim.jl objective `obj` (e.g.,
constructed via [`make_fit_objective`](@ref)), with exploration draws derived from `moments`
(cf. `SpectralMoments`). This is the search engine underlying `fit`, factored out
so that related fitting problems with custom losses (e.g., PhotonicTightBinding.jl) can
reuse it.

The search strategy is informed by the loss-landscape geometry typical of (sorted-
eigenvalue, least-squares) band-fitting problems: the landscape is funnel-like, with
hierarchies of progressively better band-assignment local minima, so that local exploration
around the incumbent best is far more effective than independent restarts. Concretely:
- the first trial starts from the deterministic first-moment (trace) fit `moments.c₀` (or
  from `init`, if provided);
- most subsequent trials are "hops": perturbations of the incumbent best fit, with
  amplitudes that grow slowly under stagnation and reset upon meaningful improvement (cf.
  `hop_start`);
- every `restart_every`th trial is instead a fresh moment-scaled random restart (cf.
  `moment_start`), hedging against over-commitment to a single funnel.

Returns `(best_cs, best_loss)`, returning early if the loss drops to (or below) `tol`.

## Keyword arguments
- `optimizer`, `max_multistarts`, `restart_every`, `options`, `polish`, `polish_options`,
  `verbose`: as in `fit`.
- `tol` (default, `0.0`): early-return threshold on the *loss value* (for `fit`'s
  mean-error-per-band-and-k-point semantics, pass `length(ks) * N * atol^2`).
- `init` (default, `nothing`): if provided, replaces `moments.c₀` as the deterministic
  first start. If `moments.c₀` vanishes identically — as it does for the reference-free
  moments of a design problem — a random moment-scaled draw is used instead.
"""
function multistart_fit(
    obj,
    moments::SpectralMoments;
    optimizer::Optim.AbstractOptimizer = NewtonTrustRegion(),
    max_multistarts::Integer = max(100, 25length(moments.c₀)),
    restart_every::Integer = BASIN_RESTART_EVERY,
    tol::Real = 0.0,
    options::Optim.Options = Optim.Options(; g_abstol = 1e-2, f_reltol = 1e-5),
    polish::Bool = true,
    polish_options::Optim.Options = Optim.Options(),
    verbose::Bool = false,
    init::Union{Nothing, AbstractVector{<:Real}} = nothing,
)
    best_cs = Vector{Float64}(undef, length(moments.c₀))
    best_loss = Inf
    since_improve = 0 # trials since the last meaningful improvement (sets hop amplitude)
    verbose && println("Starting multi-start optimization with $max_multistarts trials:")
    for t in 1:max_multistarts
        kind = t == 1 ? :start : (t % restart_every == 0 ? :restart : :hop)
        verbose && (print("   trial #$t "); printstyled("[", kind, "]"; color=:light_blue))
        init_cs = if kind === :start
            # deterministic first start: user-provided `init` or first-moment (trace) fit —
            # except for reference-free moments, whose `c₀` is an identically vanishing (and
            # hence useless) Hamiltonian, which we replace by a random moment-scaled draw
            isnothing(init) ? (iszero(moments.c₀) ? moment_start(moments) :
                                                    copy(moments.c₀)) :
                              convert(Vector{Float64}, init)
        elseif kind === :hop
            hop_start(moments, best_cs, since_improve)
        else # kind === :restart
            moment_start(moments)
        end
        o = optimize(obj, init_cs, optimizer, options)
        accept = o.minimum < best_loss
        # only meaningful (> BASIN_IMPROVE_RTOL, relatively) improvements reset the hop
        # amplitude; marginal ones are still kept as the incumbent best below
        meaningful = o.minimum * (1 + BASIN_IMPROVE_RTOL) < best_loss
        since_improve = meaningful ? 0 : since_improve + 1

        if verbose
            printstyled(" (loss ", round(o.minimum; sigdigits = 3), ")";
                        color = :light_black)
            accept && printstyled(" → new best"; color = :green)
            println()
        end

        accept || continue # discard local optimization; not better globally

        best_loss = o.minimum
        best_cs = o.minimizer

        if best_loss ≤ tol
            if verbose
                printstyled("   tolerance met: returning\n"; color = :green, bold = true)
            end
            break
        end
    end
    if verbose && isfinite(tol) && best_loss > tol
        printstyled(
            "   `max_multistarts` exceeded: tolerance not met\n   (consider increasing the number of tight-binding terms)\n";
            color = :yellow,
        )
    end

    # polish off the best result
    if polish
        verbose && print("Polishing off ")
        o = optimize(obj, best_cs, optimizer, polish_options)
        if o.minimum < best_loss
            best_loss = o.minimum
            best_cs = o.minimizer
            if verbose
                printstyled(
                    "(loss reduced to ", round(best_loss; sigdigits = 3), ")\n";
                    color = :green,
                )
            end
        elseif verbose
            printstyled("(no improvement)\n"; color = :yellow)
        end
    end

    return best_cs, best_loss
end

# ---------------------------------------------------------------------------------------- #
# design drivers

"""
    isolate_irrep(tbm::TightBindingModel{D, HERMITIAN}, irlab::AbstractString; kws...)
                                                    --> ParameterizedTightBindingModel{D}

Search for hopping amplitudes of `tbm` that make the band multiplet transforming as the
irrep `irlab` (e.g., `"Γ₄⁺"`) *energy-isolated*: away from the irrep's own **k**-point, no
band may attain the multiplet's energy — the multiplet's own bands included, which must
split cleanly into a set below and a set above it.

This is the design counterpart of [`fit`](@ref): the least-squares loss against a reference
spectrum is replaced by [`IrrepIsolationObjective`](@ref) — a smooth soft-minimum of the
(scale-normalized) energetic separation between the multiplet and every other band — but the
search engine, `multistart_fit`, is the same. Whether the returned model actually achieves
isolation, and by what margin, is not guaranteed by construction and should be verified
independently, on a finer **k**-mesh, with [`isolation_report`](@ref).

The function is defined as an Optim.jl extension to SymmetricTightBinding.jl: i.e., Optim.jl
must be explicitly loaded to use this function.

## Keyword arguments
- `Nk`, `ks`, `kexclude`, `split`, `β`, `λ_scale`, `σ₀`, `λ_shift`, `lasso`: forwarded to
  [`IrrepIsolationObjective`](@ref); in particular, `Nk` (default, `6`) sets the density of
  the **k**-mesh that isolation is imposed over -- or `ks` the **k**-points outright, e.g.,
  a decimated `irrfbz_path` from Brillouin.jl -- and `kexclude` (default, `0.1`) the radius
  of the neighbourhood of the irrep's **k**-point that is exempted from the constraints
  (necessarily so: the gap closes at the degeneracy itself). `lasso` (default, `nothing`)
  adds a scale-invariant ``\\ell_1`` penalty on the amplitudes, trading attained gap for a
  sparser model, exactly as in [`fit`](@ref).
- `scale` (default, `1.0`): the energy scale of the spectra explored by the multi-start
  restarts (cf. `SpectralMoments`). Only meaningful relative to `σ₀`, since the objective
  itself is scale-invariant.
- `optimizer` (default, `LBFGS()`): a *first-order* Optim.jl optimizer; the objective
  supplies an analytic (Feynman–Hellmann) gradient, but no Hessian.
- `max_multistarts` (default, `max(60, 10length(tbm))`), `restart_every`, `tol` (default,
  `-Inf`, i.e., no early return), `options`, `polish`, `verbose`: as in `multistart_fit`.
  Design searches are markedly harder than fits and generally need *many* starts: the SG 224
  example below finds nothing at all with 50 starts and succeeds with 200.
- `polish_options` (default, `Optim.Options(g_abstol=1e-6, f_reltol=1e-9)`): deliberately
  looser than Optim.jl's defaults, which the soft-minimum loss cannot meet — it would run
  to the 1000-iteration cap and spend more time than the whole search, for a fraction of a
  percent of loss.
- `init` (default, `nothing`): a warm start for the deterministic first trial, given either
  as a vector of amplitudes or as a [`ParameterizedTightBindingModel`](@ref) to take them
  from. Subsequent trials still hop around it, so a warm start biases the search without
  confining it — useful for refining an earlier result against a denser `ks`, a smaller
  `kexclude`, or a `lasso` weight, each of which changes the landscape only gradually.

## Example

```julia-repl
julia> using Crystalline, SymmetricTightBinding, Optim, Random

julia> brs = calc_bandreps(224);

julia> cbr = @composite brs[18] + brs[23]; # (4b|A₂g) + (2a|A₂), 6 bands

julia> tbm = tb_hamiltonian(cbr, [[0,0,0], [1,0,0]]);

julia> Random.seed!(1); # `Nk` divisible by 6, and many starts: cf. the keywords above

julia> ptbm = isolate_irrep(tbm, "Γ₄⁺"; Nk = 6, max_multistarts = 200);

julia> report = isolation_report(ptbm, "Γ₄⁺"; Nk = 24);

julia> report.isolated, report.bands, report.split # 2 of the 3 bands below Et, 1 above
(true, 3:5, 2)
```
"""
function isolate_irrep(
    tbm::TightBindingModel{D, HERMITIAN},
    irlab::AbstractString;
    Nk::Integer = 6,
    ks::AbstractVector = uniform_kmesh(Val(D), Nk),
    kexclude::Real = 0.1,
    split::Union{Nothing, Integer} = nothing,
    β::Real = 20.0,
    λ_scale::Real = 0.1,
    σ₀::Real = 1.0,
    λ_shift::Real = 0.1,
    lasso::Union{Nothing, Real} = nothing,
    scale::Real = 1.0,
    optimizer::Optim.AbstractOptimizer = LBFGS(),
    max_multistarts::Integer = max(60, 10length(tbm)),
    restart_every::Integer = BASIN_RESTART_EVERY,
    tol::Real = -Inf,
    options::Optim.Options = Optim.Options(; g_abstol = 1e-4, f_reltol = 1e-6),
    polish::Bool = true,
    polish_options::Optim.Options = Optim.Options(; g_abstol = 1e-6, f_reltol = 1e-9),
    verbose::Bool = false,
    init::Union{Nothing, AbstractVector{<:Real}, ParameterizedTightBindingModel{D}} = nothing,
) where D
    _fgh! = IrrepIsolationObjective(tbm, irlab; ks, kexclude, split, β, λ_scale, σ₀,
                                    λ_shift, lasso)
    obj = make_fit_objective(_fgh!)
    moments = SpectralMoments(_fgh!.cache; scale)
    init_cs = init isa ParameterizedTightBindingModel ? init.cs : init
    best_cs, _ = multistart_fit(obj, moments;
                                optimizer, max_multistarts, restart_every, tol, options,
                                polish, polish_options, verbose, init = init_cs)
    return tbm(best_cs)
end

end # module SymmetricTightBindingOptimExt