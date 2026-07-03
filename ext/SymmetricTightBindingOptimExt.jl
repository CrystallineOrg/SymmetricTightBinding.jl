module SymmetricTightBindingOptimExt

using SymmetricTightBinding
using SymmetricTightBinding: solve
using LinearAlgebra: eigen!, Hermitian, dot, norm, tr
using Optim
import SymmetricTightBinding: fit

# ---------------------------------------------------------------------------------------- #
# Define loss as sum of absolute squared error (MSE, up to scaling)

function fgh!(
    F, G, H, cs, tbm::TightBindingModel{D, S}, Em_r, ks;
    lasso::Union{Nothing,Real} = nothing
) where {D, S}
    S ≠ HERMITIAN && error("loss function can currently only handle HERMITIAN models")
    ptbm = tbm(cs)
    isnothing(G) || fill!(G, zero(eltype(G)))
    isnothing(H) || fill!(H, zero(eltype(H)))

    for (Es_r, k) in zip(eachrow(Em_r), ks)
        Hₖ = ptbm(k)
        Es, us = eigen!(Hₖ) # no Bloch phases, deliberately

        # MSE loss
        if !isnothing(F)
            F += sum(abs2∘splat(-), zip(Es_r, Es))
        end

        # loss gradient & Gauss–Newton approx. of Hessian (⋆)
        if !isnothing(G) || !isnothing(H)
            ∇Es = energy_gradient_wrt_hopping(ptbm, k, (Es, us))
            for (E_r, E, ∇E) in zip(Es_r, Es, ∇Es)
                isnothing(G) || (G .+= (-2 * (E_r - E)) .* ∇E)
                isnothing(H) || (H .+= 2 .* ∇E .* ∇E')
            end
        end
    end

    # lasso penalty term & gradient
    if !isnothing(lasso)
        !isnothing(F) && (F += lasso * sum(abs, cs))
        !isnothing(G) && (G .+= lasso .* sign.(cs))
    end

    return F
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
global optimization, seeded by spectral-moment matching: the first start is the
least-squares solution of the (linear) first-moment (trace) equations; subsequent random
starts are displaced from it with magnitudes chosen to reproduce the reference spectrum's
second moment.
The global search returns early if the mean fit error, per band and per energy, is less than
`atol`.

The function is defined as an Optim.jl extension to SymmetricTightBinding.jl: i.e., Optim.jl
must be explicitly loaded to use this function.

## Keyword arguments
- `optimizer` (default, `Optim.NewtonTrustRegion()`): a local optimizer from Optim.jl.
  First-order optimizers exploit the analytic (Feynman–Hellmann) gradient of the loss;
  second-order optimizers (e.g., `Newton()`, `NewtonTrustRegion()`) additionally exploit a
  Gauss–Newton approximation of its Hessian, `∇²F ≈ 2∑ₖ∑ₙ ∇Eₙ∇Eₙᵀ`, and thereby act as
  Gauss–Newton (line-search) or Levenberg–Marquardt-like (trust-region) least-squares
  solvers.
- `max_multistarts` (default, `100`): maximum number of multi-start iterations.
- `atol` (default, `1e-3`): threshold for early return, specifying the minimum required mean
  energetic error (averaged over bands and **k**-points).
- `verbose` (default, `false`): whether to print information on optimization progress.
- `options` (default, empty): a `Optim.Options(…)` structure of optimization options, used
  during the local optimization of the multi-start search. Defaults to
  `Optim.Options(g_abstol=1e-2, f_reltol=1e-5)` (i.e., low tolerances, suitable for the
  low precision demands of the multi-start search).
- `polish` (default, `true`): whether to polish off the multi-start optimization with a
  final local optimization step using default Optim.jl options. This is useful to ensure
  that the best candidate from the multi-start search is fully converged.
- `lasso` (defalt, `nothing`): if set to a positive number, applies a LASSO penalty to the
  hopping amplitudes, encouraging model sparsity (i.e., small hopping amplitudes to
  vanish). Setting to `nothing` disables the LASSO penalty.

## Example

As a synthetic example, we might use `fit` to recover the coefficients of a randomly
parameterized tight-binding model, using its spectrum sampled over 10 **k**-points:

```jldoctest
julia> using Crystalline, SymmetricTightBinding, Brillouin, Optim
julia> sgnum = 221;
julia> brs = calc_bandreps(sgnum);
julia> cbr = @composite brs[1] + brs[7];
julia> tbm = tb_hamiltonian(cbr);

julia> using Random; Random.seed!(123);
julia> ptbm_r = tbm(randn(length(tbm)))
4-term 6×6 ParameterizedTightBindingModel{3} over (3d|A₁g)⊕(3d|B₂g) with amplitudes:
 [-0.64573, -1.4633, -1.6236, -0.21767]

julia> kp = irrfbz_path(sgnum, directbasis(sgnum, Val(3)));
julia> ks = interpolate(kp, 10);
julia> Em_r = spectrum(ptbm_r, ks);
julia> ptbm_fit = fit(tbm, Em_r, ks)
4-term 6×6 ParameterizedTightBindingModel{3} over (3d|A₁g)⊕(3d|B₂g) with amplitudes:
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
    max_multistarts::Integer = 100,
    atol::Real = 1e-3, # minimum threshold error, per k-point & per band, averaged over both
    verbose::Bool = false,
    options::Optim.Options = Optim.Options(;
        g_abstol = 1e-2,
        f_reltol = 1e-5,
    ),
    polish::Bool = true,
    lasso::Union{Nothing,Real} = nothing,
) where D

    # let-block-capture-trick to make absolutely sure we have no closure boxing issues;
    # second-order optimizers additionally get a Gauss–Newton Hessian (cf. `fgh!`), making
    # e.g. `Newton()` act as Gauss–Newton & `NewtonTrustRegion()` as Levenberg–Marquardt
    obj = let Em_r = Em_r, ks = ks, tbm = tbm, lasso = lasso
        if optimizer isa Optim.SecondOrderOptimizer
            Optim.only_fgh!((F, G, H, cs) -> fgh!(F, G, H, cs, tbm, Em_r, ks; lasso))
        else
            Optim.only_fg!((F, G, cs) -> fgh!(F, G, nothing, cs, tbm, Em_r, ks; lasso))
        end
    end

    # --- multi-start initialization from spectral moments ---
    # the reference spectrum's moments are simple, sorting-free functions of `cs`, by
    # linearity of H(k) = ∑ᵢ cᵢhᵢ(k): the first moment ∑ₙ Eₙ(k) = tr H(k) = ∑ᵢ cᵢ tr hᵢ(k)
    # is *linear* in `cs`, and the aggregate second moment ∑ₖₙ Eₙ²(k) = ∑ₖ tr H²(k) = cᵀQ̄c
    # is quadratic, with Gram matrix Q̄ᵢⱼ = ∑ₖ tr[hᵢ(k)hⱼ(k)]. We exploit this to seed the
    # multi-start search: the deterministic first start `c₀` is the least-squares solution
    # of the (linear) first-moment equations; subsequent random starts displace `c₀` along
    # random directions `v`, with the displacement magnitude `α` solved (quadratically)
    # such that every start reproduces the reference's second moment, i.e., has the right
    # overall spectral scale
    Nᶜ = length(tbm)
    hs = [tbm[i](k) for i in 1:Nᶜ, k in ks] # hᵢ(k); Hermitian & coefficient-independent
    A  = [real(tr(hs[i, κ])) for κ in eachindex(ks), i in 1:Nᶜ]
    Q̄  = [sum(κ -> real(dot(hs[i, κ], hs[j, κ])), eachindex(ks)) for i in 1:Nᶜ, j in 1:Nᶜ]
    m₁ = vec(sum(Em_r; dims = 2))  # ∑ₙ Eₙʳ(k), for each k
    M₂ = sum(abs2, Em_r)           # ∑ₖₙ [Eₙʳ(k)]²
    c₀ = norm(A) > 0 ? A \ m₁ : zeros(Nᶜ)

    # multi-start optimization
    tol = length(ks) * tbm.N * atol^2 # sum of absolute squares tolerance
    best_cs = Vector{Float64}(undef, length(tbm))
    best_loss = Inf
    verbose && println("Starting multi-start optimization with $max_multistarts trials:")
    for t in 1:max_multistarts
        verbose && print("   trial #$t")
        init_cs = if t == 1
            copy(c₀) # deterministic first start: exact first-moment (trace) fit
        else
            v = randn(Nᶜ)
            # solve (c₀+αv)ᵀQ̄(c₀+αv) = M₂ for α; if no real root exists (trace fit alone
            # already overshoots the second moment), scale the displacement to M₂ instead.
            # NB: a negative root is fine; it merely flips the (sign-symmetric) direction v
            q₂ = max(dot(v, Q̄, v), eps())
            q₁ = 2 * dot(c₀, Q̄, v)
            disc = q₁^2 - 4 * q₂ * (dot(c₀, Q̄, c₀) - M₂)
            α = disc ≥ 0 ? (-q₁ + sqrt(disc)) / (2q₂) : sqrt(M₂ / q₂)
            c₀ .+ α .* v
        end
        o = optimize(obj, init_cs, optimizer, options)
        accept = o.minimum < best_loss
        
        if verbose
            mse_loss = o.minimum
            if !isnothing(lasso)
                mse_loss -= lasso * sum(abs, o.minimizer)
            end
            mean_err = round(mse_loss / (tbm.N * length(ks)); sigdigits = 3)
            printstyled(" (mean err ", mean_err, ")"; color = :light_black)
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
    if verbose && best_loss > tol
        printstyled(
            "   `max_multistarts` exceeded: tolerance not met\n   (consider increasing the number of tight-binding terms)\n";
            color = :yellow,
        )
    end

    # polish off the best result
    if polish
        verbose && print("Polishing off ")
        o = optimize(obj, best_cs, optimizer)
        o.minimum < best_loss && (best_loss = o.minimum; best_cs = o.minimizer)
        if verbose
            printstyled(
                "(mean error ",
                round(o.minimum / (tbm.N * length(ks)); sigdigits = 3), ")\n";
                color = :green
            )
        end
    end

    return tbm(best_cs)
end

end # module SymmetricTightBindingOptimExt