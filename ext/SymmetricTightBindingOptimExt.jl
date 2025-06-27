module SymmetricTightBindingOptimExt

using SymmetricTightBinding
using SymmetricTightBinding: solve
using LinearAlgebra: eigen!, Hermitian
using Optim
import SymmetricTightBinding: fit

# ---------------------------------------------------------------------------------------- #
# Define loss as sum of absolute squared error (MSE, up to scaling)

function fg!(
    F, G, cs, tbm::TightBindingModel, Em_r, ks;
    lasso::Union{Nothing,Real} = nothing
)
    ptbm = tbm(cs)
    if !isnothing(G)
        fill!(G, zero(eltype(G)))
    end

    for (Es_r, k) in zip(eachrow(Em_r), ks)
        H = Hermitian(ptbm(k))
        Es, us = eigen!(H) # no Bloch phases, deliberately

        # MSE loss (possibly with lasso penalty)
        if !isnothing(F)
            F += sum(abs2∘splat(-), zip(Es_r, Es))
            if !isnothing(lasso)
                F += lasso * sum(abs, cs)
            end
        end

        # gradient of loss
        if !isnothing(G)
            ∇Es = energy_gradient_wrt_hopping(ptbm, k, (Es, us))
            for (E_r, E, ∇E) in zip(Es_r, Es, ∇Es)
                G .+= (-2 * (E_r - E)) .* ∇E
                if !isnothing(lasso)
                    G .+= lasso .* sign.(cs) # lasso penalty gradient
                end
            end
        end
    end
    return F
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
global optimization.
The global search returns early if the mean fit error, per band and per energy, is less than
`atol`.

The function is defined as an Optim.jl extension to SymmetricTightBinding.jl: i.e., Optim.jl
must be explicitly loaded to use this function.

## Keyword arguments
- `optimizer` (default, `Optim.LBFGS()`): a local optimizer from Optim.jl, capable of
  exploiting gradient information.
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
    optimizer::Optim.FirstOrderOptimizer = LBFGS(),
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

    # let-block-capture-trick to make absolutely sure we have no closure boxing issues
    _fg! = let Em_r = Em_r, ks = ks, tbm = tbm, lasso = lasso
        (F, G, cs) -> fg!(F, G, cs, tbm, Em_r, ks; lasso)
    end

    # multi-start optimization
    tol = length(ks) * tbm.N * atol^2 # sum of absolute squares tolerance
    best_cs = Vector{Float64}(undef, length(tbm))
    best_loss = Inf
    init_hopping_scale = sum(Em_r) / length(Em_r) * 0.25
    verbose && println("Starting multi-start optimization with $max_multistarts trials:")
    for t in 1:max_multistarts
        verbose && print("   trial #$t")
        init_cs = randn(length(tbm))
        init_cs .*= init_hopping_scale # TODO: Improve guess; could likely do much better
        o = optimize(Optim.only_fg!(_fg!), init_cs, optimizer, options)
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
        o = optimize(Optim.only_fg!(_fg!), best_cs, optimizer)
        o.minimum > best_loss && (best_loss = o.minimum; best_cs = o.minimizer)
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