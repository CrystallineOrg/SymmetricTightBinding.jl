
using SymmetricTightBinding
using Optim
using LinearAlgebra: eigen!, Hermitian

# ---------------------------------------------------------------------------------------- #
# Define loss as sum of absolute squared error (MSE, up to scaling)

function fg!(F, G, cs, tbm::TightBindingModel, Em_r, ks, μᴸ; λ = 1)
    ptbm = tbm(cs)
    if !isnothing(G)
        fill!(G, zero(eltype(G)))
    end

    for (Es_r, k) in zip(eachrow(Em_r), ks)
        H = Hermitian(ptbm(k))
        Es, us = eigen!(H) # no Bloch phases, deliberately
        Esᴸ = @view Es[1:μᴸ]     # longitudinal bands
        Esᵀ = @view Es[μᴸ+1:end] # regular, transverse bands
        # MSE loss
        if !isnothing(F)
            F += sum(abs2∘splat(-), zip(Es_r, Esᵀ); init=zero(F)) # regular loss
            F += λ * sum(E -> max(zero(E), E)^2, Esᴸ)             # longitudinal loss
        end

        # gradient of MSE loss
        if !isnothing(G)
            ∇Es = energy_gradient_wrt_hopping(ptbm, k, (Es, us))
            ∇Esᴸ = @view ∇Es[1:μᴸ]
            ∇Esᵀ = @view ∇Es[μᴸ+1:end]
            @assert size(Esᵀ) == size(∇Esᵀ) == size(Es_r)
            for (E_r, E, ∇E) in zip(Es_r, Esᵀ, ∇Esᵀ) # regular loss gradient
                G .+= (-2 * (E_r - E)) .* ∇E
            end
            for (E, ∇E) in zip(Esᴸ, ∇Esᴸ)            # longitudinal loss gradient
                if E > 0
                    G .+= (2λ * E) .* ∇E
                end
            end
        end
    end
    return F
end

"""
    fit(tbm::TightBindingModel{D},
        freqs_r::AbstractMatrix{<:Real},
        ks::AbstractVector{<:ReciprocalPointLike{D}},
        kws...)                                  --> ParameterizedTightBindingModel{D}

Fit the hopping amplitudes of a tight-binding model `tbm` to the reference frequencies `freqs_r`,
assumed sampled over **k**-points `ks`. `freqs_r[i,n]` denotes the band frequency at `ks[i]` in
band `n` (and bands are assumed energetically sorted).

Fitting is performed using a local optimizer (configurable via `optimizer` from Optim.jl)
with mean-squared error loss. The local optimizer is used as a basis for a "multi-start"
global optimization.
The global search returns early if the mean fit error, per band and per frequency, is less than
`atol`.

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


## Notes
The frequencies are provided by the user but the energies are used internally to do the fitting.
The tight-binding model energies (E) are compared to squared frequencies (ω²), so the provided
frequencies are squared before fitting.
```
"""
function photonic_fit(
    tbm::TightBindingModel{D},
    freqs_r::AbstractMatrix{<:Real},
    ks::AbstractVector{<:ReciprocalPointLike{D}};
    optimizer::Optim.FirstOrderOptimizer = LBFGS(),
    atol::Real = 1e-3, # minimum threshold error, per k-point & per band, averaged over both
    max_multistarts::Integer = 100,
    verbose::Bool = false,
    loss_penalty_weight::Real = LOSS_PENALTY_WEIGHT,
    options::Optim.Options = Optim.Options(;
        g_abstol = 1e-2,
        f_reltol = 1e-5,
        ),
    polish::Bool = true,
) where D
    # convert frequencies to energies and sort them
    Em_r = freqs_r .^ 2
    sort!(Em_r; dims = 2)

    μᴸ = tbm.N - size(Em_r, 2) # number of longitudinal bands
    # let-block-capture-trick to make absolutely sure we have no closure boxing issues
    _fg! = let tbm = tbm, Em_r = Em_r, ks = ks, μᴸ = μᴸ, λ = loss_penalty_weight
        (F, G, cs) -> fg!(F, G, cs, tbm, Em_r, ks, μᴸ; λ)
    end

    # multi-start optimization
    n_fit = size(Em_r, 2) # number of bands to fit
    tol = length(ks) * n_fit * atol^2 # sum of absolute squares tolerance
    best_cs = Vector{Float64}(undef, length(tbm))
    best_loss = Inf
    init_hopping_scale = sum(Em_r) / length(Em_r) * 0.25
    init_cs = randn(length(tbm)) .* init_hopping_scale
    since_last_improvement = 0
    verbose && println("Starting multi-start optimization with $max_multistarts trials:")
    for t in 1:max_multistarts
        verbose && print("   trial #$t")
        o = optimize(Optim.only_fg!(_fg!), init_cs, optimizer, options)
        accept = o.minimum < best_loss
        
        if verbose
            mean_err = round(sqrt(o.minimum / (n_fit * length(ks))); sigdigits = 3)
            printstyled(" (mean err ", mean_err, ")"; color = :light_black)
            accept && printstyled(" → new best"; color = :green)
            println()
        end

        if accept
            best_loss = o.minimum
            best_cs = o.minimizer
            since_last_improvement = 0
            if best_loss ≤ tol
                if verbose
                    printstyled("   tolerance met: returning\n"; color = :green, bold = true)
                end
                break
            end
        end

        # a simple basin-hopping exploration strategy
        since_last_improvement += 1
        step_scale = since_last_improvement.^(1/4) * 0.5 / length(tbm)
        init_cs = best_cs .+ step_scale .* randn(length(tbm)) .* abs.(best_cs)
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
                round(sqrt(o.minimum / (n_fit * length(ks))); sigdigits = 3), ")\n";
                color = :green
            )
        end
    end
        
    return tbm(best_cs)
end