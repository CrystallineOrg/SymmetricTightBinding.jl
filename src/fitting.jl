
using SymmetricTightBinding
using SymmetricTightBinding: ReciprocalPointLike
using Optim
using PythonCall: pyconvert

# ---------------------------------------------------------------------------------------- #
# Define loss as sum of absolute squared error (MSE, up to scaling)

# MSE loss
function longitudinal_loss(ks, tbm, cs, μᴸ; λ = 1)
    L = zero(Float64)
    for k in ks
        Es = spectrum(tbm(cs), k)
        Es_extra = @view Es[1:μᴸ] # extra bands

        L += sum(E -> max(zero(E), E)^2, Es_extra) # penalty for extra bands above 0
    end
    return λ*L
end

function grad_longitudinal_loss!(ks, tbm, cs, μᴸ, G = zeros(Float64, length(cs)); λ = 1)
    ptbm = tbm(cs)
    for k in ks
        Es = spectrum(ptbm, k)
        ∇Es = energy_gradient_wrt_hopping(ptbm, k)
        # gradient of the penalty for positive extra bands
        for i in 1:μᴸ
            if Es[i] > 0
                G .+= (2λ * Es[i]) .* ∇Es[i]
            end
        end
    end
    return G
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
- `options` (default, empty): a `Optim.Options(…)` structure of optimization options.
- `max_multistarts` (default, `150`): maximum number of multi-start iterations.
- `atol` (default, `1e-3`): threshold for early return, specifying the minimum required mean
  energetic error (averaged over bands and **k**-points).
- `verbose` (default, `false`): whether to print information on optimization progress.


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
    options::Optim.Options = Optim.Options(),
    max_multistarts::Integer = 150,
    atol::Real = 1e-3, # minimum threshold error, per k-point & per band, averaged over both
    verbose::Bool = false,
    loss_penalty_weight::Real = LOSS_PENALTY_WEIGHT,
) where D
    # convert frequencies to energies and sort them
    Em_r = freqs_r .^ 2
    Em_r = sort(Em_r; dims=2)

    # let-block-capture-trick to make absolutely sure we have no closure boxing issues
    μᴸ = tbm.N - size(Em_r, 2) # number of extra bands
    loss_closure = let Em_r = Em_r, ks = ks, tbm = tbm, μᴸ = μᴸ, λ = loss_penalty_weight
        cs -> loss(Em_r, ks, tbm, cs; start = μᴸ+1) + longitudinal_loss(ks, tbm, cs, μᴸ; λ)
    end
    grad_loss_closure! = let Em_r = Em_r, ks = ks, tbm = tbm, λ = loss_penalty_weight
        (G, cs) -> begin
            fill!(G, zero(eltype(G)))
            grad_loss!(Em_r, ks, tbm, cs, G; start = μᴸ+1)
            grad_longitudinal_loss!(ks, tbm, cs, μᴸ, G; λ)
        end
    end

    # multi-start optimization
    n_fit = size(Em_r, 2) # number of bands to fit
    tol = length(ks) * n_fit * atol^2 # sum of absolute squares tolerance
    best_cs = Vector{Float64}(undef, length(tbm))
    best_loss = Inf
    for t in 1:max_multistarts
        init_cs = randn(length(tbm))
        o = optimize(loss_closure, grad_loss_closure!, init_cs, optimizer, options)
        o.minimum > best_loss && continue # discard local optimization; not better globally

        if verbose
            println(
                "   Loss updated (trial $t): mean error = ",
                round(sqrt(o.minimum / (n_fit * length(ks))); sigdigits = 3),
            )
        end

        best_loss = o.minimum
        best_cs = o.minimizer

        if best_loss ≤ tol
            verbose && printstyled("      tolerance met: returning\n"; color = :green)
            break
        end
    end
    if verbose && best_loss > tol
        printstyled(
            "      `max_multistarts` exceeded: tolerance not met\n Consider increasing the number of tight-binding terms\n";
            color = :yellow,
        )
    end

    return tbm(best_cs)
end