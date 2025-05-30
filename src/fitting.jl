
using SymmetricTightBinding
using SymmetricTightBinding: ReciprocalPointLike
using Optim
using PythonCall: pyconvert

# ---------------------------------------------------------------------------------------- #
# Define loss as sum of absolute squared error (MSE, up to scaling)

# MSE loss
function photonic_loss(Em_r, ks, tbm, cs)
    L = zero(Float64)
    n_extra = length(tbm[1].axis) - size(Em_r, 2) # number of extra bands
    if n_extra > 0
        for (Es_r, k) in zip(eachrow(Em_r), ks)
            Es = spectrum(tbm(cs), k)
            Es_fit = Es[1:n_extra] # bands to fit
            Es_extra = Es[n_extra+1:end] # extra bands, not to fit

            # fit loss
            L += sum(abs2, (E_r - E for (E_r, E) in zip(Es_r, Es_fit)))

            # penalty for extra bands above 0
            penalty = sum(abs2, max.(Es_extra, 0.0))
            λ = 0.1 # penalty weight
            L += λ * penalty
        end
        return L
    end

    for (Es_r, k) in zip(eachrow(Em_r), ks)
        Es = spectrum(tbm(cs), k)
        L += sum(abs2, (E_r - E for (E_r, E) in zip(Es_r, Es)))
    end
    return L
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
The energies and frequencies are related by the equation `E = ω²`, where `E` is the energy and
`ω` is the frequency. In order to compare the spectrum of the tight-binding model with
the frequencies, the square root of the energies should be taken.
```
"""
function fit(
    tbm::TightBindingModel{D},
    freqs_r::PythonCall.Core.Py,
    ks::PythonCall.Core.Py;
    optimizer::Optim.FirstOrderOptimizer = LBFGS(),
    options::Optim.Options = Optim.Options(),
    max_multistarts::Integer = 150,
    atol::Real = 1e-3, # minimum threshold error, per k-point & per band, averaged over both
    verbose::Bool = false,
) where D

    # fix types coming from PythonCall and transform frequencies to energies
    Em_r = pyconvert(Matrix, freqs_r) .^ 2 # square the frequencies to get energies
    ks = pyconvert(Vector{ReciprocalPointLike{D}}, ks) # convert to a vector of k-points

    # let-block-capture-trick to make absolutely sure we have no closure boxing issues
    loss_closure = let Em_r = Em_r, ks = ks, tbm = tbm
        cs -> photonic_loss(Em_r, ks, tbm, cs)
    end
    grad_loss_closure! = let Em_r = Em_r, ks = ks, tbm = tbm
        (G, cs) -> grad_loss!(Em_r, ks, tbm, cs, G)
    end

    # multi-start optimization
    tol = length(ks) * tbm.N * atol^2 # sum of absolute squares tolerance
    best_cs = Vector{Float64}(undef, length(tbm))
    best_loss = Inf
    for t in 1:max_multistarts
        init_cs = randn(length(tbm))
        o = optimize(loss_closure, grad_loss_closure!, init_cs, optimizer, options)
        o.minimum > best_loss && continue # discard local optimization; not better globally
        if verbose
            println(
                "   Loss updated (trial $t): mean error = ",
                round(sqrt(o.minimum / (tbm.N * length(ks))); sigdigits = 3),
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