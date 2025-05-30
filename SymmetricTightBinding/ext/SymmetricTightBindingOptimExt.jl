module SymmetricTightBindingOptimExt

using SymmetricTightBinding
using Optim

# ---------------------------------------------------------------------------------------- #
# Define loss as sum of absolute squared error (MSE, up to scaling)

# MSE loss
function loss(Em_r, ks, tbm, cs)
    L = zero(Float64)
    for (Es_r, k) in zip(eachrow(Em_r), ks)
        Es = spectrum(tbm(cs), k)
        L += sum(abs2, (E_r - E for (E_r, E) in zip(Es_r, Es)))
    end
    return L
end

# gradient of MSE loss
function grad_loss!(Em_r, ks, tbm, cs, G = zeros(Float64, length(cs)))
    ptbm = tbm(cs)
    fill!(G, zero(eltype(G)))
    for (Es_r, k) in zip(eachrow(Em_r), ks)
        Es = spectrum(ptbm, k)
        ∇Es = energy_gradient_wrt_hopping(ptbm, k)
        for (E_r, E, ∇E) in zip(Es_r, Es, ∇Es)
            G .+= (E_r - E) .* ∇E
        end
    end
    G .*= -2
    return G
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

## Keyword arguments
- `optimizer` (default, `Optim.LBFGS()`): a local optimizer from Optim.jl, capable of
  exploiting gradient information.
- `options` (default, empty): a `Optim.Options(…)` structure of optimization options.
- `max_multistarts` (default, `150`): maximum number of multi-start iterations.
- `atol` (default, `1e-3`): threshold for early return, specifying the minimum required mean
  energetic error (averaged over bands and **k**-points).
- `verbose` (default, `false`): whether to print information on optimization progress.

## Example

As a synthetic example, we might use `fit` to recover the coefficients of a randomly
parameterized tight-binding model, using its spectrum sampled over 10 **k**-points:

```jldoctest
julia> using Crystalline, SymmetricTightBinding, Brillouin
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
function SymmetricTightBinding.fit(
    tbm::TightBindingModel{D},
    Em_r::AbstractMatrix{<:Real},
    ks::AbstractVector{<:SymmetricTightBinding.ReciprocalPointLike{D}};
    optimizer::Optim.FirstOrderOptimizer = LBFGS(),
    options::Optim.Options = Optim.Options(),
    max_multistarts::Integer = 150,
    atol::Real = 1e-3, # minimum threshold error, per k-point & per band, averaged over both
    verbose::Bool = false,
) where D

    # let-block-capture-trick to make absolutely sure we have no closure boxing issues
    loss_closure = let Em_r = Em_r, ks = ks, tbm = tbm
        cs -> loss(Em_r, ks, tbm, cs)
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
            "      `max_multistarts` exceeded: tolerance not met\n";
            color = :yellow,
        )
    end

    return tbm(best_cs)
end

end # module SymmetricTightBindingOptimExt