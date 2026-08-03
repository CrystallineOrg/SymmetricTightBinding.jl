# [Fitting models to band structures](@id fitting)

A common task is to determine the hopping amplitudes of a symmetry-constrained tight-binding model so that its bands reproduce a *reference* band structure -- e.g., one obtained from a first-principles (DFT) calculation, an experiment, or a classical physics calculation to which a mapping is sought (e.g., photonic or acoustic crystals). SymmetricTightBinding.jl provides the [`fit`](@ref) function for this purpose, as an Optim.jl extension.

!!! note "Optim.jl extension"
    `fit` is implemented as a package extension and only becomes available once [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl) is loaded (`using Optim`).

Because [`tb_hamiltonian`](@ref) already fixes the *form* of every allowed Hamiltonian term (enforcing spatial symmetry, time-reversal, and hermiticity), fitting only has to determine a handful of real amplitudes -- one per symmetry-distinct hopping term. This makes the problem low-dimensional and the resulting model automatically symmetry-consistent, regardless of how noisy or incomplete the reference data is.

## A synthetic reference

To demonstrate the workflow, we generate a synthetic "reference" band structure from a known set of amplitudes and then try to recover it -- pretending, for the moment, that we only have the band energies. (In a real application, the reference energies might instead come from DFT, experiment, or some other wave-equation calculation.)

We base the model on graphene: the (2b|A₁) elementary band representation of plane group 17 (*p*6*mm*), with nearest- and next-nearest-neighbor hopping:

```@example fitting
using Crystalline, SymmetricTightBinding

sgnum = 17
brs = calc_bandreps(sgnum, Val(2))
cbr = @composite brs[5] # (2b|A₁)
tbm = tb_hamiltonian(cbr, [[0, 0], [1, 0], [1, 1]]) # on-site, NN, NNN
```

Next, we pick a set of "ground-truth" amplitudes and evaluate the associated band structure along a high-symmetry **k**-path, using [Brillouin.jl](https://github.com/thchr/Brillouin.jl). This spectrum will play the role of our reference data:

```@example fitting
using Brillouin

Rs = directbasis(sgnum, Val(2))
kpi = interpolate(irrfbz_path(sgnum, Rs), 100)

using Random
Random.seed!(2)
cs_ref = randn(length(tbm)) # ground-truth amplitudes (unknown to the fit)
ptbm_ref = tbm(cs_ref)

Em_ref = spectrum(ptbm_ref, kpi) # reference energies: a (100 × 2) matrix
nothing # hide
```

The reference energies `Em_ref[i, n]` give the energy of band `n` at the `i`th **k**-point, with bands sorted in ascending energy -- the format `fit` expects.

## Fitting

With the model and the reference in hand, fitting is a single call:

```@example fitting
using Optim # invoke to load extension
ptbm_fit = fit(tbm, Em_ref, kpi)
```

The returned [`ParameterizedTightBindingModel`](@ref) carries the fitted amplitudes. In this simple example, we can expect near-exact recovery of the original amplitudes:

```@example fitting
using LinearAlgebra: norm
norm(ptbm_ref.cs - ptbm_fit.cs)
```

We can confirm the recovery visually by overlaying the fitted bands (crimson, dashed) on the reference band structure (gray, solid):

```@example fitting
using GLMakie

Em_fit = spectrum(ptbm_fit, kpi)
plot(
    kpi, Em_ref, Em_fit;
    color = [:gray, :crimson], linestyle = [nothing, :dash], linewidth = [5, 2]
)
```

## Fitting to non-exactly representable data

Real band structures are usually never exactly mappable to finite-range tight-binding models. To emulate this, we can create a model with a few small long-range hopping terms and then attempt to fit its spectrum to a shorter-range model. Because the model form is symmetry-constrained, the fit remains rather good and simply returns a spectrally close symmetry-consistent model (in a least-squares sense).

```@example fitting
tbm_long = tb_hamiltonian(cbr, [[0, 0], [1, 0], [1, 1], [1,2], [2,2]]) # terms 6-8 are additional to onsite+NN+NNN
Random.seed!(2)
ptbm_long_ref = tbm_long(vcat(randn(5), randn(length(tbm_long)-5)*.05)) # small amplitudes in last three terms
Em_long_ref = spectrum(ptbm_long_ref, kpi)

ptbm_short_fit = fit(tbm, Em_long_ref, kpi) # fit to 5-term "short-range" model (onsite, NN, NNN)
Em_short_fit = spectrum(ptbm_short_fit, kpi)

plot(
    kpi, Em_long_ref, Em_short_fit;
    color = [:gray, :crimson], linestyle = [nothing, :dash], linewidth = [5, 2]
)
```

The fitted hopping amplitudes of the short-range model are in rough agreement with the ground-truth amplitudes of the long-range model, but not exact. The deviation is primarily a result of the short-range amplitudes compensating for the absence of longer-range terms.

```@example fitting
hopping_scale = sum(abs, ptbm_long_ref.cs[1:5]) / 5
(ptbm_short_fit.cs - ptbm_long_ref.cs[1:5]) ./ hopping_scale * 100 # relative deviation (%)
```

### An infinite-range reference

The example above is still representable in principle -- the reference merely has more terms than the fitting model. A sharper test is a reference that *no* finite-range model can represent. Consider a 1D chain with two orbitals per site, coupled across every distance with an exponentially decaying amplitude ``t_n = t_0 e^{-\gamma n}``. Summing the geometric series gives a closed-form Bloch Hamiltonian, ``H(k) = f(k)\sigma_x``, and hence exact bands ``\pm f(k)``:

```@example fitting
brs_1d = calc_bandreps(2, Val(1)) # 1D space group 2
cbr_1d = @composite brs_1d[1] + brs_1d[1] # two orbitals per site

tₙ(n, γ = 0.5, t₀ = 1.0) = t₀ * exp(-γ * abs(n)) # hopping amplitude across n cells
f(k, γ = 0.5, t₀ = 1.0) = t₀ * (1 - exp(-2γ)) / (1 - 2exp(-γ) * cospi(2k) + exp(-2γ)) # = ∑ₙ tₙ exp(i2πnk)

kpi_1d = interpolate(irrfbz_path(2, directbasis(2, Val(1))), 100)
Em_1d = (v = f.(only.(kpi_1d)); [-v v])
nothing # hide
```

We now fit truncations of this chain, keeping hoppings only out to `nmax` cells. `tb_hamiltonian` orders its terms block-major, so the two intra-orbital blocks come first and the inter-orbital couplings last; slicing off the former leaves exactly the terms of the model above:

```@example fitting
function truncated_chain(nmax)
    tbm = tb_hamiltonian(cbr_1d, [[n] for n in 0:nmax])
    return tbm[2(nmax + 1) + 1 : end] # keep only the inter-orbital couplings
end

ptbm_2 = fit(truncated_chain(2), Em_1d, kpi_1d)
ptbm_6 = fit(truncated_chain(6), Em_1d, kpi_1d)

faxp = plot(
    kpi_1d, Em_1d, spectrum(ptbm_2, kpi_1d), spectrum(ptbm_6, kpi_1d);
    color = [:gray, :cornflowerblue, :crimson],
    linestyle = [nothing, :dash, :dash], linewidth = [5, 2, 2],
    label = ["Reference", "Fit: 2 terms", "Fit: 6 terms"]
)
axislegend(; framevisible=false)
faxp # hide
```

Two cells (blue) are too short-ranged to resolve the sharp peak of ``f`` at Γ, and the bands acquire spurious oscillations -- and even a spurious touching -- across the rest of the zone. Six cells (crimson) already track the exact bands closely. Notably, the fit is never told that the amplitudes decay exponentially, yet recovers ``t_n = t_0e^{-\gamma n}`` to three digits:

```@example fitting
round.([ptbm_6.cs tₙ.(0:6)]; sigdigits = 3) # fitted amplitudes vs. exact tₙ
```


## Choosing the number of terms

The set of hopping ranges passed to [`tb_hamiltonian`](@ref) determines how many free amplitudes the fit has to work with. Too few, and the model cannot represent the reference; too many, and one runs the risk of overfitting (and, in addition, a more challenging fitting convergence due to a preponderance of local minima). A useful diagnostic is to fit models of increasing range and monitor the residual error:

```@example fitting
Em_fits = Matrix{Float64}[]
for range in ([[0, 0]], [[0, 0], [1, 0]], [[0, 0], [1, 0], [1, 1]])
    tbm_range = tb_hamiltonian(cbr, range)
    ptbm_range_fit = fit(tbm_range, Em_ref, kpi)
    Em_fit = spectrum(ptbm_range_fit, kpi)
    rms = norm(Em_fit - Em_ref) / sqrt(length(Em_ref))
    println("$(length(tbm_range)) terms: RMS error = ", round(rms; sigdigits = 3))
end
```

Here the reference was generated with next-nearest-neighbor hopping, so the error drops sharply once the model includes terms in the `[1, 1]` range and only then reaches (numerical) zero.

### Sparsity-encouraged fitting

If we want to *encourage sparsity* -- letting the fit decide which longer-range terms are actually (or merely primarily) needed -- we can set the `lasso` keyword to a finite value. This adds an ``\ell_1`` penalty on the amplitudes, tending to drive weakly-supported terms toward zero:

```@example fitting
ptbm_sparse = fit(tbm, Em_ref, kpi; lasso = 1e-1)
round.(ptbm_sparse.cs; sigdigits = 3)
```

(In this example every term genuinely contributes to the reference, so none is driven all the way to zero -- the penalty merely biases the amplitudes slightly. The effect is most useful when the model includes ranges that the data does not actually require.)

## Under the hood

`fit` performs a moment-seeded, basin-hopping multi-start minimization of a sorted-eigenvalue least-squares loss:
- the loss landscape of band-fitting problems is *funnel-like*, with hierarchies of progressively better band-assignment local minima, so the search perturbs the incumbent best fit ("hops") far more often than it restarts from scratch;
- the deterministic first trial, and the occasional fresh restarts, are seeded using the reference spectrum's first two moments (exploiting the linearity ``H(\mathbf{k}) = \sum_i c_i h_i(\mathbf{k})`` for cheap, sorting-free handles on the coefficient scale);
- second-order optimizers use a Gauss–Newton approximation of the Hessian, ``\nabla^2 F \approx 2\sum_{\mathbf{k}, n} \nabla E_n \nabla E_n^\top``, so the default `NewtonTrustRegion()` acts as a Levenberg–Marquardt least-squares solver.

The search engine is exposed separately as [`multistart_fit`](@ref), which minimizes any Optim.jl objective (built with [`make_fit_objective`](@ref)) against a moment structure. This lets related fitting problems with custom losses -- for example photonic band-fitting, which penalizes longitudinal modes differently -- reuse the same machinery.

For workloads that evaluate the same model over the same **k**-points many times with varying amplitudes (fitting, above all), the coefficient-independent term matrices ``h_i(\mathbf{k})`` are tabulated once up front in a [`TightBindingCache`](@ref); `fit` builds and shares one internally.

## Fitting API

```@docs
fit
multistart_fit
make_fit_objective
```

([`TightBindingCache`](@ref) is documented in the [API](@ref) reference.)
