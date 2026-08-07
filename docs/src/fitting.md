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
    color = [:gray, :crimson], linestyle = [nothing, :dash], linewidth = [5, 3]
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
    color = [:gray, :crimson], linestyle = [nothing, :dash], linewidth = [5, 3]
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
    color = [:gray, :royalblue, :crimson],
    linestyle = [nothing, :dash, :dash], linewidth = [5, 3, 3],
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

If we want to *encourage sparsity* -- letting the fit decide which longer-range terms are actually (or merely primarily) needed -- we can set the `lasso` keyword to a positive value ``\lambda``. This adds an [``\ell_1`` penalty](https://en.wikipedia.org/wiki/Lasso_(statistics)) ``\lambda\sum_i|c_i|`` to the loss, which shrinks weakly-supported amplitudes and, when successful, drives them to zero outright.

The setting this targets is the one described above: data produced by something other than a tight-binding model -- a DFT or other wave-equation calculation -- for which no finite-range parameterization is exact. There are then no "true" amplitudes to recover, and the goal shifts to finding the *simplest* model that still tracks the data closely. Since we cannot know in advance how many hopping ranges the data will support, the practical approach is to include generously many and let the penalty prune those that do not earn their place.

!!! details "Example: pruning an over-specified model"
    As a synthetic stand-in for such data, we generate a reference from a 20-term model whose five leading amplitudes are ``\mathcal{O}(1)`` and whose remaining fifteen form a small, long-ranged tail. We then fit it with a 15-term model: fewer terms than the reference, but still far more than we would like to end up with.

    ```@example fitting
    Rs_20  = [[0,0], [1,0], [1,1], [1,2], [2,2], [0,3], [1,3], [2,3], [3,3], [1,4], [2,4]]
    tbm_20 = tb_hamiltonian(cbr, Rs_20)      # reference model
    tbm_15 = tb_hamiltonian(cbr, Rs_20[1:8]) # fitting model
    (length(tbm_20), length(tbm_15))
    ```

    ```@example fitting
    Random.seed!(2)
    cs_20 = vcat(randn(5), 0.08 .* randn(length(tbm_20)-5) .* [0.9^i for i in 1:length(tbm_20)-5])
    Em_20 = spectrum(tbm_20(cs_20), kpi)
    nothing # hide
    ```

    Fitting the 15-term model twice, with and without the penalty:

    ```@example fitting
    Random.seed!(3)
    ptbm_dense = fit(tbm_15, Em_20, kpi)              # unpenalized
    ptbm_lasso = fit(tbm_15, Em_20, kpi; lasso = 1.0) # ℓ₁-penalized
    round.([ptbm_dense.cs ptbm_lasso.cs]; digits = 3) # amplitudes, side by side
    ```

    The unpenalized fit puts amplitude on every term available to it, while the penalized fit switches a third of them off entirely:

    ```@example fitting
    (dense = count(>(1e-3), abs.(ptbm_dense.cs)),
     lasso = count(>(1e-3), abs.(ptbm_lasso.cs)))
    ```

    The spectral price for this is modest -- measured as an RMS energy error relative to the bandwidth of the reference bands:

    ```@example fitting
    bandwidth = maximum(Em_20) - minimum(Em_20)
    rel_rms(ptbm) = 100norm(spectrum(ptbm, kpi) - Em_20) / sqrt(length(Em_20)) / bandwidth
    (dense = round(rel_rms(ptbm_dense); sigdigits = 2), # RMS error, % of bandwidth
     lasso = round(rel_rms(ptbm_lasso); sigdigits = 2))
    ```

    Both models track the reference to well under 1% of its bandwidth, so the sparser one is arguably the better answer: it is very nearly as faithful, but involves only two thirds as many hopping terms. A visual comparison is instructive also:

    ```@example fitting
    faxp = plot(
        kpi, Em_20, spectrum(ptbm_dense, kpi), spectrum(ptbm_lasso, kpi),
        color = [:gray, :royalblue, :crimson],
        linestyle = [nothing, :dash, :dash], linewidth = [5, 3, 3],
        label = ["Reference", "Fit (dense)", "Fit (lasso)"]
    )
    axislegend(; framevisible = false, position = :cb)
    faxp # hide
    ```

!!! warning "`lasso` is a blunt instrument"
    The [Lasso](https://en.wikipedia.org/wiki/Lasso_(statistics)) penalty's sparsifying tendency is not systematic. Whether a given term is driven to zero depends on ``\lambda``, on the model, and on which local minimum the search happens to settle in -- and the number of surviving terms need not even decrease monotonically with ``\lambda``. Treat `lasso` as a knob worth exploring rather than a dependable model-selection procedure.

## An example with real data

As a demonstration on non-synthetic data, we fit the band structure of face-centered cubic lead (Pb) near the Fermi level, using the reference data of [Wannier90's tutorial 2](https://github.com/wannier-developers/wannier90/tree/develop/tutorials/tutorial02) (four bands along L–Γ–X–U–Γ, digitized from its `lead.pdf` figure and included [here](https://github.com/CrystallineOrg/SymmetricTightBinding.jl/blob/main/docs/src/assets/lead-wannier90-bands.dat)).

!!! details "Example: the band structure of lead"
    The four bands originate from the *sp*³ manifold at the 4a Wyckoff position of space group 225 (*Fm*-3*m*), which decomposes into the (4a|A₁g) ⊕ (4a|T₁ᵤ) band representations:

    ```@example fitting
    using DelimitedFiles

    sgnum = 225
    brs_pb = calc_bandreps(sgnum, Val(3))
    cbr_pb = @composite brs_pb[end-9] + brs_pb[end-2] # (4a|A₁g) ⊕ (4a|T₁ᵤ)
    ```

    The reference **k**-path is a uniform sampling of the L–Γ–X–U–Γ path, which we reconstruct explicitly (379 points in total, with 101, 116, 42, and 123 points per segment, end-points included):

    ```@example fitting
    using Crystalline: SVector

    kvs = Dict(:Γ => SVector(0.0, 0, 0),   :L => SVector(0.5, 0.5, 0.5),
               :X => SVector(0.5, 0, 0.5), :U => SVector(0.625, 0.25, 0.625))
    labs, Ns = [:L, :Γ, :X, :U, :Γ], [101, 116, 42, 123]
    ks = mapreduce(vcat, 1:4) do i
        seg = range(kvs[labs[i]], kvs[labs[i+1]], Ns[i])
        i == 4 ? seg : seg[1:end-1] # drop the end-point shared with the next segment
    end
    Gs = primitivize(dualbasis(directbasis(sgnum, Val(3))), centering(sgnum, 3))
    kpi_pb = KPathInterpolant([ks], [Dict(1 => :L, 101 => :Γ, 216 => :X, 257 => :U, 379 => :Γ)],
                              Gs, Ref(Brillouin.LATTICE))

    Em_pb = readdlm(joinpath(pkgdir(SymmetricTightBinding), "docs", "src", "assets",
                             "lead-wannier90-bands.dat"), Float64; comments = true)
    size(Em_pb) # 379 k-points × 4 bands (eV)
    ```

    We then build a model with hoppings along the [100] direction out to two cells and fit it. Of the resulting 12 terms, 3 carry negligible amplitude (identified by a `lasso`-penalized trial fit) and are dropped, leaving a 9-term model:

    ```@example fitting
    tbm_pb = tb_hamiltonian(cbr_pb, [[0,0,0], [1,0,0], [2,0,0]])[[1:6..., 8, 9, 11]]
    Random.seed!(1)
    ptbm_pb = fit(tbm_pb, Em_pb, kpi_pb; max_multistarts = 50)
    nothing # hide
    ```

    ```@example fitting
    E_F = 5.27 # Fermi level (eV)
    faxp = plot(
        kpi_pb, fill(E_F, length(kpi_pb), 1), Em_pb, spectrum(ptbm_pb, kpi_pb);
        ylabel = "Energy (eV)", label = [rich("E", subscript("F")), "Reference", "Fit"],
        color = [:gray50, :gray, :crimson], linewidth = [1, 5, 3],
        linestyle = [:dash, nothing, :dash]
    )
    axislegend(; framevisible = false, position = (:center, :top))
    faxp # hide
    ```

    The agreement is good across the entire path, from just 9 amplitudes. Note, however, that the reference data is itself Wannier-interpolated, i.e., is effectively already a (long-ranged) tight-binding model: a genuinely *ab initio* reference would generally be harder to match this closely.

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
