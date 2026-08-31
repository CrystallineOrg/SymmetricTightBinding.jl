# [Designing models to an objective](@id design)

[Fitting](@ref fitting) answers the question "which hopping amplitudes reproduce *this* band structure?". Often, though, there is no reference band structure to aim at -- only a property that the bands ought to have: a degeneracy at a particular energy, a gap at a particular filling, a flat band, a particular Chern number. These are optimization problems too; they merely swap the least-squares loss for a loss that measures the desired property.

SymmetricTightBinding.jl provides one such *design* objective out of the box, [`isolate_irrep`](@ref), and exposes the pieces it is assembled from so that other objectives can be posed in the same way.

!!! note "Optim.jl extension"
    `isolate_irrep` is implemented as a package extension and only becomes available once [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl) is loaded (`using Optim`). The objective itself ([`IrrepIsolationObjective`](@ref)) lives in the main package and needs no optimizer to evaluate.

## An isolated three-fold degeneracy

Consider space group 224 (*Pn*-3*m*) and the band representations (4b|A₂g) ⊕ (2a|A₂), giving a six-band model:

```@example design
using Crystalline, SymmetricTightBinding

sgnum = 224
brs = calc_bandreps(sgnum)
cbr = @composite brs[18] + brs[23]
tbm = tb_hamiltonian(cbr, [[0, 0, 0], [1, 0, 0]]) # on-site + nearest neighbour
SymmetryVector(cbr)
```

At Γ, the six bands split as 2Γ₂⁺ + Γ₄⁺ + Γ₁⁻: a single three-dimensional irrep, Γ₄⁺, and three singlets. Γ₄⁺ is thus a symmetry-enforced three-fold degeneracy for *any* choice of amplitudes -- but a generically uninteresting one, since the singlet bands will usually cross its energy elsewhere in the Brillouin zone.

What we want instead is for the degeneracy to be *energy-isolated*: writing ``E_t`` for the energy of the Γ₄⁺ triplet, the constant-energy surface ``E = E_t`` should meet the bands nowhere but at the degeneracy point itself. The triplet's *own* three bands are the ones that usually spoil this: they emanate from ``E_t`` at Γ and, left to themselves, wander back across it elsewhere in the zone. Isolation therefore requires them to split -- two of them settling below ``E_t`` and one above, or vice versa -- and to stay put. Whether amplitudes achieving this exist at all is a question about this particular band representation; that they do is what we set out to show here, constructively.

```@example design
using Optim, Random

Random.seed!(1)
ptbm = isolate_irrep(tbm, "Γ₄⁺"; Nk = 6, max_multistarts = 200)
nothing # hide
```

This search takes several minutes: design problems are markedly harder than fits, and the default of `max(60, 10length(tbm))` starts is not always enough -- with 50 starts, this one finds nothing at all. Pass `verbose = true` to watch it work.

## Verifying the outcome

Nothing about the search *guarantees* isolation: it is a nonconvex optimization over a **k**-mesh, so the result must be checked independently -- and on a finer mesh than the one it optimized over. That is what [`isolation_report`](@ref) is for:

```@example design
report = isolation_report(ptbm, "Γ₄⁺"; Nk = 24)
```

The triplet is found at bands 3:5, with `report.split == 2` of its bands below ``E_t`` and one above -- so ``E_t`` sits at a filling of four, with two singlets and two triplet bands beneath it. That the filling had to be *even* is forced by symmetry alone: at both M and X every irrep of this model is two-dimensional, so the bands come in degenerate pairs there and no odd filling can be gapped. Note also that the loss is invariant under ``H \to -H``, so a 2-below/1-above solution and its 1-below/2-above mirror are equally good; which appears depends on the seed.

Plotting the bands along a high-symmetry path, with ``E_t`` marked, shows what the optimizer settled on:

```@example design
using Brillouin, GLMakie

kpi = interpolate(irrfbz_path(sgnum, directbasis(sgnum, Val(3))), 100)
faxp = plot(
    kpi, spectrum(ptbm, kpi), fill(report.Et, length(kpi), 1);
    color = [:black, :crimson], linewidth = [3, 1], linestyle = [nothing, :dash],
    label = ["Bands", rich("E", subscript("t"))]
)
axislegend(; framevisible = false, position = :rb)
faxp # hide
```

The dashed line threads the band structure without touching it anywhere except at Γ, where the three bands of the triplet meet it -- two arriving from below, one from above. This is the whole point of the exercise: ``E_t`` is a chemical potential at which the material is a semimetal with a single, symmetry-enforced three-fold node and no other states at the Fermi level.

## How the objective is posed

The formulation that makes this tractable is to state isolation in terms of *sorted* band indices. Let the multiplet occupy bands `m:m+d-1` at Γ. Demanding that `p` of its bands fall below ``E_t`` and `d-p` above is then the same as demanding a clean separation at the filling ``\nu = m+p-1``:

```math
E_\nu(\mathbf{k}) < E_t < E_{\nu+1}(\mathbf{k})
\qquad \text{for every } \mathbf{k} \neq \mathbf{k}_0,
```

Because the bands are energy-sorted, this single pair of inequalities per **k**-point implies the separation of every other band for free: no band tracking, no connectivity analysis. The split must be proper, ``1 \leq p \leq d-1``; ``p = 0`` or ``p = d`` would leave ``E_t`` at a band edge of the multiplet, touched from one side along its entire constant-energy surface. Unless `split` fixes `p`, the objective minimizes over the admissible values of `p` as well as over the amplitudes.

!!! note "The excluded neighbourhood of the node"
    The requirement above provably fails as ``\mathbf{k} \to \mathbf{k}_0``, where the multiplet's bands converge on ``E_t``, so **k**-points within `kexclude` of ``\mathbf{k}_0`` (in reduced coordinates) are exempted from it. They are *not* dropped altogether: the bands *outside* the multiplet stay a finite distance from ``E_t`` right up to the degeneracy, and pushing them away there is exactly what keeps the node clean. Every sampled **k**-point therefore also carries ``E_{m-1}(\mathbf{k}) < E_t < E_{m+d}(\mathbf{k})`` -- "no more than the multiplet's own ``d`` bands may reach ``E_t``" -- so each contributes at most four constraints.

    This makes `kexclude` a genuine modelling choice rather than a numerical detail: it sets the scale below which the *node* may be approached, and hence what the attained gap means. Sampling `ks` along high-symmetry lines that keep clear of ``\mathbf{k}_0`` -- a decimated `irrfbz_path` from Brillouin.jl, say -- is a workable and often more economical alternative to a uniform mesh, since band extrema overwhelmingly sit on those lines; the trade is that the outside-band constraints then go unsampled near ``\mathbf{k}_0``.

Each constraint is a signed gap ``\delta_j`` that we want positive and large. [`IrrepIsolationObjective`](@ref) therefore minimizes a smooth (log-sum-exp) surrogate for the negated smallest of them,

```math
L = -\text{softmin}_\beta\big(\delta_j/\sigma\big)
  = \beta^{-1}\log \sum_j e^{-\beta\delta_j/\sigma},
```

which -- exactly as the intuition that one should "push the nearby bands away hardest" suggests -- weights each constraint by ``e^{-\beta\delta/\sigma}``, so that only the near-touching ones matter. The normalization ``\sigma`` is the spectrum's root-mean-square energy spread over `ks`, without which the loss could be improved indefinitely by simply scaling every amplitude up.

Three implementation points are worth spelling out:

- **Locating the multiplet.** Which bands carry Γ₄⁺ is not fixed: it changes as the amplitudes move. The objective therefore identifies it afresh at every evaluation, by grouping the degenerate bands at the irrep's **k**-point and matching their summed symmetry eigenvalues against the irrep's characters ([`locate_multiplet`](@ref)). This is cheap, because the matrices converting eigenvectors into symmetry eigenvalues depend only on the band representation and on **k**, never on the amplitudes, and are cached in the [`IrrepTarget`](@ref).
- **Gradients.** The loss has an analytic gradient: ``\partial_c\delta_j`` follows from the Feynman–Hellmann theorem (cf. [`energy_gradient_wrt_hopping`](@ref)), and ``\partial_c\sigma`` in closed form from the model's spectral moments, without diagonalizing anything. `isolate_irrep` accordingly defaults to a first-order optimizer (`LBFGS()`); no Hessian is provided.
- **Gauge fixing.** The loss is invariant under both an overall rescaling of the amplitudes and an overall energy shift. Two small quadratic penalties (`λ_scale`, `λ_shift`) pin those two flat directions, which changes which representative of a scale-and-offset family is returned, but not which models are optimal.

The search itself is the same moment-seeded, basin-hopping [`multistart_fit`](@ref) engine that [`fit`](@ref) uses; only the loss differs. One default does change: the final polish runs at `polish_options = Optim.Options(g_abstol=1e-6, f_reltol=1e-9)` rather than Optim.jl's own tolerances, which a soft-minimum loss cannot meet -- left at the defaults, the polish runs to its 1000-iteration cap and can cost more time than the entire search before it, for a fraction of a percent of loss.

### Warm starts

The landscape changes only gradually with the **k**-sampling, with `kexclude`, and with a `lasso` weight, so a solution found under one setting is usually an excellent starting point for the next. `init` takes either a vector of amplitudes or a whole [`ParameterizedTightBindingModel`](@ref) to lift them from, and seeds the deterministic first trial with it; later trials still hop around it, so the search is biased rather than confined. A cheap coarse search followed by a warm-started refinement against a denser mesh is generally a far better use of a time budget than one long cold run:

```@example design
Random.seed!(2)
ptbm_fine = isolate_irrep(tbm, "Γ₄⁺"; Nk = 12, max_multistarts = 20, init = ptbm)
isolation_report(ptbm_fine, "Γ₄⁺"; Nk = 24)
```

### Sparser models

As with [`fit`](@ref), a `lasso` weight adds an ``\ell_1`` penalty on the amplitudes, here normalized by the number of terms and by ``\sigma`` so that it stays scale-invariant (an unnormalized penalty would be minimized by shrinking every amplitude at once). It asks for the *simplest* model that still isolates the multiplet, and, as always, trades attained gap for sparsity. Switching the penalty on is itself a gradual change of landscape, so it is a natural place to warm-start from the model we already have:

```@example design
Random.seed!(1)
ptbm_sparse = isolate_irrep(tbm, "Γ₄⁺"; Nk = 6, max_multistarts = 30, lasso = 0.4,
                            init = ptbm)
report_sparse = isolation_report(ptbm_sparse, "Γ₄⁺"; Nk = 24)
nonzero(p) = count(>(1e-3), abs.(p.cs))
(dense  = (terms = nonzero(ptbm),        relgap = round(report.relgap; digits = 3)),
 sparse = (terms = nonzero(ptbm_sparse), relgap = round(report_sparse.relgap; digits = 3)))
```

Here two of the eighteen terms are switched off for a relative gap essentially unchanged in the third digit -- a modest gain, and a reminder that, exactly as for [`fit`](@ref), `lasso` is a blunt instrument: which terms it prunes depends on ``\lambda_1``, on the warm start, and on which local minimum the search settles in, and the count of surviving terms need not decrease monotonically with ``\lambda_1``. Treat it as a knob worth sweeping, not a model-selection procedure.

!!! warning "The objective is a proxy, not a proof"
    Maximizing the minimum gap over finitely many **k**-points is a surrogate for a statement about the whole Brillouin zone, and the search is a heuristic over a nonconvex landscape. A negative result means only that this model, sampling, and search failed to find isolation -- not that it is impossible. Always confirm a positive result with [`isolation_report`](@ref) on a finer sampling than the search used, and preferably also visually.

    The sharpest version of this trap is a **k**-sampling that misses a high-symmetry point. A symmetry-enforced degeneracy there can rule out a filling outright, and a mesh that skips it will happily report a gap that symmetry forbids -- as `Nk = 8` does for the hexagonal model of the previous section, whose K point sits at ``(\tfrac13, \tfrac13)``. Prefer an `Nk` divisible by 6, so that both the ½- and ⅓-type points are sampled.

## Posing your own objective

The seam is deliberately narrow: any function `_fgh!(F, G, H, cs)` following Optim.jl's convention -- return the loss when `F` is non-`nothing`, fill `G` in place when it is not -- can be handed to [`make_fit_objective`](@ref) and then to [`multistart_fit`](@ref). Because a design problem has no reference spectrum to take moments of, the exploration scales for the multi-start search instead come from a reference-free `SpectralMoments`, which asks only for the energy `scale` of the spectra to explore.

As a minimal illustration, here is a wholly different objective, on a wholly different model: maximize the *indirect* gap at filling one -- with no symmetry input whatsoever -- for a three-band plane group 17 (*p*6*mm*) model. The objective supplies no derivatives at all, so it is minimized with a derivative-free optimizer:

```@example design
using LinearAlgebra: eigvals!
using SymmetricTightBinding: spectralmoments

brs₂ᴰ = calc_bandreps(17, Val(2))
cbr₂ᴰ = @composite brs₂ᴰ[12] + brs₂ᴰ[8]                # (1a|E₁) ⊕ (1a|A₁): 3 bands
tbm₂ᴰ = tb_hamiltonian(cbr₂ᴰ, [[0, 0], [1, 0]])
cache = TightBindingCache(tbm₂ᴰ, uniform_kmesh(Val(2), 12))

function gap_loss(F, G, H, cs)
    (isnothing(G) && isnothing(H)) || error("no derivatives are provided")
    Em = reduce(hcat, eigvals!(cache(cs, κ)) for κ in eachindex(cache.ks))
    Ē = sum(Em) / length(Em)
    σ = sqrt(sum(abs2, Em) / length(Em) - Ē^2)
    return -(minimum(Em[2, :]) - maximum(Em[1, :])) / σ # −(indirect gap)/σ
end

Random.seed!(1)
cs, loss = multistart_fit(
    make_fit_objective(gap_loss), spectralmoments(cache; scale = 1.0);
    optimizer = NelderMead(), max_multistarts = 25, tol = -Inf, polish = false
)
-loss # the largest attained indirect gap, in units of the spectrum's RMS spread
```

Note the `tol = -Inf`: `multistart_fit` returns early once the loss drops to `tol`, whose default of `0.0` suits a least-squares fit but would stop a design search the moment it turned negative.

```@example design
kpi₂ᴰ = interpolate(irrfbz_path(17, directbasis(17, Val(2))), 100)
plot(kpi₂ᴰ, spectrum(tbm₂ᴰ(cs), kpi₂ᴰ))
```

## Design API

```@docs
isolate_irrep
```

([`IrrepIsolationObjective`](@ref), [`isolation_report`](@ref), [`IrrepTarget`](@ref), [`locate_multiplet`](@ref), and [`uniform_kmesh`](@ref) are documented in the [API](@ref) reference.)
