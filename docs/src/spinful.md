# [Spinful models](@id spinful)

So far, we have built models from *spinless* (i.e., integer spin) band representations. Spin–orbit coupled electrons instead require *spinful* (i.e., half-integer spin) band representations that involve the double space group and its representations. In Crystalline.jl, these are obtained by passing `spinful = Val(true)` to `bandreps`; from there, building a model works as usual, via [`tb_hamiltonian`](@ref).

Relative to the spinless setting, the main new ingredient is that time-reversal symmetry, ``\mathcal{T}``, squares to ``\mathcal{T}^2 = −1`` for spinful particles (as opposed to ``\mathcal{T}^2 = +1`` for spinless particles), which pairs states into Kramers partners. Below, we first explain how Kramers partners appear in a model, using the simplest possible example, and then build a spin–orbit coupled model of graphene.

!!! note "Spinful models are currently 3D only"
    Crystalline.jl does not yet provide spinful band representations in 1D and 2D. For now, a 2D model can instead be built from a 3D space group with the same in-plane symmetry, as we do for graphene below.

## Kramers partners as pseudospin

We start with a single orbital in the unit cell of the space group *P*1 (⋕1), which has no symmetry except lattice translations and time-reversal symmetry:

```@example spinful
using Crystalline, SymmetricTightBinding
brs = bandreps(1; spinful = Val(true)) # spinful band representations of P1 (⋕1)
pin_free!(brs, [1 => [0, 0, 0]])       # place the free 1a position at the origin
cbr = @composite brs[1]                # (1a|AˢAˢ), 2 bands
```

Since there are no spatial symmetries, the 2 bands of (1a|AˢAˢ) simply correspond to the spin-up and spin-down states of the orbital. The superscript ˢ denotes a spinful irrep, and "AˢAˢ" denotes a pair of 1D irreps "Aˢ" paired by time reversal.

The two orbitals are Kramers partners, i.e., time-reversal conjugates. Time reversal acts on them via ``\mathcal{T} = \mathrm{i}\sigma_y K`` (with ``K`` denoting complex conjugation), exactly as on the up and down states of a spin-½. We can think of them as pseudospin ↑ and ↓ states; when, as in all examples on this page, the pair comes from a single spinless orbital, they are just spin ↑ and ↓. Crystalline.jl returns the unitary part of ``\mathcal{T} = \boldsymbol{\Gamma} K``, i.e., ``\boldsymbol{\Gamma} = \mathrm{i}\sigma_y``, via `timereversal_unitary`[^1]:

```@example spinful
timereversal_unitary(brs[1].siteir)
```

[^1]: More generally, for a site with a 2*n*-dimensional site-symmetry irrep, orbitals *i* and *i*+*n* are Kramers partners, and `timereversal_unitary` returns ``\boldsymbol{\Gamma} = \mathrm{i}\sigma_y \otimes \boldsymbol{1}_n``. For spinless particles, ``\boldsymbol{\Gamma} = \boldsymbol{1}``.

Next, we build a model with on-site terms and hoppings along the first lattice vector ``\mathbf{a}_1``:

```@example spinful
tbm = tb_hamiltonian(cbr, [[0, 0, 0], [1, 0, 0]])
```

Terms 1 and 2 are the on-site energy and an ordinary, spin-independent hopping: they are diagonal in pseudospin and equal for ↑ and ↓. Terms 3--5 are spin-dependent hoppings, i.e., spin–orbit coupling: they are proportional to ``\sin(2\pi k_1)`` times ``\sigma_y``, ``\sigma_z``, and ``\sigma_x``, respectively. Time-reversal symmetry, ``\mathbf{H}(-\mathbf{k}) = \boldsymbol{\Gamma} \mathbf{H}^*(\mathbf{k}) \boldsymbol{\Gamma}^\dagger = \mathrm{i}\sigma_y \mathbf{H}^*(\mathbf{k}) (\mathrm{i}\sigma_y)^\dagger``, requires spin-dependent terms to be odd in **k**, and spin-independent terms to be even.

As a result, the spin–orbit coupling terms vanish at time-reversal invariant momenta (TRIMs), i.e., at **k**-points where ``\mathbf{k} \equiv -\mathbf{k}``, such as ``k_1 = 0`` and ``k_1 = 1/2``. This is simply Kramers' theorem: at TRIMs, every band is (at least) doubly degenerate:

```@example spinful
ptbm = tbm([0, 1, 0.3, 0.2, 0.1])
ks = range(-1/2, 1/2, 201)
kvs = [[k₁, 0, 0] for k₁ in ks]
Es = spectrum(ptbm, kvs)

using GLMakie
update_theme!(linewidth = 4) # hide
axis = (; xlabel = rich("k", font=:italic) * subscript("1"), ylabel = "Energy", xticks = [-1/2, 0, 1/2])
lines(ks, Es[:, 1]; color = :royalblue, axis)
lines!(ks, Es[:, 2]; color = :royalblue)
current_figure() # hide
```

### Adding inversion symmetry

If we add inversion symmetry, moving to the space group *P*``\overline{1}`` (⋕2), the spin–orbit coupling terms of this model are forbidden:

```@example spinful
brs′ = bandreps(2; spinful = Val(true))
cbr′ = @composite brs′[end] # (1a|AgˢAgˢ), 2 bands
tb_hamiltonian(cbr′, [[0, 0, 0], [1, 0, 0]])
```

Inversion acts trivially on the (1a|AgˢAgˢ) orbitals (they are even, "g", under inversion), so inversion requires ``\mathbf{H}(-\mathbf{k}) = \mathbf{H}(\mathbf{k})``. Combined with time reversal, ``\mathbf{H}(-\mathbf{k}) = \mathrm{i}\sigma_y \mathbf{H}^*(\mathbf{k}) (\mathrm{i}\sigma_y)^\dagger``, this requires ``\mathbf{H}(\mathbf{k}) = \sigma_y \mathbf{H}^*(\mathbf{k}) \sigma_y``. Writing ``\mathbf{H}(\mathbf{k}) = h_0(\mathbf{k})\boldsymbol{1} + \mathbf{h}(\mathbf{k})\cdot\boldsymbol{\sigma}`` and using ``\sigma_y \boldsymbol{\sigma}^* \sigma_y = -\boldsymbol{\sigma}``, we find ``\mathbf{h}(\mathbf{k}) = \mathbf{0}``: so ``\mathbf{H}(\mathbf{k}) = h_0(\mathbf{k})\boldsymbol{1}``, and the two bands are degenerate at every **k**-point.

This holds for any spinful model with both inversion and time-reversal symmetry, and every band is then doubly degenerate at every **k**-point. Intuitively, inversion maps a state at (**k**, ↑) to (−**k**, ↑) (spin is unchanged by inversion), and time reversal maps that to (**k**, ↓): so the two have the same energy. More precisely, the combination of inversion and time reversal is an antiunitary symmetry that maps **k** to itself and squares to −1, so Kramers' theorem applies at every **k**-point. Spin–orbit coupling is not forbidden by inversion in general (see also the Kane--Mele term below)[^2], but it cannot remove this requirement of double degeneracy at every **k**-point: that requires breaking inversion or time reversal.

[^2]: For example, a model with two orbitals at the 1a position, one even and one odd under inversion (`@composite brs′[end] + brs′[end-1]`, i.e., ``(1a|A_g^sA_g^s) + (1a|A_u^sA_u^s)``), has a nearest-neighbor spin–orbit term ``\propto \sin(2\pi k_1)`` that couples the two orbitals while flipping pseudospin. Still, every band remains doubly degenerate.

## Kane--Mele model

```@setup spinful
ENV["COLUMNS"] = 400 # print the wide terms below without truncation
```

As a more involved example, we next build a spinful version of graphene: the Kane--Mele model[^3]. The spinless graphene model of the [tutorial](@ref) is built from *p*<sub>*z*</sub> orbitals at the 2b Wyckoff position of plane group *p*6*mm* (⋕17), transforming as the site-symmetry irrep A<sub>1</sub>. Here, we use instead the 3D space group *P*6*mm* (⋕183), which has the same in-plane symmetry.[^4] Adding spin, the *p*<sub>*z*</sub>↑ and *p*<sub>*z*</sub>↓ orbitals transform jointly as the 2D site-symmetry irrep E₁ˢ:

```@example spinful
sgnum = 183
brs = bandreps(sgnum; spinful = Val(true))
```

The 2b position of *P*6*mm* is [1/3, 2/3, *z*], with a free *z*-coordinate, which we set to zero:

```@example spinful
pin_free!(brs, [3 => [0, 0, 0]]) # set the free `z` coordinate of `position(brs[3])` to zero
cbr = @composite brs[3]
```

The model has 4 orbitals, ordered by pseudospin and then by site:

```
orbital index:  1  2  3  4
   pseudospin:  ↑  ↑  ↓  ↓
         site:  1  2  1  2
```

Here, the pseudospin states are spin ↑ and ↓ along *z*. The Hamiltonian thus has a 2×2 block form in pseudospin: matrix elements ``H_{ij}(\mathbf{k})`` in diagonal blocks (``i, j ≤ 2`` or ``i, j > 2``) conserve pseudospin, while those in off-diagonal blocks flip it. We include on-site terms and hoppings up to next-nearest neighbors (the first 6 terms):

```@example spinful
tbm = tb_hamiltonian(cbr, [[0, 0, 0], [1, 0, 0]])[1:6]
```

Since there are no hoppings along *z*, the model describes decoupled graphene layers. The terms are:

1. On-site energy.
2. Nearest-neighbor hopping that flips pseudospin: the Rashba spin–orbit coupling. Its complex coefficients (e.g., 1.155 = 2/√3) encode the bond-direction dependence of the Rashba coupling.
3. Nearest-neighbor hopping that conserves pseudospin: the usual graphene hopping.
4. Next-nearest-neighbor hopping that conserves pseudospin, equal for ↑ and ↓: the usual next-nearest-neighbor hopping.
5. Next-nearest-neighbor hopping that flips pseudospin: a longer-range Rashba-like term.
6. Next-nearest-neighbor hopping that conserves pseudospin, with opposite signs for ↑ and ↓ (and for the two sites): the Kane--Mele spin–orbit coupling.

!!! note "Why does the usual hopping (term 3) carry factors of ±i?"
    Term 3 is the ordinary, spin-independent graphene hopping, but written in a basis where the ↑ and ↓ orbitals of site 2 carry extra phases relative to those of site 1. These phases come from how the orbitals of site 2 are defined: as the images of the orbitals of site 1 under a symmetry operation that maps site 1 to site 2. For our band representation, this operation is:
    ```@example spinful
    g = cosets(group(brs[3]))[2] # maps site 1, [1/3, 2/3, 0], to site 2, [2/3, 1/3, 0]
    ```
    The rightmost matrix is its SU(2) part, i.e., its action on spin: ``-\mathrm{i}\sigma_z``. The operation leaves the *p*<sub>*z*</sub> orbital itself unchanged, so the ↑ and ↓ orbitals of site 2 simply acquire phases −i and +i relative to those of site 1. Collecting these phases in a diagonal matrix ``\mathbf{U} = \mathrm{diag}(1, -\mathrm{i}, 1, +\mathrm{i})``, with one entry per orbital (recalling the orbital order above, (1↑, 2↑, 1↓, 2↓)), term 3 equals ``\mathbf{U}^\dagger [\boldsymbol{1}_2 \otimes \mathbf{H}^{\text{spinless}}(\mathbf{k})] \mathbf{U}``, where ``\boldsymbol{1}_2 \otimes \mathbf{H}^{\text{spinless}}(\mathbf{k})`` is the spinless graphene hopping of the [tutorial](@ref), "stacked" for ↑ and ↓ pseudospins, and where ``\mathbf{U}`` rotates to our orbital frame. Only matrix elements between site 1 and site 2 pick up ``\pm\mathrm{i}`` phases; since the phases amount to a choice of basis, they do not change the spectrum.

Time reversal acts on the full model by ``\mathbf{H}(-\mathbf{k}) = \boldsymbol{\Gamma}\mathbf{H}^*(\mathbf{k})\boldsymbol{\Gamma}^\dagger``, with ``\boldsymbol{\Gamma}`` a Kronecker product of ``\mathrm{i}\sigma_y`` over the pseudospin and the identity over the sites: here, for two sites, ``\boldsymbol{\Gamma} = \mathrm{i}\sigma_y \otimes \boldsymbol{1}_2``. (For a model with several band representations, ``\boldsymbol{\Gamma}`` is additionally block-diagonal over the featured band representations.) We can verify this for a random set of amplitudes:

```@example spinful
using LinearAlgebra
Γ = kron(timereversal_unitary(brs[3].siteir), I(2))
ptbm = tbm(randn(6))
k = [0.1, 0.2, 0.0]
H₊, H₋ = copy(ptbm(k)), copy(ptbm(-k)) # `copy`, since `ptbm` reuses its output buffer
H₋ ≈ Γ * conj(H₊) * Γ'
```

### Band structures

Since the model is independent of ``k_3``, we restrict the **k**-path to the ``k_3 = 0`` plane, via [Brillouin.jl](https://github.com/thchr/Brillouin.jl):

```@example spinful
using Brillouin
Rs = directbasis(sgnum)
kp = irrfbz_path(sgnum, Rs)
kp = KPath(kp.points, [[:Γ, :M, :K, :Γ]], kp.basis, kp.setting) # in-plane part of path
kpi = interpolate(kp, 200)
nothing # hide
```

With only the usual graphene hopping (term 3), we recover graphene's Dirac cones at K, now with every band doubly degenerate (gray). Adding the Kane--Mele term (term 6) opens a gap at K (red), turning the model into a quantum spin Hall insulator[^3]:

```@example spinful
Es_graphene = spectrum(tbm([0, 0, 1, 0, 0, 0]), kpi)    # nearest-neighbor coupling only
Es_KM       = spectrum(tbm([0, 0, 1, 0, 0, 0.05]), kpi) # + Kane-Mele SOC
plot(kpi, Es_graphene, Es_KM; color = [:gray, :crimson])
```

Both band structures remain doubly degenerate everywhere: terms 3 and 6 are also allowed in inversion-symmetric, free-standing graphene, so the argument from *P*``\overline{1}`` applies.
The Rashba term (term 2), in contrast, breaks the *z*-mirror (and hence inversion) and splits the bands -- except at the TRIMs Γ and M, where Kramers' theorem keeps them degenerate. At K, which is not a TRIM, the 4 bands split into two singlets (K₄ˢ and K₅ˢ) and a doublet (K₆ˢ):

```@example spinful
ptbm = tbm([0, 0.05, 1, 0, 0, 0.05]) # Kane-Mele and Rashba SOC
plot(kpi, spectrum(ptbm, kpi); annotations = collect_irrep_annotations(ptbm))
```

[^3]: C.L. Kane & E.J. Mele, *Quantum Spin Hall Effect in Graphene*, [Phys. Rev. Lett. **95**, 226801 (2005)](https://doi.org/10.1103/PhysRevLett.95.226801).

[^4]: Free-standing graphene is also symmetric under the mirror *z* → −*z*, giving the space group *P*6/*mmm* (⋕191). *P*6*mm* lacks this mirror, as e.g., for graphene on a substrate or in a perpendicular electric field. This is what allows the Rashba term below.

```@setup spinful
delete!(ENV, "COLUMNS")
```
