# SymmetricTightBinding.jl

SymmetricTightBinding.jl enables the construction and manipulation of tight-binding models. The main novelty - and principal strength - of SymmetricTightBinding.jl is that each model is associated with, and specified by, a set of band representations.

Put more simply, SymmetricTightBinding.jl will automatically parameterize all possible tight-binding Hamiltonians that are compatible with a selection of orbitals with specified local symmetries (i.e., transforming as specific site symmetry irreps), each situated at specified positions in the unit cell (i.e., at specific Wyckoff positions).

The underlying physics is that the Bloch Hamiltonian of a Wannierizable set of bands must transform under under a site-symmetry induced representation (also called band representation) $D(g)$ for operations $g$ in the associated space group. That is, the Bloch Hamiltonian $\mathbf{h}(\mathbf{k})$ is symmetric in the sense:

```math
\mathbf{D}_\mathbf{k}(g) \mathbf{h}_{\mathbf{k}} \mathbf{D}_\mathbf{k}^\dagger(g) = \mathbf{h}_{g\mathbf{k}}
```

where $\mathbf{D}_{\mathbf{k}}(g)$ is the (momentum-)block-diagonal part of the Fourier transformed band representation $D(g)$.

## Installation

The package is registered in the Julia General registry and can be installed from the `pkg>` command line (entered by pressing `]` in the Julia REPL):
```jl
pkg> add SymmetricTightBinding, Crystalline
```

SymmetricTightBinding.jl is designed to work as a companion package to [Crystalline.jl](https://github.com/thchr/Crystalline.jl); so we add Crystalline.jl in the above as well.

## Tutorial

As a first step, we load both Crystalline.jl and SymmetricTightBinding.jl into our current Julia session:

```@example basic-use
using Crystalline, SymmetricTightBinding
```


As our first example, we'll build the tight-binding model of graphene. Once we've done that, we'll explore how to create related tight-binding models in the same symmetry setting.

### Graphene

Graphene is a two-dimensional material, with two carbon atoms arranged in a honeycomb lattice. For the present purposes, the important aspects is its crystal symmetry: the lattice has 6-fold rotation symmetry and associated in-plane mirror symmetry. In the language of crystallography, its symmetry is that of [plane group](https://en.wikipedia.org/wiki/Wallpaper_group) [*p*6*mm*](https://en.wikipedia.org/wiki/Wallpaper_group#Group_p6m_(*632)) (here, specified by its Hermann-Mauguin label).
This plane group has a conventional numbering assigned---namely, plane group ⋕17. The mapping between the Hermann-Mauguin label and the conventional number can e.g., be determined using Crystalline.jl's `iuc` or looked up in online tables, such as the [Bilbao Crystallographic Server](http://www.cryst.ehu.es/).

Using Crystalline, we can build the maximal *band representations* in plane group ⋕17:

```@example basic-use
sgnum = 17 # space group number of p6mm
brs = calc_bandreps(sgnum, Val(2)) # `Val(2)` specifies the dimensionality (here, 2D)
```

The top row of the output lists the possible positions that a symmetrically placed orbital can reside, specified as a [Wyckoff position](https://en.wikipedia.org/wiki/Wyckoff_positions) label (e.g., 1a, 2b, 3c). In the second row, the possible local symmetry that an orbital placed there can have (e.g., A₁, A₂, B₁, …) are listed, specified in [Mulikken notation](https://en.wikipedia.org/wiki/List_of_character_tables_for_chemically_important_3D_point_groups). The remaining rows contain information about the projection of each band representation to band symmetries at high-symmetry **k**-points and which is not needed in the present context.

Graphene's two *p*<sub>*z*</sub> orbitals sit at the 2b Wyckoff position: though odd (i.e., changing sign) under mirror in the out-of-plane direction, the *p*<sub>*z*</sub> orbital is even (i.e., invariant) under all in-plane symmetries (rotations and mirrors). The associated site-symmetry irrep is the A₁ site-symmetry irrep of the 2b Wyckoff position. In the above tables, this is the fifth column of `brs`, which we may select by:

```@example basic-use
brs[1]
```

To construct a tight-binding model, we must construct a `CompositeBandRepresentation`: this is necessary because we may generally be interested in building models for multiple band representations (say, for placing orbitals at multiple distinct Wyckoff positions). We can construct such a composite representation by `@composite a*brs[i] + b*brs[j] + …` which will contain `a` times the `brs[i]` band representation and so on. Here, we just need the (2b|A₁) representation once:

```@example basic-use
cbr = @composite brs[5]
```

With this in hand, we can finally use SymmetricTightBinding.jl. In particular, we may use [`tb_hamiltonian`](@ref). First, we build the nearest-neighbor tight-binding Hamiltonian:

```@example basic-use
tbm = tb_hamiltonian(cbr)
```

The output lists the "basis terms" of the tight-binding Bloch Hamiltonian, each implicitly parameterized by a free on-site energy or hopping amplitude. The notation `𝕖(δᵢ)` is introduced for brevity, a short-hand for the complex momentum-dependent exponential $\mathrm{e}^{\mathrm{i}\mathbf{k}\cdot\boldsymbol{\delta}_i}$. Here $\boldsymbol{\delta}_i$ denotes a hopping vector; in turn, these vectors are expressed above as `δᵢ`, given in the basis of the primitive direct lattice $\{\mathbf{a}_i\}$. I.e., a term like `δ₁ = [-1/3, 1/3]` really means $\boldsymbol{\delta}_1 = -\tfrac{1}{3}\mathbf{a}_1 + \tfrac{1}{3}\mathbf{a}_2$.


### Visualization

It's often helpful to visualize model terms, and SymmetricTightBinding.jl facilitates this via a Makie extension. We may e.g., use GLMakie.jl to plot the second tight-binding term:

```@example basic-use
using GLMakie
Rs = directbasis(sgnum, Val(2)) # a direct lattice basis, to allow plotting in a Cartesian setting
plot(tbm[2], Rs)
```

Here, red markers indicate "source" sites while blue markers indicat "drains"; electrons hop from sources to drains, as also indicated by the arrowheads. The visualization (and the internal representation of `tbm`) includes only the hoppings for a *single* unit cell, such that tiling unit cells do not lead to counting hoppings multiple times.

We might want to go beyond nearest-neighbor in our tight-binding model. To do so, we must provide `tb_hamiltonian` with a second argument that gives a set of possible direct-lattice vector separations of sources and drains (in addition to an intra-lattice term). It is enough to include a representative direct lattice vector; if e.g., `[1,0]` and `[0,1]` are symmetry-related, the latter will be automatically included by providing the former. 
For the graphene example, we might include direct lattice separations `[0,0]` (default, if a second argument is not provided) and `[1,0]`:

```@example basic-use
tbm = tb_hamiltonian(cbr, [[0,0], [1,0]])
```

And we now have four terms. We can visualize `tbm[3]` and `tbm[4]` as before, individually, or we can visualize all terms at once:

```@example basic-use
plot(tbm, Rs)
```

### Model evaluation & band structures

To evaluate the tight-binding model, we must specify a set of (real) hopping amplitudes. To associate coefficients `c₁, c₂, c₃, c₄` to each of the basis terms of the model `tbm`, we can invoke it as a functor to create a [`ParameterizedTightBindingModel`](@ref):

```@example basic-use
cs = [0, 1, 0, 0] # set c₂ = 1 and all other cᵢ to zero
ptbm = tbm(cs)
```

We can then evaluate the parameterized model `ptbm` at **k**-point, again using `ptbm` as a functor:

```@example basic-use
k = ReciprocalPoint(1/2, 0) # the M point
h = ptbm(k)
```

Note that the **k**-point coordinates must be given in the basis of the primitive reciprocal lattice vectors $\{\mathbf{b}_i\}$ (the dual lattice to $\{\mathbf{a}_i\}$), i.e., the `k` variable above corresponds to the point $\mathbf{k} = \tfrac{1}{2}\mathbf{b}_1 + 0\mathbf{b}_2$.

We will usually be more interested in the overall behavior of the model across the Brillouin zone than its behavior at any single **k**-point. E.g., we might be interested in the band structure along high-symmetry lines of the Brillouin zone. To quickly build such a path, we leverage [Brillouin.jl](https://github.com/thchr/Brillouin.jl)'s `irrfbz_path`:
```@example basic-use
using Brillouin
kp = irrfbz_path(sgnum, Rs)
kpi = interpolate(kp, 200) # aim for 200 interpolations points
```

Next, to obtain the band structure along the interpolated **k**-path, we use SymmetricTightBinding.jl's `spectrum` function and plot the result using the Brillouin.jl's Makie.jl extension:

```@example basic-use
Es = spectrum(ptbm, kpi); # a 200×2 Matrix
plot(kpi, Es)
```

### Band symmetry

Since the theory behind SymmetricTightBinding.jl is anchored in symmetry analysis, the package also provides several tools to analyze band symmetry.

For instance, we can label the previously constructed band structure with the little group irrep labels at high-symmetry **k**-points:
```@example basic-use
plot(kpi, Es; annotations=collect_irrep_annotations(ptbm))
```

Similarly, we can analyze the compatibility respecting bands contained in `ptbm` via [`collect_compatible`](@ref). Here, since we our model contains only a single band representation - and one which is intrinsically connected - such an analysis has only possible answer (the connected band representation itself):

```@example basic-use
collect_compatible(ptbm)
```

We can easily set up a more interesting situation, however, by incorporating more band representations (i.e., more orbitals) into our model. E.g., below, we add three *s*-like orbitals placed at the 3c position (edges of the hexagonal unit cell) to the usual graphene model, and pick a reasonably large hybridization between the graphene and 3c orbitals:
```@example basic-use
cbr′ = @composite brs[5] + brs[1]
tbm′ = tb_hamiltonian(cbr′)
ptbm′ = tbm′([-4, -0, -0.1, 0.0, 1.0, -1.0, 1.0])
plot(kpi, spectrum(ptbm′, kpi); annotations=irrep_annotations(ptbm))
```
The band structures features two connected groups of bands. We can obtain the same result (via a compatibility-analysis involving only the high-symmetry **k**-points) via `collect_compatible`:
```@example basic-use
collect_compatible(ptbm′)
```

!!! todo
    I suspect that the above is an interesting example - I *think* it should be an example where the connected band groupings do not actually decompose to the original EBRs - but the case is clearly bugged presently (no doubt due to issue CrystallineOrg/SymmetricTightBinding#65).