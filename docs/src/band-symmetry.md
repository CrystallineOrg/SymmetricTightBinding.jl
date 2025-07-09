# Band symmetry and topology

Since the theory behind SymmetricTightBinding.jl is anchored in symmetry analysis, the package naturally provides several tools to analyze band symmetry, as well as to use this information to apply the frameworks of topological quantum chemistry and symmetry indicators to analyze band topology.

To explore these tools, we first re-build the graphene model previously explored in the [tutorial](@ref) section, instantiating the variables `tbm` (model), `ptbm` (coefficient-parameterized model), `kpi` (interpolated **k**-path, via Brillouin.jl), and `Es` (band structure of `ptbm`, evaluated over `kpi`).

!!! details "Setup details"
    ```@example band-symmetry
    using Crystalline, SymmetricTightBinding
    using Brillouin, GLMakie           # for k-space path and plotting
    sgnum = 17                         # plane group p6mm
    brs = calc_bandreps(sgnum, Val(2)) # band representations
    cbr = @composite brs[5]            # (2b|A₁) EBR
    tbm = tb_hamiltonian(cbr)          # tight-binding model (nearest neigbors)
    ptbm = tbm([0, 1])                 # zero self-energy, nonzero nearest-neighbor hopping
    Rs = directbasis(sgnum, Val(2))    # (conventional) direct lattice basis
    kp = irrfbz_path(sgnum, Rs)        # high-symmetry k-path
    kpi = interpolate(kp, 100)         # interpolated k-path
    Es = spectrum(ptbm, kpi)           # band structure over `kpi`
    nothing # hide
    ```

## Annotating little group irrep labels
To annotate a band structure plot with the little group irrep labels at high-symmetry **k**-points, we can use [`collect_irrep_annotations`](@ref) in combination with the `annotations` keyword argument of the Makie `plot` extension:
```@example band-symmetry
plot(kpi, Es; annotations = collect_irrep_annotations(ptbm))
```

## Collecting compatibility respecting band groups
Similarly, we can analyze the compatibility respecting bands contained in `ptbm` via [`collect_compatible`](@ref):

```@example band-symmetry
collect_compatible(ptbm)
```

`collect_compatible` returns a list of symmetry vectors, from lowest-energy band grouping to highest, each aggregating the symmetry content of a minimal set of compatibility-respecting bands.
Here, since our model contains only a single band representation -- which is additionally an intrinsically connected one -- such a list can have only one possible element: the only possible band groupings is the original band representation. We can verify this by comparing with the symmetry vector of the band representation used to build `tbm`:
```@example band-symmetry
SymmetryVector(CompositeBandRep(ptbm))
```

We can set up a more interesting situation by incorporating more band representations (i.e., more orbitals) into our model. E.g., below, we add three *s*-like orbitals placed at the 3c Wyckoff position (edges of the hexagonal unit cell; i.e., a kagome-like lattice) to the usual graphene model. First, we look at a situation without hybridization and with the bands of the two orbitals sets overlapping:
```@example band-symmetry
cbr′ = @composite brs[3] + brs[5] # (2a|A₁) + (3c|B₂)
tbm′ = tb_hamiltonian(cbr′)
ptbm′ = tbm′([2.5, 0, 0.2, 0, -1, 0])
plot(kpi, spectrum(ptbm′, kpi); annotations = collect_irrep_annotations(ptbm′))
```

Next, we turn on hybridization (controlled by the fifth term of `tbm′`):
```@example band-symmetry
ptbm′′ = tbm′([2.5, 0, 0.2, 0, -1, .5])
plot(kpi, spectrum(ptbm′′, kpi); annotations = collect_irrep_annotations(ptbm′′))
```

We can verify that neither of the band groupings in the hybridized band structure have the same content as either of the underlying band representations. In particular, the band symmetry content of the underlying band representations is:
```@example band-symmetry
SymmetryVector.([brs[3], brs[5]])
```

We can compare this against the band symmetry content of each of the hybridized bands is obtained from `collect_compatible` (using a compatibility analysis that involves *only* the high-symmetry **k**-points):
```@repl band-symmetry
ns = collect_compatible(ptbm′′)
```

In the hybridized model above, the symmetry content of each band grouping differs from any of the original band representations used to build the model: in particular, the assignment of the Γ₃ and Γ₆ irreps and the M₂ and M₄ irreps are inverted.

!!! details "Recovering an EBR decomposition"
    The new bands can still be interpreted as induced by band representations: the lowest bands correspond to ``(1a|E₁) + (1a|A₁)``, or `brs[end]+brs[end-5]`, while the upper bands are topologically fragile with a possible decomposition of the form ``(2b|A₁) + (1a|B₂) + (1a|E₂) - (3c|A₁)``, or `brs[5] + brs[end-3] + brs[end-1] - brs[1]`.
    These expansions can be obtained using [SymmetryBases.jl](https://github.com/thchr/SymmetryBases.jl)'s `decompose` function.

## Band topology

Using the extracted symmetry vectors, we can compute the associated symmetry-diagnosable band topology of each band grouping. We can obtain e.g., a coarse topological diagnosis of `TRIVIAL` (encompassing both trivial _and_ fragile phases) and `NONTRIVIAL` using Crystalline.jl's `calc_topology`:
```@repl band-symmetry
calc_topology.(ns, Ref(brs))
```

I.e., in this example, both band representations are either a trivial or a fragile phase. To resolve this distinction, we can use [SymmetryBases.jl](https://github.com/thchr/SymmetryBases.jl)'s `calc_detailed_topology`:
```@repl band-symmetry
using SymmetryBases
calc_detailed_topology.(ns, Ref(brs))
```
