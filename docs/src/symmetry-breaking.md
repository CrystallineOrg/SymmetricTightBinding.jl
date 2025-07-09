# Symmetry breaking

A frequent question in tight-binding modelling is whether -- and which -- new hoppings terms might become allowed if the overall symmetry is reduced, either by breaking spatial symmetries or time-reversal symmetry. Such terms might e.g., break degeneracies or enable topological phase transitions.

SymmetricTightBinding.jl exports `complement` as a tool to answer exactly this question. Here, we apply it to understand the effect of symmetry breaking on a 2-band model in plane group *p*4mm (⋕11).

We start by constructing our symmetry-unbroken model, picking the (2c|A₁) band representation of *p*4mm:

```@repl symmetry-break
using Crystalline, SymmetricTightBinding
brs = calc_bandreps(11, Val(2))
cbr = @composite brs[1]
tbm = tb_hamiltonian(cbr, [[0,0], [1,0]])
ptbm = tbm([0, 1, -1, 1])
```

!!! note "Interpretation of tight-binding terms"
    We can visualize the tight-binding terms using `plot`, providing also a lattice basis for the illustration:
    ```@example symmetry-break
    using GLMakie
    plot(tbm, directbasis(11, Val(2)))
    ```
    While terms 3 and 4 appear identical in this visualization, they are not (cf. [issue #75](https://github.com/CrystallineOrg/SymmetricTightBinding.jl/issues/75)): in term 3, horizontal hoppings (δ₁ and δ₂) are associated to the first Wyckoff position (at `[1/2, 0]`) and vertical hoppings to the second Wyckoff position (at `[0, 1/2]`), and vice versa for term 4.

The parameterized model has a quadratic degeneracy at M, associated with the M₅ irrep:

```@repl symmetry-break
using Brillouin, GLMakie
Rs = directbasis(11, Val(2))
kp = irrfbz_path(11, Rs)
kpi = interpolate(kp, 100);
```

```@example symmetry-break
plot(kpi, spectrum(ptbm, kpi))
```

Next, we may study which terms become allowed if we break either time-reversal symmetry or the mirror symmetries of the system. In the latter case, this corresponds to lowering the plane group symmetry to *p*4 (⋕10):
```@repl symmetry-break
Δtbm_mirror = complement(tbm, 10)                   # maintain TR, break mirror symmetry
Δtbm_tr = complement(tbm, 11; timereversal = false) # maintain spatial symmetries, break TR
```

!!! todo
    Is the above for `Δtbm_tr` wrong? It looks wrong relative to the result below. Maybe the effect of breaking either is equivalent, but the effect of breaking both is "more than the sum of their parts" - would make sense.

We can also break both simultaneously:
```@repl symmetry-break
Δtbm = complement(tbm, 10; timereversal = false)
```