# Symmetry breaking

A frequent question in tight-binding modelling is whether -- and which -- new hoppings terms might become allowed if the overall symmetry is reduced, either by breaking spatial symmetries or time-reversal symmetry. Such terms might e.g., break degeneracies or enable topological phase transitions.

SymmetricTightBinding.jl exports `subduced_complement` as a tool to answer exactly this question. Here, we apply it to understand the effect of symmetry breaking on a 2-band model in plane group *p*4mm (⋕11).

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
Rs = directbasis(11, Val(2));
kp = irrfbz_path(11, Rs);
kpi = interpolate(kp, 100);
```

```@example symmetry-break
plot(kpi, spectrum(ptbm, kpi))
```

We can study whether any additional terms become allowed if we reduce the symmetry. For instance, we might breaj the 4-fold rotational symmetry, reducing the plane group symmetry from *p*4mm (#11) to *p*2mm (#6):

```@repl symmetry-break
Δtbm_C₄ = subduced_complement(tbm, 6) # break 4-fold rotation sym.
```

This allows three additional terms. Conversely, we could have also tried to break time-reversal or mirror symmetry (in the latter case, reducing the plane group symmetry to *p*4 (⋕10)). However, both of these cases allow no new terms [^1]:

```@repl symmetry-break
Δtbm_m  = subduced_complement(tbm, 10)                       # break mirror
Δtbm_tr = subduced_complement(tbm, 11; timereversal = false) # break TR
```

[^1]: For mirror symmetry-breaking, the absence of new terms is a result of looking only at a limited set of hopping orbits (in the original model `tb_hamiltonian(cbr, [[0,0], [1,0]])`): by including longer-range hopping orbits, we would eventually find new mirror-symmetry-broken terms. This is not so for time-reversal breaking, however: in *p*4mm, mirror symmetry and hermicity jointly impose an effective time-reversal symmetry.

However, by breaking both mirror and time-reversal symmetry simultaneously, additional terms do appear:

We can also break both symmetry simultaneously:

```@repl symmetry-break
Δtbm_mtr = subduced_complement(tbm, 10; timereversal = false) # break mirror & TR
```

That the result is effecitvely "more than the sum of the parts" of breaking time-reversal and mirror symmetry individually is merely a reflection of the fact that the term is only allowed when _both_ time-reversal symmetry and mirror symmetry is broken (or, put differently, the term is not invariant under either symmetry, and so forbidden in the presence of either).

!!! note "Subgroup relationships"
    To determine which symmetry-reductions are possible -- or, equivalently, which subgroups a particular group might have -- use Crystalline.jl's `maximal_subgroups(num(tbm))`.
    Note, however, that `subduced_complement` only allows subgroup-relationships that do not involve a change of unit cell volume; i.e., the subgroup relationship cannot be associated with a change of translational symmetry.

We can build a new "total" model, incorporating both the original terms as well as any additional symmetry-breaking terms by using `vcat`. For instance, we could incorporate the mirror-and-time-reversal symmetry breaking term into the original model:

```@repl symmetry-break
tbm′ = vcat(tbm, Δtbm_mtr)
```

And we can then verify that the original band degeneracy at M is split when the mirror-and-time-reversal-breaking term is nonzero:
```@repl symmetry-break
ptbm′ = tbm′([0, 1, -1, 1 #= original terms =#,
              0.1         #= symmetry breaking =#])
```

```@example symmetry-break
plot(kpi, spectrum(ptbm′, kpi))
```

!!! warning "Symmetry analysis in a symmetry-broken setting"
    While the "total" models `tbm′` and `ptbm′` work well for e.g., band structure purposes, they are _not_ amenable to symmetry analysis in the symmetry-reduced setting. This is because the associated band representations in `tbm′`, which are used to infer the "ingredients" of the symmetry analysis, still refer to the original group's symmetry, i.e., to *p*11 rather than *p*10.

    This may change in future versions of SymmetricTightBinding.jl, depending on available time.

