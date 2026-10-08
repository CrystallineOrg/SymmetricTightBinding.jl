# Symmetry breaking

A frequent question in tight-binding modelling is whether -- and which -- new hoppings terms might become allowed if the overall symmetry is reduced, either by breaking spatial symmetries or time-reversal symmetry. Such terms might e.g., break degeneracies or enable topological phase transitions.

SymmetricTightBinding.jl exports `subduced_complement` as a tool to answer exactly this question. Here, we apply it to understand the effect of symmetry breaking on a 2-band model in plane group *p*4mm (⋕11).

We start by constructing our symmetry-unbroken model, picking the (2c|A₁) band representation of *p*4mm:

```@repl symmetry-break
using Crystalline, SymmetricTightBinding
brs = bandreps(11, Val(2))
cbr = @composite brs[1]
Rs = [[0,0], [1,0]]
tbm = tb_hamiltonian(cbr, Rs)
ptbm = tbm([0, 1, -1, 1, 0])
```

!!! note "Interpretation of tight-binding terms"
    We can visualize the tight-binding terms using `plot`, providing also a lattice basis for the illustration:
    ```@example symmetry-break
    using GLMakie
    plot(tbm[1:4], directbasis(11, Val(2)))
    ```
    We show just the first four terms (the omitted 5th term is a longer-range hopping).

The parameterized model has a quadratic degeneracy at M, associated with the M₅ irrep:

```@repl symmetry-break
using Brillouin, GLMakie
kp = irrfbz_path(11, directbasis(11, Val(2)));
kpi = interpolate(kp, 100);
```

```@example symmetry-break
plot(kpi, spectrum(ptbm, kpi))
```

We can study whether any additional terms become allowed if we reduce the symmetry. For instance, we might break the 4-fold rotational symmetry, reducing the plane group symmetry from *p*4mm (#11) to *p*2mm (#6):

```@repl symmetry-break
Δtbm_C₄ = subduced_complement(tbm, Rs, 6) # break 4-fold rotation sym.
```

This allows four additional terms.

!!! note "Why must `Rs` be passed to `subduced_complement`?"
    `tbm` does not carry information about the `Rs` that it was built with: only the resulting terms. These terms need not include hoppings for every element of `Rs` in general, since some `Rs[i]` hoppings may be forbidden in the parent group (here, *p*4mm). On symmetry reduction, such terms may nevertheless be allowed in the subgroup.
    
    `Rs` must therefore be provided and must be the same as the one `tbm` was built with: a larger `Rs` would also return longer-range terms that *p*4mm allows (a subset is also possible, resulting then in a corresponding subset of the full complement of subgroup terms). Likewise, `tbm` should be the full model `tb_hamiltonian(cbr, Rs)`, not a subset of its terms.

Conversely, we could have also tried to break time-reversal or mirror symmetry (in the latter case, reducing the plane group symmetry to *p*4 (⋕10)). Each of these allows just one new term, involving the longer-range hoppings of term 5 [^1]:

```@repl symmetry-break
Δtbm_m  = subduced_complement(tbm, Rs, 10)                       # break mirror
Δtbm_tr = subduced_complement(tbm, Rs, 11; timereversal = false) # break TR
```

[^1]: Including longer-range hoppings in `Rs` generally allows further symmetry-breaking terms; e.g., for a model built with `Rs = [[0,0], [1,0], [1,1]]`, breaking mirror symmetry allows two new terms.

We can also break both mirror and time-reversal symmetries simultaneously:

```@repl symmetry-break
Δtbm_mtr = subduced_complement(tbm, Rs, 10; timereversal = false) # break mirror & TR
```

That the result is effectively "more than the sum of the parts" of breaking time-reversal and mirror symmetry individually is merely a reflection of the fact that the additional terms are only allowed when _both_ time-reversal symmetry and mirror symmetry are broken (or, put differently, these terms are not invariant under either symmetry, and so forbidden in the presence of either).

!!! note "Subgroup relationships"
    To determine which symmetry-reductions are possible -- or, equivalently, which subgroups a particular group might have -- use Crystalline.jl's `maximal_subgroups(num(tbm))`.
    Note, however, that `subduced_complement` only allows subgroup-relationships that do not involve a change of unit cell volume; i.e., the subgroup relationship cannot be associated with a change of translational symmetry.

We can build a new "total" model, incorporating both the original terms as well as any additional symmetry-breaking terms by using `vcat`. For instance, we could incorporate the mirror-and-time-reversal symmetry breaking term into the original model:

```@repl symmetry-break
tbm′ = vcat(tbm, Δtbm_mtr)
```

And we can then verify that the original band degeneracy at M is split when the mirror-and-time-reversal-breaking term is nonzero:
```@repl symmetry-break
ptbm′ = tbm′([0, 1, -1, 1, 0 #= original terms =#,
              0.1, 0, 0, 0   #= symmetry breaking =#])
```

```@example symmetry-break
plot(kpi, spectrum(ptbm′, kpi))
```

!!! warning "Symmetry analysis in a symmetry-broken setting"
    While the "total" models `tbm′` and `ptbm′` work well for e.g., band structure purposes, they are _not_ amenable to symmetry analysis in the symmetry-reduced setting. This is because the associated band representations in `tbm′`, which are used to infer the "ingredients" of the symmetry analysis, still refer to the original group's symmetry, i.e., to *p*11 rather than *p*10.

    This may change in future versions of SymmetricTightBinding.jl, depending on available time.

