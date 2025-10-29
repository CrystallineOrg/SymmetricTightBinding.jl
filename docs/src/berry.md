# [Berry curvature and Chern numbers](@id berry)

SymmetricTightBinding.jl implements some functionalities to compute the Berry curvature and associated Chern numbers of parameterized tight-binding models via [`berrycurvature`](@ref), [`chern`](@ref), and [`chern_fukui`](@ref). We illustrate the functionality here, with the canonical example of the Haldane model.

## Haldane model

The Haldane model [^1], featuring both an inversion-breaking mass term and a time-reversal breaking next-nearest-neighbor complex hopping, has *p*3 symmetry (plane group ⋕13). The model is built from two *s*-like orbitals (i.e., with site-symmetry irrep A) placed at the 1b and 1c Wyckoff positions.[^2]

[^1]: F.D.M. Haldane, *Model for a Quantum Hall Effect without Landau Levels: Condensed-Matter Realization of the "Parity Anomaly"*, [Phys. Rev. Lett. **61**, 2015 (1988)](https://doi.org/10.1103/PhysRevLett.61.2015).

[^2]: While it is natural to think of the Haldane model as associated with 6-fold symmetry and a graphene-like model, and hence with plane groups *p*6 or *p*6mm, the presence of the staggered mass term reduces the model's symmetry to *p*3. Thus, the orbitals do not associate with the 2b Wyckoff position of *p*6(mm), as is usually the case for graphene-like models, but with the two distinct Wyckoff positions that the *p*6(mm) 2b Wyckoff position subduces to when inversion symmetry is broken: this is precisely the 2b and 2c positions in *p*3.

```@example berry
using Crystalline, SymmetricTightBinding
sgnum = 13 # p3
brs = calc_bandreps(sgnum, Val(2); timereversal = false)
brs[4], brs[1] = brs[1], brs[4]; # see below
# (↑) We manually swap the positions of the (1b|A) and (1c|A) EBRs in `brs` above.
#     Without this adjustment, the model we construct below would place the 1b
#     position in the "second orbital slot" and the 1c position in the "first slot",
#     with the ordering being inherited from the sorting in `brs`. But Haldane places
#     the 1b position as the first orbital (A sites) and the 1c position as the
#     second orbital (B sites): aligning the two simply makes the comparison easier
cbr = @composite brs[1] + brs[4] # (1b|A) + (1c|A)

tbm = tb_hamiltonian(cbr, [[0,0], [0,1]])
```

The terms in `tbm` form a basis for many possible Hamiltonians, including for the Haldane model. By comparing term by term with Haldane's expressions, the correct parameterization can be determined to be:

```@example berry
haldane_model(t₁, m, t₂, ϕ) = tbm([m, t₂*cos(ϕ), t₂*sin(ϕ), -m, t₂*cos(ϕ), -t₂*sin(ϕ), t₁, 0])
```

This realizes the Haldane Hamiltonian with nearest-neighbor hopping $t_1$, a staggered mass term $m$, and a complex next-nearest-neighbor hopping $t_2 \exp(\pm\mathrm{i}\phi)$ with Haldane's zero-flux pattern. The model is gapless for $m/t_2 = 3\sqrt{3}|\sin\phi|$ and otherwise gapped when $|t_2 / t_1| < 1/3$.

The Berry curvature is nonzero at generic **k**-points and generic values of `t₁, m, t₂, ϕ`, as we can verify with the [`berrycurvature`](@ref) method:

```@example berry
ptbm = haldane_model(1.0, 0.1, 0.1, π/2)
berrycurvature(ptbm, [.2, .3], 1)
```

We can use this to e.g., visualize the Berry curvature distribution over the Brillouin zone:

```@example berry
# compute the Berry curvature over the parallelipiped Brillouin zone
Gs = dualbasis(directbasis(sgnum, Val(2)))
ks = range(-0.5, 0.5, 100) # k-points' range in reciprocal basis coordinates
k12s = (ReciprocalPoint(k1, k2) for k1 in ks, k2 in ks)
Ω = [berrycurvature(ptbm, k12, 1) for k12 in k12s] # Berry curvatures of first band

# use Brillouin.jl and GLMakie.jl to visualize the Berry curvature
using Brillouin, GLMakie
c = wignerseitz(Gs)

kxys = [cartesianize(k12, Gs) for k12 in k12s] # k-points in Cartesian coordinates
kxs = getindex.(kxys, 1)
kys = getindex.(kxys, 2)

f = Figure()
ax = Axis(f[1,1]; aspect=Makie.DataAspect())

Ωz = Ω ./ Crystalline.Bravais.volume(Gs) # transform from reduced to Cartesian setting
maxΩ = maximum(abs, Ωz)
colormap_kws = (; colormap = Reverse(:RdBu), colorrange=(-maxΩ, maxΩ))
for s1 in -1:1, s2 in -1:1 # map across adjacent parallelipiped BZs
    G = s1 * Gs[1] + s2 * Gs[2]
    kxys′ = kxys .+ Ref(G)
    surface!(ax, getindex.(kxys′, 1), getindex.(kxys′, 2), zeros(size(Ωz));
             color = Ωz, shading = NoShading, colormap_kws...)
end
plot!(c)
ax.limits = (-1.45π, 1.45π, -1.35π, 1.35π)
ax.xlabel = rich("k", subscript("x"); font=:italic)
ax.ylabel = rich("k", subscript("y"); font=:italic)

cb = Colorbar(f[1,2]; label="Berry curvature", colormap_kws...)
f
```

Of course, we can also compute the Chern number $C_n = \frac{\mathrm{i}}{2\pi} \int_{\text{BZ}} \nabla_{\mathbf{k}} \times \langle u_{n\mathbf{k}} | \nabla_{\mathbf{k}} | u_{n\mathbf{k}} \rangle \, \mathrm{d}^2 \mathbf{k}$ by integrating (i.e., summing up) the Berry curvature over the Brillouin zone. The convenience method [`chern`](@ref) does exactly this:

```@example berry
Nk = 21 # number of k-point samples per BZ direction
chern(ptbm, 1, Nk) # Chern number of first band
```

Generally, however, it will be faster and more accurate to use the [`chern_fukui`](@ref) method, which uses the approach of Fukui *et al.* [^3] to obtain a manifestly quantized Chern number:

```@example berry
chern_fukui(ptbm, 1, Nk)
```

The Fukui method is also applicable in non-Abelian settings, i.e., for degenerate band multiplets.

[^3]: Fukui, Hatsugai, & Suzuki, *Chern Numbers in Discretized Brillouin Zone: Efficient Method of Computing (Spin) Hall Conductances*, [J. Phys. Soc. Jpn. **74**, 1674 (2005)](https://doi.org/10.1143/JPSJ.74.1674).

## Phase diagrams
We can use e.g., the `chern_fukui` approach to reproduce the classical phase diagram of the Haldane model. In particular, we show below the Chern number of the Haldane model with $t_2 = t_1 / 5$ for varying $M / t_2$ and hopping phase $\phi$:

```@example berry
const t₁ = 1.0
const t₂ = t₁ / 5
haldane_model(m_div_t₂, ϕ) = haldane_model(t₁, m_div_t₂ * t₂, t₂, ϕ)

Nk = 21
N = 151
Mdivt2s = range(-3√3, 3√3, N)
ϕs = range(-π, π, N)
Cs = [SymmetricTightBinding.chern_fukui(haldane_model(Mdivt2, ϕ), 1, Nk) for ϕ in ϕs, Mdivt2 in Mdivt2s]

f, ax, p = contourf(ϕs, Mdivt2s, Cs; levels = (-1:2) .- 1/2, colormap=Reverse(:RdBu_5))
lines!(ϕs, 3√3*sin.(ϕs); color=:gray, linestyle=:dash, linewidth=4) # add analytical phase boundaries
lines!(ϕs, 3√3*sin.(.-ϕs); color=:gray, linestyle=:dash, linewidth=4)

ax.limits = (-π, π, -3√3, 3√3)
ax.xticks = ([-π, 0, π], ["-π", "0", "π"])
ax.yticks = ([-3√3, 0, 3√3], ["-3√3", "0", "3√3"])
ax.xlabel = "Hopping phase (TR breaking), ϕ"
ax.ylabel = "Mass term (inversion breaking), M/t₂"
ax.title = "Haldane model phase diagram"

cb = Colorbar(f[1,2], p)
cb.ticks = -1:1
cb.label = "Chern number"

f
```

We could also have obtained a similar-looking phase diagram by using topological quantum chemistry:

```@example berry
νs = Matrix{Int}(undef, N, N)
F = smith(stack(brs)) # Smith decomposition of `brs`; precomputed for efficiency
for (i, ϕ) in enumerate(ϕs), (j, Mdivt2) in enumerate(Mdivt2s)
    ptbm = haldane_model(Mdivt2, ϕ)
    n = first(collect_compatible(ptbm)) # get first symmetry vector
    ν = only(symmetry_indicators(n, F))
    νs[i, j] = ν
end

f, ax, p = contourf(ϕs, Mdivt2s, νs; levels = (0:3) .- 1/2, colormap=Reverse(:sunset))
lines!(ϕs, 3√3*sin.(ϕs); color=:gray, linestyle=:dash, linewidth=4) # add analytical phase boundaries
lines!(ϕs, 3√3*sin.(.-ϕs); color=:gray, linestyle=:dash, linewidth=4)

ax.limits = (-π, π, -3√3, 3√3)
ax.xticks = ([-π, 0, π], ["-π", "0", "π"])
ax.yticks = ([-3√3, 0, 3√3], ["-3√3", "0", "3√3"])
ax.xlabel = "Hopping phase (TR breaking), ϕ"
ax.ylabel = "Mass term (inversion breaking), M/t₂"
ax.title = "Haldane model phase diagram"

cb = Colorbar(f[1,2], p)
cb.ticks = 0:2
cb.label = "Symmetry indicator ν"

f
```

!!! note "Computing phase maps more efficiently"
    Phase maps can be constructed much more efficiently than by the brute-force grid-exploration adopted above (which scales with the phase map area). E.g., the [PhaseMap](https://github.com/greschd/PhaseMap) Python package uses a recursive exploration strategy that scales with the phase map boundary length. An example of interfacing with PhaseMap is provided in [`examples/phasemap.jl`](https://github.com/CrystallineOrg/SymmetricTightBinding.jl/blob/main/examples/phasemap.jl).

## Berry curvature in three-dimensional models

!!! note
    The `berrycurvature` method is also applicable to 3D models. An interface for the `chern` and `chern_fukui` methods, however, has not yet been added. However, it is still possible to do so manually, e.g., to analyze the existence of Weyl points. For instance, one may compute the Chern number associated with a given $k_3$ slice by:

    ```jl
    k3 = 0.1 # pick a k₃ plane
    Nk = 101 # number of integration points
    k12s = range(-0.5, 0.5, Nk+1)[2:end] # integration points in plane
    
    # compute flux of Ω₃ through k₃ plane (for first band)
    Φ = sum(berrycurvature(ptbm, [k1, k2, k3], 1)[3] for k1 in k12s, k2 in k12s; init=0.0) / Nk^2

    # normalize by 2π to obtain associated Chern number
    C₃ = Φ / (2π)
    ```