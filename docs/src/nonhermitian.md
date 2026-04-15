# [Non-Hermitian tight-binding models](@id nonhermitian)

By default, the models returned by `tb_hamiltonian` are Hermitian. In addition to Hermitian models, however, it is also possible to return anti-Hermitian or generically non-Hermitian (neither Hermitian or anti-Hermitian) models. Here, we showcase the latter.

## Hatano--Nelson model

The archetypical non-Hermitian model is the 1D Hatano--Nelson model, consisting of a single site and nearest-neighbor hoppings, in a setting with only time-reversal symmetry.
It is simple to build this model with SymmetricTightBinding.jl:

```@example hatano-nelson
using Crystalline, SymmetricTightBinding
brs = calc_bandreps(1, 1) # EBRs in plane group 1, with time-reversal symmetry
pin_free!(brs, [1=>[0]])  # the 1a Wyckoff position in plane group 1 has a free parameter: set to 0 for definiteness
cbr = @composite brs[1]   # single-site model
tbm = tb_hamiltonian(cbr, [[1]], NONHERMITIAN) # nearest neighbor hoppings
```

The model is very simple: two different hopping terms, corresponding to right- and left-directed hops. 
It is the absence of hermiticity that allows the hopping amplitudes to differ in the two directions, in clear contrast to the Hermitian case:

```@example hatano-nelson
tb_hamiltonian(cbr, [[1]], HERMITIAN)
```

Of course, the non-Hermitian model reduces to the Hermitian model when the left- and right-directed hopping amplitudes are equal.
However, when the two amplitudes are unequal, the Hatano--Nelson model features nontrivial spectral winding in the sense that its spectrum traces out a finite-area spectral loop in the complex plane as its momentum is varied across one loop.
We can see this by visualizing the complex energy as we vary $k$ from -1/2 to 1/2:

```@example hatano-nelson
ptbm = tbm([0.8, 1.2]) # a model with 0.8 hopping amplitude to right, 1.2 to the left

using GLMakie
update_theme!(linewidth = 4)

ks = range(-1/2, 1/2, 500) # 500 sampling points
Es = spectrum(ptbm, ks) # 500×1 matrix
Es = vec(Es)            # convert to vector
lines(real(Es), imag(Es); axis = (; autolimitaspect = 1))
```

The loop is associated with a quantized spectral winding $\nu = (2\pi \mathrm{i})^{-1}\oint \mathrm{d}k\, \partial_k \log E(k) = \pm 1$ when the two hopping amplitudes are unequal.


### Breaking time-reversal symmetry

We can also create models that do not assume time-reversal symmetry. In our context, this allows additional hopping terms, differing only from the time-reversal symmetric Hatano--Nelson terms by having overall imaginary prefactor:

```@example hatano-nelson
brs_notr = calc_bandreps(1, 1; timereversal=false)  # EBRs in plane group 1, without time-reversal symmetry
pin_free!(brs_notr, [1=>[0]])
tbm_notr = tb_hamiltonian((@composite brs_notr[1]), [[0], [1]], NONHERMITIAN) # on-site terms _and_ nearest-neighbor hoppings
```


## Two-band Hatano--Nelson-like model with exceptional points

We can generalize the single-band Hatano--Nelson model from above simply by constructing a two-band model, with two orbitals placed at the unit cell center.
To do so, we simply include two copies of the same EBR in the composite band representation provided to `tb_hamiltonian`: each copy is treated as a separate physical orbital.

```@example hatano-nelson
cbr = @composite 2brs[1]
tbm = tb_hamiltonian(cbr, [[0], [1]], Val(NONHERMITIAN))
summary(tbm)
```

The resulting model has 10 free terms: the first 6 are self-couplings (self-energies and intercell hoppings between the same physical orbital); the last 4 are hoppings between the two distinct orbitals. We restrict the model to this subset in the interest of simplicity:

```@example hatano-nelson
tbm = tbm[7:10]
```

The four terms span a model of the kind:

```math
\mathbf{H}(k) = 
\begin{bmatrix}
    0 & t_1 + t_2 \mathrm{e}^{2\pi\mathrm{i}k} \\
    t_1' + t_2' \mathrm{e}^{-2\pi\mathrm{i}k} & 0
\end{bmatrix},
```

with intercell hoppings $t_2^{(\prime)}$ and intracell hoppings $t_1^{(\prime)}$.
A specific instance of this model can be created from `tbm` via `tbm([t₁, t₂, t₁′, t₂′])`.

An simple special case -- but interesting, as we will see -- is equal inter- and intracell hoppings from orbital 2 to orbital 1 ($2\rightarrow 1$ hopping), i.e., $t_1 = t_2 = 1$, fully suppressed intercell $1 \rightarrow 2$ hopping, i.e., $t_2' = 0$ and free intracell $1 \rightarrow 2$ hopping, i.e., $t_1' = t$.
With this restriction, the model features an exceptional point (band degeneracy without a complete associated eigenfunction basis), occuring when $t_1 + t_2 \mathrm{e}^{2\pi\mathrm{i}k} = 1 + \mathrm{e}^{2\pi\mathrm{i}k} = 0$, i.e., when $\mathrm{e}^{2\pi\mathrm{i}k} = -1 \Leftrightarrow k = \pm 1/2$ (the BZ edge).
At the exceptional point, the Bloch Hamiltonian is defective in the sense that it is similar to a Jordan block $\big[\begin{smallmatrix} 0 & 1 \\ 0 & 0\end{smallmatrix}\big]$:

```@example hatano-nelson
t = .5 # intracell 1→2 hopping amplitude t₁′
t₁, t₂, t₁′, t₂′ = 1, 1, t, 0
ptbm = tbm([t₁, t₂, t₁′, t₂′])
ptbm([1/2]) # evaluate the Bloch Hamiltonian at k = 1/2
```

We can visualize the resulting spectrum of the model over the Brillouin zone to learn more:

```@example hatano-nelson
ks = range(-1/2, 1/2, 500)
Es = spectrum(ptbm, ks)
Es_re = real(Es) # real parts
Es_im = imag(Es) # imaginary parts

faxp = lines(ks, Es_re[:,1]; color=:royalblue, label = "Re " * rich("E", font=:italic) * subscript("1"))
lines!(ks, Es_im[:,1]; color=:royalblue,  label = "Im " * rich("E", font=:italic) * subscript("1"), linestyle=:dash)
lines!(ks, Es_re[:,2]; color=:firebrick2, label = "Re " * rich("E", font=:italic) * subscript("2"))
lines!(ks, Es_im[:,2]; color=:firebrick2, label = "Im " * rich("E", font=:italic) * subscript("2"), linestyle=:dash)
faxp.axis.xlabel = "Momentum " * rich("k", font=:italic)
faxp.axis.ylabel = "Energy " * rich("E", font=:italic)
axislegend(faxp.axis; framevisible=false)
xlims!(-1/2, 1/2)
faxp # hide
```

The spectrum is degenerate at $k = \pm 1/2$ as expected, generally complex, and exhibiting both time-reversal symmetry $E_1(k) = E_2(-k)^*$ and ``accidental'' particle-hole symmetry $E_1(k) = -E_1(k)$ (resulting from our restriction to a small set of hopping hoppings terms).


### Exceptional points with PT symmetry

Exceptional points are especially interesting in contexts where the Hamiltonian is not only non-Hermitian but also PT-symmetric (inversion and time).
The previous Hatano--Nelson-like model, however, is T-symmetric (by default, `calc_bandreps` assumes time-reversal symmetry, and this assumption is propagated via `brs` and `cbr` to `tb_hamiltonian`) but inversion-broken, and so lacks PT-symmetry.

We can build a variant, however, that breaks both P and T but retains PT symmetry.
To do so, first construct the terms of a time-reversal model, starting now with a set of time-reversal broken EBRs:

```@example PT-symmetry
using Crystalline, SymmetricTightBinding # hide
brs = calc_bandreps(1, 1; timereversal=false) # a single EBR, as before, but now without assumption of time-reversal
pin_free!(brs, [1=>[0]]) # as before, pin free parameters of the EBR's Wyckoff position
cbr = @composite 2brs[1]
tbm = tb_hamiltonian(cbr, [[0], [1]], Val(NONHERMITIAN))
summary(tbm)
```

The resulting model has no less than 20 possible terms: simply due to being a fully unconstrained problem -- lack both hermitian, spatial, and time-reversal symmetry.
We pick a small subset of these terms, with the aim of building a simple model. In particular, we retain imaginary onsite terms and the terms in our previous reduced model:

```@example PT-symmetry
onsite_terms  = [2, 8]
hopping_terms = [13, 15, 17, 19]
tbm = tbm[vcat(onsite_terms, hopping_terms)]
```

The span of these terms result in a Bloch Hamiltonian:

```math
\mathbf{H}(k) = 
\begin{bmatrix}
    \mathrm{i}\gamma_1 & t_1 + t_2 \mathrm{e}^{2\pi\mathrm{i}k} \\
    t_1' + t_2' \mathrm{e}^{-2\pi\mathrm{i}k} & \mathrm{i}\gamma_2
\end{bmatrix}.
```

Choosing $\gamma_1 = -\gamma_2 = \gamma$, $t_1 = t_1'$, and $t_2 = t_2'$ we obtain a PT symmetric model:

```math
\mathbf{H}(k) = 
\begin{bmatrix}
    \mathrm{i}\gamma & t_1 + t_2 \mathrm{e}^{2\pi\mathrm{i}k} \\
    t_1 + t_2 \mathrm{e}^{-2\pi\mathrm{i}k} & -\mathrm{i}\gamma
\end{bmatrix}.
```

The spectrum of the model is $E_\pm(k) = \sqrt{|t_1+t_2\mathrm{e}^{2\pi\mathrm{i}k}|^2 - \gamma^2}$, which is degenerate -- in fact, exceptional -- when $|t_1+t_2\mathrm{e}^{2\pi\mathrm{i}k}| = \gamma$ (assuming $\gamma>0$). The spectrum is qualitatively distinct before and after the exceptional point:
- When $|t_1+t_2\mathrm{e}^{2\pi\mathrm{i}k}| > \gamma$: $E_\pm(k)$ is real (PT-unbroken phase).
- When $|t_1+t_2\mathrm{e}^{2\pi\mathrm{i}k}| < \gamma$: $E_\pm(k)$ is imaginary (spontaneously broken PT-symmetry).

We can see this readily by constructing the associated Hamiltonian from `tbm`, noting that `tbm([γ₁, γ₂, t₁, t₂, t₁′, t₂′])` corresponds to the general model and `tbm([γ, -γ, t₁, t₂, t₁, t₂])` to the PT-symmetric one:

```@example PT-symmetry
γ, t₁, t₂ = 0.5, 1, 1
ptbm = tbm([γ, -γ, t₁, t₂, t₁, t₂])

# calculate spectrum
ks = range(-1/2, 1/2, 500)
Es = spectrum(ptbm, ks)
Es_re = real(Es) # real parts
Es_im = imag(Es) # imaginary parts
Es_re = sort(Es_re; dims=2) # necessary to explicitly sort for visualization, due to intrinsic
Es_im = sort(Es_im; dims=2) # difficult of sorting floating-point rounded complex numbers

# plot spectrum
using GLMakie # hide
update_theme!(linewidth = 4) # hide
f = Figure()
ax = Axis(f[1,1])
lines!(ks, Es_re[:,1]; color=:royalblue,  label="Re " * rich("E", font=:italic) * subscript("−"))
lines!(ks, Es_im[:,1]; color=:royalblue,  label="Im " * rich("E", font=:italic) * subscript("−"), linestyle=:dash)
lines!(ks, Es_re[:,2]; color=:firebrick2, label="Re " * rich("E", font=:italic) * subscript("+"))
lines!(ks, Es_im[:,2]; color=:firebrick2, label="Im " * rich("E", font=:italic) * subscript("+"), linestyle=:dash)
ax.xlabel = "Momentum " * rich("k", font=:italic)
ax.ylabel = "Energy " * rich("E", font=:italic)
axislegend(ax; framevisible=false)
xlims!(-1/2, 1/2)
f # hide
```

## Non-Hermitian SSH model
***WIP***

```@example nonhermitian-ssh
using Crystalline, SymmetricTightBinding # hide
# (1b|A′) ⊕ (1a|A′) in 1D SG 2 (inversion symmetry); with intra-cell hoppings & onsite terms
brs = calc_bandreps(2, Val(1))
cbr = @composite brs[1] + brs[3]
tbm = tb_hamiltonian(cbr, [[0], [2]], Val(NONHERMITIAN))

# retain only inter-orbital (offdiagonal) terms for simplicity
tbm = tbm[5:8]
ptbm = tbm([.5, 1, 1, .5]) # antisymmetric hopping pattern

# calculate spectrum
ks = range(-1/2, 1/2, 500)
Es = spectrum(ptbm, ks)
Es_re = real(Es) # real parts
Es_im = imag(Es) # imaginary parts
Es_re = sort(Es_re; dims=2) # necessary to explicitly sort for visualization, due to intrinsic
Es_im = sort(Es_im; dims=2) # difficult of sorting floating-point rounded complex numbers

# plot spectrum
using GLMakie # hide
update_theme!(linewidth = 4) # hide
f = Figure()
ax = Axis(f[1,1])
lines!(ks, Es_re[:,1]; color=:royalblue,  label="Re " * rich("E", font=:italic) * subscript("−"))
lines!(ks, Es_im[:,1]; color=:royalblue,  label="Im " * rich("E", font=:italic) * subscript("−"), linestyle=:dash)
lines!(ks, Es_re[:,2]; color=:firebrick2, label="Re " * rich("E", font=:italic) * subscript("+"))
lines!(ks, Es_im[:,2]; color=:firebrick2, label="Im " * rich("E", font=:italic) * subscript("+"), linestyle=:dash)
ax.xlabel = "Momentum " * rich("k", font=:italic)
ax.ylabel = "Energy " * rich("E", font=:italic)
axislegend(ax; framevisible=false)
xlims!(-1/2, 1/2)
f # hide
```

## A more complicated example: exceptional lines in p4

We can also construct more complicated examples where symmetry plays a role.
Consider for example a non-Hermitian model on a 2D lattice with p4 symmetry, obtained by placing *s*-like orbitals at the two symmetry-related edges of the unit cell (i.e., a (2c|A) orbital):

```@example nonhermitian-p4
using Crystalline, SymmetricTightBinding # hide
brs = calc_bandreps(10, Val(2))
cbr = @composite brs[1] # pick the (2c|A) EBR
tbm_H  = tb_hamiltonian(cbr, [[0,0], [1,0]], Val(HERMITIAN))
tbm_NH = tb_hamiltonian(cbr, [[0,0], [1,0]], Val(NONHERMITIAN))
```

It is instructive to visualize both the Hermitian and non-Hermitian models and compare the involved hopping terms:

```@example nonhermitian-p4
using GLMakie # hide
update_theme!(linewidth = 4) # hide
plot(tbm_H)
```

```@example nonhermitian-p4
plot(tbm_NH)
```

The second and third terms of the non-Hermitian model clearly break Hermiticity (unless of equal amplitude).
We can verify that the model hosts exceptional lines in this case:

```@example nonhermitian-p4
tbm = tbm_NH[2:3] # retain only terms 2 and 3
ptbm_NH = tbm([1.2, 0.8]) # 1.2 vs. 0.8 hopping asymmetry
ks = range(-1/2, 1/2, 201)
Em = [spectrum_single_k(ptbm_NH, [k1, k2]) for k1 in ks, k2 in ks]

# sort the real and imaginary parts of the energies for plotting purposes
Em_re = [sort(real(_Es)) for _Es in Em]
Em_im = [sort(imag(_Es)) for _Es in Em]
Es_re_1 = getindex.(Em_re, 1) # band 1
Es_im_1 = getindex.(Em_im, 1)
Es_re_2 = getindex.(Em_re, 2) # band 2
Es_im_2 = getindex.(Em_im, 2)

# plotting configurations
E_label = rich("E", font=:italic) * subscript("±")
axis_args = (;
    aspect = (1,1,.75),
    xautolimitmargin = (0,0), yautolimitmargin = (0,0),
    xlabel = rich("k"; font=:italic) * subscript("1"), xlabeloffset = 10,
    ylabel = rich("k"; font=:italic) * subscript("2"), ylabeloffset = 10,
    zlabeloffset = 35,
    xticks = [-1/2, 0, 1/2], xticklabelsvisible = false,
    yticks = [-1/2, 0, 1/2], yticklabelsvisible = false
)
colorbar_args = (; vertical = false, height = 8, ticklabelsize = 12)
shininess = 1.0
lims_re = extrema(Iterators.flatten(Em_re))
lims_im = extrema(Iterators.flatten(Em_im))

# plot the real and imaginary energy surfaces across the BZ
f = Figure(size=(600,280), figure_padding=(20,0,0,5))
rowgap!(f.layout, 0)

ax1 = Axis3(f[2,1]; zlabel="Re "*E_label, axis_args...)
p1 = surface!(ks, ks, Es_re_1; colormap=:PiYG_8, colorrange=lims_re, shininess)
surface!(     ks, ks, Es_re_2; colormap=:PiYG_8, colorrange=lims_re, shininess)
zlims!(ax1, lims_re)
Colorbar(f[1,1], p1; ticks=([.9999*lims_re[1], 0, .9999*lims_re[2]], ["min", "0", "max"]), colorbar_args...)

ax2 = Axis3(f[2,2]; zlabel="Im "*E_label, axis_args...)
p2 = surface!(ks, ks, Es_im_1; colormap=:RdBu_8, colorrange=lims_im, shininess)
surface!(     ks, ks, Es_im_2; colormap=:RdBu_8, colorrange=lims_im, shininess)
zlims!(ax2, lims_re)
Colorbar(f[1,2], p2; ticks=([.9999*lims_im[1], 0, .9999*lims_im[2]], ["min", "0", "max"]), colorbar_args...)
f # hide
```