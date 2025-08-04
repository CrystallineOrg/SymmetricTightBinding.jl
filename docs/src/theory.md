# Theory Notes

This package heavily relays on [representation theory of groups](https://en.wikipedia.org/wiki/Representation_theory_of_finite_groups) and in [band theory](https://en.wikipedia.org/wiki/Electronic_band_structure) of crystals. Almost all this theory was introduced before, and, can be found in [Bradley & Cracknell](https://academic.oup.com/book/54787) and, later, developed by [Bradlyn *et al.*](https://www.nature.com/articles/nature23268). Here, we aim to make a practical introduction to the main concepts and deduce the essential functions and relations that we needed for the implementation of this package. Additionally, we generalize some of the results previously derived and make them more accessible to the general public.

## Table of contents

- [Theory Notes](#theory-notes)
  - [Table of contents](#table-of-contents)
  - [Introduction](#introduction)
  - [Transformation properties of orbitals](#transformation-properties-of-orbitals)
    - [Transformation properties of induced Bloch functions](#transformation-properties-of-induced-bloch-functions)
  - [Build a tight-binding Hamiltonian from a set of symmetric orbitals](#build-a-tight-binding-hamiltonian-from-a-set-of-symmetric-orbitals)
    - [Transformation properties of the Bloch states](#transformation-properties-of-the-bloch-states)
      - [Transformation properties under lattice translations](#transformation-properties-under-lattice-translations)
      - [Transformation properties under symmetry operations](#transformation-properties-under-symmetry-operations)
  - [Appendix A](#appendix-a)
    - [Transformation properties within Convention 2](#transformation-properties-within-convention-2)
    - [Bloch Hamiltonian under Convention 2](#bloch-hamiltonian-under-convention-2)
      - [Bloch states under Convention 2](#bloch-states-under-convention-2)
        - [Transformation under lattice translations](#transformation-under-lattice-translations)
        - [Transformation properties under symmetry operations](#transformation-properties-under-symmetry-operations-1)
    - [Conversions between Convention 1 and 2](#conversions-between-convention-1-and-2)


## Introduction

The introduction of [Topological Quantum Chemistry](https://academic.oup.com/book/54787) (TQC) made a link between trivial insulators and atomic limits. It states that if a set of isolated bands can be described by a set of isolated — atomic-like — orbitals, the set must be topologically trivial. This link is determined by, first, analyzing all band symmetries of this "atomic-like" orbitals. Then, the band's set under study will be non-trivial if it doesn't fit in that list.

The analysis of the band symmetries of the isolated orbitals can be performed by placing localized, symmetric orbitals at some high-symmetry points $𝐪_α$ — [Wyckoff position](https://en.wikipedia.org/wiki/Wyckoff_positions) — with some internal symmetry — corresponding to a particular site-symmetry irrep $ρ$. Those orbitals can be labeled as $ϕ_{αi}(𝐫)$, where $i$ runs in the dimension of the irrep $ρ$, or just by $(𝐪_α|ρ)$. By applying the Fourier transform, the induced Bloch functions can be obtain as:

```math
φ_{αi,𝐤}(𝐫) = \sum_𝐭 e^{i𝐤·𝐭} ϕ_{αi}(𝐫-𝐭)
```

How $φ_{αi𝐤}(𝐫)$ transform under symmetries will define a [band representation](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.23.2824).

The idea behind this package is to use these orbitals to build a tight-binding model that respect the symmetries and topology of the (trivial) set of bands whose band representation is matched by such set of orbitals.

To do so, it is necessary to access all possible types of symmetry-independent orbitals and high-symmetry points in every space group. Luckily, this was tabulated by [Bradley & Cracknell](https://academic.oup.com/book/54787) and implemented in the Julia package [Crystalline.jl](https://github.com/thchr/Crystalline.jl). This package is going to depend on the former package to access that information.

Now that the basic framework has been stablish, in the following sections, we are going to deduce how those functions will transform and how they can be used to build a symmetric tight-binding model.

## Transformation properties of orbitals

Let us assume we have identified a set of orbitals that describe the band structure of a trivial set of bands. This can be achieved through a TQC analysis of the bands, and several tools exist to perform such decomposition. For example, the package [PhotonicTightBinding.jl](https://github.com/AntonioMoralesPerez/PhotonicTightBinding.jl) provides this functionality for photonic bands.

> [!NOTE]
> For instance, Graphene's two *p*<sub>*z*</sub> orbitals sit at the 2b Wyckoff position. Although these orbitals are odd (i.e., changing sign) under mirror in the out-of-plane direction, they are even (i.e., invariant) under all in-plane symmetries, including rotations and mirrors. The corresponding site-symmetry irrep is the A₁ irrep of the 2b Wyckoff position. Thus, these orbitals can be expressed as (2b|A₁). However, take into account that not all band representations can be induced from a set of atomic (electronic) orbitals. Some might correspond to a hybridization or a complex mixture of them.

Let us denote by $α$ as the site in the Wyckoff position where the orbital is located, and $i$ labels the number of orbitals associated to that site. Then, the orbital $i$ at site $α$ can be denoted as $ϕ_I(𝐫)$, where we introduce the compound index $I=(α, i)$. The complete set of orbitals that will be needed to describe the system is obtained by considering all orbitals at all sites and all lattice translations of them, i.e., $\{ϕ_I(𝐫-𝐭)\}_{I𝐭}$, where $𝐭$ is a lattice translation vector.

We are going to focus on a particular site $𝐪_1$, whose orbitals $ϕ_{1i}(𝐫)$ will transform under a particular site-symmetry representation $ρ$ of the site-symmetry group $G_{𝐪_1}$. Then, for $h ∈ G_{𝐪_1}$, this function will transform as:

```math
h ϕ_{1i}(𝐫) = [ρ(h)]_{ji} ϕ_{1j}(𝐫)
```

Since the orbitals are localized at a Wyckoff position, there exist a coset decomposition of the space group $G$ that relates each site in the Wyckoff position, i.e., $𝐪_α = g_α 𝐪_1$ with $g_α ∈ G$.

> [!NOTE]
> The set of $\{g_α\}$, in combination with translations $T$, will generate a decomposition of $G$ with respect to $G_𝐪$:
> ```math
> G = \bigcup_α g_α (G_{𝐪_1} \ltimes T)
> ```

Thus, each function in the unit cell can be built from the ones at site $𝐪_1$ as follows:

```math
ϕ_{αi}(𝐫) = g_α ϕ_{1i}(𝐫) = ϕ_{1i}(g_α^{-1} 𝐫)
```

By extension, translated counterparts can be defined by:

```math
\{E|𝐭\} ϕ_I(𝐫) = ϕ_I(𝐫-𝐭)
```

The aforementioned coset decomposition also have an interesting implication: for any operation $g = \{R|𝐯\} ∈ G$, there is an unique choice of $β$ for each $α$ such that $g g_α = \{E|𝐭_{βα}\} g_β h$, for some $h ∈ G_{𝐪_1}$ and $𝐭_{βα} = g 𝐪_α - 𝐪_β$. The formal proof of this statement is out of the scope of this notes and can be found in this [article](https://www.nature.com/articles/nature23268). An intuitive picture of this statement is represented by the following figure:

![Coset decomposition](./figures/coset_decomposition.png)

Taking into consideration the definitions of the transformed orbitals and the previous decomposition, we deduce that the orbitals transform under the induced representation $ρ_G$ according to:

```math
ρ_G(g) ϕ_{αi}(𝐫-𝐭) = g \{E|𝐭\} ϕ_{αi}(𝐫) \\
= \{E|R𝐭\} g ϕ_{αi}(𝐫) \\
= \{E|R𝐭\} \{E|𝐭_{βα}\} g_β h g_α^{-1} ϕ_{αi}(𝐫) \\
= \{E|R𝐭 + 𝐭_{βα}\} g_β h ϕ_{1i}(𝐫) \\
= \sum_j \{E|R𝐭 + 𝐭_{βα}\} g_β [ρ(h)]_{ji} ϕ_{1j}(𝐫) \\
= \sum_j [ρ(h)]_{ji} \{E|R𝐭 + 𝐭_{βα}\} ϕ_{βj}(𝐫) \\
= \sum_j [ρ(h)]_{ji} ϕ_{βj}(𝐫 - R𝐭 - 𝐭_{βα})
```

In principle, we could use the complete set of orbitals — $\{ϕ_I(𝐫-𝐭)\}$, with all degrees of freedom $I$ and all lattice translations $𝐭$ — to build a tight-binding model. However, it is more practical (and usual) to use the translational invariance of this orbitals to define a Fourier transform, and use their Fourier transformed functions as a basis — we are going to label such functions as induced Bloch functions. By doing so, instead of working with $\dim(I) \times N$ orbitals, where $\dim(I)$ is the number of sites plus the number of orbitals at each site and $N$ is the number of unit cells; you can consider $\dim(I)$ functions evaluated at $N$ points inside the Brillouin zone.

However, when defining a Fourier transform, there is a gauge freedom which leads to different, so-called, "conventions". This choice has important implications on the representations of the symmetry operations and, even, in the representation of the Hamiltonian. Here, we are going to focus on one convention, and we are going to discuss changes and similarities with another convention in [Appendix A](#appendix-a).

### Transformation properties of induced Bloch functions

Using the translational invariance of the orbitals, we can formally define a Fourier transform of them. This functions will not be solution to any Schrödinger-like problem, so instead of calling them Bloch states we pick up the term of *induced* Bloch functions, or just Bloch functions.

As mentioned before, there is a gauge freedom on the choice of the Fourier transform. Here, we are going to choose the following one:

```math
φ_{I,𝐤}(𝐫) ≡ \sum_𝐭 e^{i𝐤·(𝐭+𝐪_α)} ϕ_I(𝐫-𝐭)
```

The main reason behind this choice is due to fact that, with this gauge choice, the 𝐤-space dependence of the space group transformations' representations enters as a global phase, as we will see. This is really convenient for computation purposes and that is why we picked it. However, this convention enforced that the Bloch functions are not periodic in reciprocal space:

```math
φ_{I,𝐤+𝐆} = \sum_𝐭 e^{i(𝐤+𝐆)·(𝐭+𝐪_α)} ϕ_I(𝐫-𝐭) \\
= \sum_𝐭 e^{i𝐆·(𝐭+𝐪_α)} e^{i𝐤·(𝐭+𝐪_α)} ϕ_I(𝐫-𝐭) \\
= e^{i𝐆·𝐪_α} φ_{I,𝐤}
```

This implies that if the orbital is located at a non-integer position in the unit-cell, i.e, located at positions that are integer combinations of lattice vectors, the phase factor will differ from unity and the Bloch function will gain a phase. This yields that, in general, Bloch functions are not periodic under reciprocal lattice translations within this convention. This has interesting implications in the computation of some parts of this package, such as, the representation of symmetry operations or symmetry eigenvalues.

Since this functions are derived from the orbitals, the transformations properties of this functions can be obtained. How this functions transform under symmetry operations will conform what is usually called a *band representation*. In particular, this band representation will be:

```math
g φ_{iα,𝐤}(𝐫) = \frac{1}{\sqrt{N}} \sum_𝐭 e^{i𝐤·(𝐭+𝐪_α)} g ϕ_{iα}(𝐫-𝐭) \\
= \frac{1}{\sqrt{N}} \sum_{𝐭,j} e^{i𝐤·(𝐭+𝐪_α)} [ρ(h)]_{ji} ϕ_{jβ}(𝐫-R𝐭-𝐭_{βα}) \\
= \frac{1}{\sqrt{N}} \sum_{𝐭',j} [ρ(h)]_{ji} e^{i𝐤·R^{-1}(𝐭'+𝐪_β-𝐯)} ϕ_{jβ}(𝐫-𝐭') \\
= e^{-i([R^{-1}]^T 𝐤)·𝐯} \frac{1}{\sqrt{N}} \sum_{𝐭',j} [ρ(h)]_{ji} e^{i([R^{-1}]^T 𝐤)·(𝐭'+𝐪_β)} ϕ_{jβ}(𝐫-𝐭') \\
= e^{-i([R^{-1}]^T 𝐤)·𝐯} \sum_j [ρ(h)]_{ji} φ_{jβ,[R^{-1}]^T 𝐤}(𝐫),
```
where we have defined $𝐭' = R𝐭 + 𝐭_{βα} ⇒ 𝐭 = R^{-1} (𝐭'-𝐭_{βα})$, and we have used the following property: $𝐤·(R 𝐫) = (R^T 𝐤)·𝐫$. Finally, if the define the action of a symmetry operation $g = \{R|𝐯\}$ on a reciprocal space vector 𝐤 as: $g𝐤 ≡ [R^{-1}]^T 𝐤$, we can rewrite the previous relation as:

```math
g φ_{I,𝐤}(𝐫) = e^{-i(g 𝐤)·𝐯} \sum_J [ρ(h)]_{JI} φ_{J,g 𝐤}(𝐫)
```

This relation will be crucial in the implementation of the package, since it states the band representation of the system. If a tight-binding model is built from this set of functions, its band structure will, by construction, replicate the one of the original system. In other words, the tight-binding model will inherit all symmetries of the system, forcing the same degeneracies the system has, and exhibit the same symmetry-indicated topology.

For the sake of simplicity, we are going to define a matrix $D_𝐤(g)$, whose entries will be conformed by the previous operation, i.e., $[D_𝐤(g)]_{JI} = e^{-i(g 𝐤)·𝐯} [ρ(h)]_{JI}$, where remember that: $I = (α,i)$, $J = (β,j)$ and $𝐭_{βα} = g 𝐪_α - 𝐪_β$. Then, we can rewrite the previous relation as:

```math
\boxed{g φ_{I,𝐤}(𝐫) = \sum_J [D_𝐤(g)]_{JI} φ_{J,g 𝐤}(𝐫)}
```

It is important to notice that the dependence on 𝐤 of the representation $D_𝐤$ is a global phase factor. This is really convenient for computational purposes when imposing the symmetry constraints in the Hamiltonian.

Then, our next objective is to build a tight-binding model that uses this functions as basis and replicates the band structure of the system. We perform this construction in the next section.

## Build a tight-binding Hamiltonian from a set of symmetric orbitals

Second quantization rephrases quantum mechanics in terms of fields and occupation numbers. Instead of tracking individual particles, we describe how many particles occupy each quantum state. This is ideal for many-body physics and that's why we are going to implement it here.

In order to do so, we need to introduce a creation and annihilation operators. Since we want to use the basis of orbitals previously introduced, we can define them as:

```math
\ket{ϕ_{I,𝐭}} ≡ ĉ^†_{I,𝐭} \ket{\text{vac}}
```

Then, the most general tight-binding Hamiltonian can be written using those operators as:

```math
Ĥ = \sum_{IJ,𝐭𝐭'} h_{IJ,𝐭-𝐭'} ĉ^†_{I,𝐭} ĉ_{J,𝐭'}
```

This Hamiltonian reads that the probability of "hopping" from an orbital $\ket{ϕ_{J,𝐭'}}$ to an orbital $\ket{ϕ_{I,𝐭}}$ is given by the amplitude term $h_{IJ,𝐭-𝐭'}$. Notice that we assumed that the hopping amplitude only depends on the relative distance between both orbitals. This implies that the Hamiltonian will be translational invariant, as it should be. In the following, we are going to refer to that distance as $𝐑 = 𝐭-𝐭'$. Realize that it must be a lattice translation. Using this definition we can rewrite the previous Hamiltonian as:

```math
Ĥ = \sum_{IJ,𝐑𝐭'} h_{IJ,𝐑} ĉ^†_{I,𝐑+𝐭'} ĉ_{J,𝐭'}
```

In order to be consistent with the previous choice of the Fourier transform, we obtain that the creation operator in reciprocal space must be related to $ĉ^†_{I,𝐭}$ as:

```math
\ket{φ_{I,𝐤}} = â_{I,𝐤}^† \ket{\text{vac}} \\
= \frac{1}{\sqrt{N}} \sum_𝐭 e^{i𝐤·(𝐭+𝐪_α)} \ket{ϕ_{I,𝐭}} = \frac{1}{\sqrt{N}} \sum_𝐭 e^{i𝐤·(𝐭+𝐪_α)} ĉ_{I,𝐭}^† \ket{\text{vac}} \\
⇒ â_{I,𝐤}^† = \frac{1}{\sqrt{N}} \sum_𝐭 e^{i𝐤·(𝐭+𝐪_α)} ĉ_{I,𝐭}^†
```

Notice that since $φ_{I,𝐤+𝐆}(𝐫) = e^{i𝐆·𝐪_α} φ_{I,𝐤}(𝐫)$, we also have that $â_{I,𝐤+𝐆}^† = e^{i𝐆·𝐪_α} â_{I,𝐤}^†$, consistently.

Considering this, we can rewrite the tight-binding Hamiltonian in reciprocal space as:

```math
Ĥ = \frac{1}{N} \sum_{IJ,𝐑𝐭'} h_{IJ,𝐑} \sum_{𝐤𝐤'} e^{i𝐤·(𝐑+𝐭'+𝐪_α)} e^{-i𝐤'·(𝐭'+𝐪_β)} â_{I,𝐤}^† â_{J,𝐤'} \\
= \frac{1}{N} \sum_{IJ,𝐑𝐭',𝐤𝐤'} h_{IJ,𝐑} e^{i𝐤·(𝐑+𝐪_α)} e^{i(𝐤-𝐤')·𝐭'} e^{-𝐤'·𝐪_β} â_{I,𝐤}^† â_{J,𝐤'} \\
= \sum_{IJ,𝐑,𝐤} h_{IJ,𝐑} e^{i𝐤·(𝐑+𝐪_α-𝐪_β)} â_{I,𝐤}^† â_{J,𝐤},
```
where we have used the property of exponential functions: $\sum_{𝐭'} e^{i(𝐤-𝐤')·𝐭'} = N δ_{𝐤𝐤'}$.

Finally, if we define $h_{IJ,𝐤} = \sum_𝐑 h_{IJ,𝐑} e^{i𝐤·(𝐑+𝐪_α-𝐪_β)}$, we obtain the usual expression for a tight-binding Hamiltonian in reciprocal space:

```math
Ĥ = \sum_{IJ,𝐤} h_{IJ,𝐤} â_{I,𝐤}^† â_{J,𝐤}
```

As shown, the hopping amplitude in reciprocal space is computed from a summatory of the real space hopping amplitudes for all lattice translations 𝐑. Usually, such summation is cut with some "arbitrary" (subjective) cutoff. One common approach is to just consider a certain number of nearest neighbors. Since we are interested in building a symmetry-constrained tight-binding model, as the symmetry-related terms might not be equal to the number of $n$-th nearest neighbors, we are going to consider a different approach. As we are going to develop later, our strategy will be focus on asking the user to provide a set of 𝐑-vectors where to look, at least, for hopping terms, and then, search for all symmetry-related terms starting from that initial, but potentially not complete, set.

As can be seen, the tight-binding Hamiltonian is diagonal in reciprocal space. This is due to the assumption that the Hamiltonian must be translational invariant. Then, it is natural to define what is usually called the *Bloch Hamiltonian* $Ĥ_𝐤$, which consist on the diagonal blocks in reciprocal space of the tight-binding Hamiltonian:

```math
Ĥ_𝐤 = \sum_{IJ} h_{IJ,𝐤} â_{I,𝐤}^† â_{J,𝐤}
```

Some general properties must be fulfilled independent of its representation, such as its periodicity in reciprocal space. However, as shown above, the creation and annihilation operators are not periodic under reciprocal lattice translations and we also have that:

```math
h_{IJ,𝐤+𝐆} = \sum_𝐑 h_{IJ,𝐑} e^{i(𝐤+𝐆)·(𝐑+𝐪_β-𝐪_α)} \\
= e^{i𝐆·(𝐪_β-𝐪_α)} \sum_𝐑 h_{IJ,𝐑} \cancel{e^{i𝐆·𝐑}} e^{i𝐤·(𝐑+𝐪_β-𝐪_α)} \\
= e^{i𝐆·(𝐪_β-𝐪_α)} h_{IJ,𝐤}
```

However, all these phase factors cancel out in the Bloch Hamiltonian so it is translational invariant in reciprocal space, as can be seeing:

```math
Ĥ_{𝐤+𝐆} = \sum_{IJ} h_{IJ,𝐤+𝐆} â_{I,𝐤+𝐆}^† â_{J,𝐤+𝐆} \\
= \sum_{IJ} e^{i𝐆·(𝐪_β-𝐪_α)} h_{IJ,𝐤} e^{i𝐆·𝐪_α} â_{I,𝐤}^† e^{-i𝐆·𝐪_β} â_{J,𝐤} \\
= \sum_{IJ} h_{IJ,𝐤} â_{I,𝐤}^† â_{J,𝐤} = Ĥ_𝐤
```

Then, we obtain the important translational invariance in reciprocal space of the Bloch Hamiltonian: $Ĥ_{𝐤+𝐆} = Ĥ_𝐤$. This property allow us to just consider the first Brillouin zone when we examine the Bloch Hamiltonian.

The Bloch Hamiltonian can be expressed as a matrix by:

```math
Ĥ_𝐤 = Â_𝐤^† H_𝐤 Â_𝐤,
```
where $Â_𝐤^† = [ â_{1,𝐤}^†, â_{2,𝐤}^†, … ]$ is a row vector collecting all creation operators, similarly with $Â_𝐤$, and $H_𝐤$ is a complex matrix which each entry is defined by: $[H_𝐤]_{IJ} ≡ h_{IJ,𝐤}$. The matrix $H_𝐤$ is the one we are going to use in our package to compute the eigenvectors and eigenvalues for each 𝐤-point.

Notice that this matrix $H_𝐤$ is strongly dependent on the Fourier transformation picked. As proved above, within this convention, this matrix is not invariant under reciprocal lattice translations. However, this does not hold under other conventions as exposed in [Appendix A](#appendix-a). This property is not suitable for computing some fundamental properties such as the symmetry eigenvalues, but it will have some computational advantages when encoding the matrix representation $H_𝐤$ in the package.

Before proceeding, we are going to deduce the constraints that the symmetries of the system impose on the matrix $H_𝐤$. This will ensure that the model replicates the symmetry and (symmetry-indicated) topology of the system. For that purpose, first, we are going to deduce how the creation and annihilation operators transform under the symmetry operations. Let us start with the creation operator:

```math
ĝ â_{I,𝐤}^† ĝ^{-1} \ket{\text{vac}} = ĝ â_{I,𝐤}^† \ket{\text{vac}} = ĝ \ket{φ_{I,𝐤}} \\
= \sum_J [D_𝐤(g)]_{JI} \ket{φ_{J,g𝐤}} = \sum_J [D_𝐤(g)]_{JI} â_{J,g𝐤}^† \ket{\text{vac}} \\
⇒ ĝ â_{I,𝐤}^† ĝ^{-1} = \sum_J [D_𝐤(g)]_{JI} â_{J,g𝐤}^†
```

Since the symmetry operations $ĝ$ are unitary, i.e., $ĝ^{-1} = ĝ^†$, we can easily deduce the transformation properties of the annihilation operator from the creation one, and it reads as:

```math
ĝ â_{I,𝐤} ĝ^{-1} = \sum_J [D_𝐤^*(g)]_{JI} â_{J,g𝐤}
```

Considering this two transformation properties of the operators, we can deduce the set of relations that the symmetry operations will enforce in the Bloch Hamiltonian. The invariance of the Hamiltonian under symmetry operations reads as:

```math
Ĥ = ĝ Ĥ ĝ^{-1}
```

Expanding the Hamiltonian in terms of the creation and annihilation operator basis leads us to:

```math
\sum_{IJ,𝐤} â_{I,𝐤}^† h_{IJ,𝐤} â_{J,𝐤} = \sum_{IJ,𝐤} ĝ â_{I,𝐤}^† h_{IJ,𝐤} â_{J,𝐤} ĝ^{-1} \\
= \sum_{IJ,𝐤} ĝ â_{I,𝐤}^† ĝ^{-1} h_{IJ,𝐤} ĝ â_{J,𝐤} ĝ^{-1} \\
= \sum_{IJ,𝐤,I'J'} [D_𝐤(g)]_{I'I} â_{I',g𝐤}^†  h_{IJ,𝐤} [D_𝐤^*(g)]_{J'J} â_{J',g𝐤} \\
= \sum_{𝐤,I'J'} â_{I',g𝐤}^† [D_𝐤(g) H_𝐤 D_𝐤^†(g)]_{I'J'} â_{J',g𝐤} \\
⇒ \boxed{H_{g𝐤} = D_𝐤(g) H_𝐤 D_𝐤^†(g)}
```

This symmetry constraints strongly restrict the functional form of $H_𝐤$. Rather than being a completely general Hermitian (or anti-Hermitian) matrix, $H_𝐤$ must now lie in the subspace of matrices that fulfill the previous constraints. This ensures that the model preserves all symmetries and reproduces the correct degeneracies and connectivity of the original band structure.

Additionally, as exposed above, the 𝐤-dependence on the representation matrices of operations $D_𝐤$ is only a global phase factor, so it can be dropped in the previous relation. This is really practical in the implementation of the package since the 𝐤-dependence on the previous relation will be just located at the matrix $H_𝐤$ making it easier to encode in non-symbolic programming languages as Julia.

As stated previously, we are interested on diagonalizing this matrix and find the eigenvectors and eigenvalues associated to it at each 𝐤-point. Those eigenvectors will correspond to a vector of coefficients, associated to the basis set we built the Bloch Hamiltonian on, and will describe the Bloch state of the system at a particular 𝐤-point and energy. In the following section, we will elaborate on this topic and will develop how this Bloch states will transform under the symmetry operations of the system. A sanity check will be to compare the band representation of both, the real system's band structure and the tight-binding model's band structure, which should be equal.

### Transformation properties of the Bloch states

Until now we have focus on building a symmetry-constrained Hamiltonian to model the band structure of a physical system. Now, we shift our attention to analyze the properties of the Bloch states of the model, which must replicate the band structure of the physical system.

Let us start by defining the eigenvalue problem from where we start:

```math
H_𝐤 𝐰_{n,𝐤} = E_{n𝐤} 𝐰_{n,𝐤},
```
where $\{E_{n𝐤}\}$ is the set of eigenvalues (energies) at each 𝐤-point and $\{𝐰_{n,𝐤}\}$ the set of eigenvectors associated to them. Each eigenvector is a vector of coefficients which will correspond to a particular Bloch state in the basis used for describing the Bloch Hamiltonian, i.e.:

```math
\ket{ψ_{n𝐤}} = \sum_I w_{I,n𝐤} \ket{φ_{I,𝐤}} = \frac{1}{\sqrt{N}} \sum_{I,𝐭} w_{I,n𝐤} e^{i𝐤·(𝐭+𝐪_α)} \ket{ϕ_{I,𝐭}}
```

Now that the Bloch states have been defined, we can deduce their transformation properties. First, we are going to analyze their transformation properties under lattice translations in real and reciprocal space, and later, we will analyze their transformation properties under symmetry operation of the space group of the crystal.

#### Transformation properties under lattice translations

Here we are going to analyze how the Bloch states of the model transform under lattice translations in real and in reciprocal space. They should transform as Bloch functions transform under these transformations.

Firstly, let us start with lattice translations in real space:

```math
ψ_{n𝐤}(𝐫+𝐑) = \braket{𝐫+𝐑|ψ_{n𝐤}} \\
= \frac{1}{\sqrt{N}} \sum_{I,𝐭} w_{I,n𝐤} e^{i𝐤·(𝐭+𝐪_α)} \braket{𝐫+𝐑|ϕ_{n,𝐭}} \\
= \frac{1}{\sqrt{N}} \sum_{I,𝐭} w_{I,n𝐤} e^{i𝐤·(𝐭+𝐪_α)} \braket{𝐫|ϕ_{n,𝐭-𝐑}} \\
= e^{i𝐤·𝐑} \frac{1}{\sqrt{N}} \sum_{I,𝐭} w_{I,n𝐤} e^{i𝐤·(𝐭-𝐑+𝐪_α)} \braket{𝐫|ϕ_{n,𝐭-𝐑}} \\
= e^{i𝐤·𝐑} \braket{𝐫|ψ_{n𝐤}} = e^{i𝐤·𝐑} ψ_{n𝐤}(𝐫)
```

The Bloch states transform as Bloch functions under translations in real space, as expected.

Secondly, let us analyze how they transform under reciprocal lattice translations. Remind that, within this convention, the matrix representation $H_𝐤$ is not periodic under reciprocal lattice translations. This implies that $𝐰_{n,𝐤}$ will not be either, but the eigenvalues $E_{n𝐤}$ must be periodic since does are the energies associated to each Bloch state — independent of the basis chosen to represent the Hamiltonian. Let us analyze this odd behavior:

```math
H_{𝐤+𝐆} 𝐰_{n,𝐤+𝐆} = E_{n,𝐤+𝐆} 𝐰_{n,𝐤+𝐆} \\
⇒ \sum_J h_{IJ,𝐤+𝐆} w_{Jn,𝐤+𝐆} = E_{n,𝐤} w_{In,𝐤+𝐆} \\
⇒ \sum_J e^{i𝐆·(𝐪_β-𝐪_α)} h_{IJ,𝐤} w_{Jn,𝐤+𝐆} = E_{n,𝐤} w_{In,𝐤+𝐆} \\
⇒ \sum_J h_{IJ,𝐤} e^{i𝐆·𝐪_β} w_{Jn,𝐤+𝐆} = E_{n,𝐤} e^{i𝐆·𝐪_α} w_{In,𝐤+𝐆}
```

Then, this implies that the eigenvectors gain a phase factor we translated in reciprocal space such that:

```math
w_{In,𝐤+𝐆} = e^{-i𝐆·𝐪_α} w_{In,𝐤}
```

To make this easier, we ca define a diagonal matrix such that $[Θ_𝐆]_{II} = e^{-i𝐆·𝐪_α}$, then the previous expression can be rewritten as:

```math
𝐰_{n,𝐤+𝐆} = Θ_𝐆 𝐰_{n,𝐤}
```

Notice that this transformation is not a simple phase factor — which is indeterminate on eigenvectors, yet it acts differently in each entry of the eigenvector. This extra factor is crucial when analyzing the invariance of the Bloch states under reciprocal lattice translations, as we will see now.

Let us deduce how Bloch states will transform under reciprocal lattice translations:

```math
\ket{ψ_{n,𝐤+𝐆}} = \frac{1}{\sqrt{N}} \sum_{I,𝐭} w_{In,𝐤+𝐆} e^{i(𝐤+𝐆)·(𝐭+𝐪_α)} \ket{ϕ_{I,𝐭}} \\
= \frac{1}{\sqrt{N}} \sum_{I,𝐭} \cancel{e^{-i𝐆·𝐪_α}} w_{In,𝐤} \cancel{e^{i𝐆·𝐭}} \cancel{e^{i𝐆·𝐪_α}} e^{i𝐤·(𝐭+𝐪_α)} \ket{ϕ_{I,𝐭}} \\
= \frac{1}{\sqrt{N}} \sum_{I,𝐭} w_{In,𝐤} e^{i𝐤·(𝐭+𝐪_α)} \ket{ϕ_{I,𝐭}} = \ket{ψ_{n𝐤}}
```

The Bloch states will be invariant under reciprocal lattice translations. This is an important feature and must remain independently of which basis is used for representing the Hamiltonian. This property is proven for another Fourier convention in [Appendix A](#appendix-a).

#### Transformation properties under symmetry operations

Here we analyze how the Bloch states transform under more complex symmetry operations $g = \{R|𝐯\}$ that might involve translations $𝐯$ and site-symmetry operations $R$. The Bloch state will transform under this operations as:

```math
ĝ \ket{ψ_{n𝐤}} = \sum_I w_{I,n𝐤} ĝ \ket{φ_{I,𝐤}} \\
= \sum_{IJ} w_{I,n𝐤} [D_𝐤(g)]_{JI} \ket{φ_{J,g𝐤}} \\
= \sum_{IJ} [D_𝐤(g)]_{JI} w_{I,n𝐤} \ket{φ_{J,g𝐤}}
```

We are particularly interested in the transformation under operations $ĝ$ in the little-group $G_𝐤$ of a particular 𝐤-point. This operations will leave invariant the particular 𝐤-point up to a lattice translation, i.e., $g 𝐤 = 𝐤 + 𝐆$. How these functions transform under those operations at each high-symmetry point will allow us to assign an irrep to each of the Bloch states at that 𝐤-point. Those should coincide with the ones obtained from the original system's band structure. The irrep could be assigned by computing the symmetry eigenvalues associated to each Bloch state. Those are compute by:

```math
\braket{ψ_{n𝐤}|ĝ|ψ_{n𝐤}} = \sum_{IJ} (w_{I,n𝐤})^* w_{J,n𝐤} \braket{φ_{I,𝐤}|ĝ|φ_{J,𝐤}} \\
= \sum_{IJJ'} (w_{I,n𝐤})^* w_{J,n𝐤} [D_𝐤(g)]_{J'J} \braket{φ_{I,𝐤}|φ_{J',g𝐤}} \\
= \sum_{IJJ'} (w_{I,n𝐤})^* w_{J,n𝐤} [D_𝐤(g)]_{J'J} \braket{φ_{I,𝐤}|φ_{J',𝐤+𝐆}} \\
= \sum_{IJJ'} (w_{I,n𝐤})^* w_{J,n𝐤} [D_𝐤(g)]_{J'J} e^{i𝐆·𝐪_{β'}} \braket{φ_{I,𝐤}|φ_{J',𝐤}} \\
= \sum_{IJJ'} (w_{I,n𝐤})^* w_{J,n𝐤} [D_𝐤(g)]_{J'J} e^{i𝐆·𝐪_{β'}} \delta_{IJ'} \\
= \sum_{IJ} (w_{I,n𝐤})^* e^{i𝐆·𝐪_α} [D_𝐤(g)]_{IJ} w_{J,n𝐤}
```
where we have used how the Bloch functions transform under reciprocal lattice translations — a property inherit from the convention choice — and their orthogonality.

Notice that this expression has a phase factor that needs to be accounted for. In other conventions this phase factor does not appears making it easier to compute. Nevertheless, we stick to the current convention due to the property of the 𝐤-dependence in the representation matrices of symmetry operations. However, it is interesting to be able to change from one convention to others. Because of that, we include some functions in the package to be able to change from one convention to another one — heavily used in the literature. The relation between these two conventions can be found in [Appendix A](#appendix-a).

Finally, it is interesting to vectorize the previous expression in order to implemented it in the package. To do so, we make use of the previous phase factor matrix $Θ_𝐤$. Making use of it, the previous expression can be written as:

```math
\boxed{\braket{ψ_{n𝐤}|ĝ|ψ_{n𝐤}} = (Θ_𝐆 𝐰_{n𝐤}) · (D_𝐤(g) 𝐰_{n𝐤})}
```

## Appendix A

In this appendix we aim to present, develop and compare two of the main conventions present on the literature for Fourier transforms. The two Fourier transform conventions we are going to analyze are:

1. **Convention 1:** $φ^{(1)}_{I,𝐤}(𝐫) ≡ \sum_𝐭 e^{i𝐤·(𝐭+𝐪_α)} ϕ_I(𝐫-𝐭)$
2. **Convention 2:** $φ^{(2)}_{I,𝐤}(𝐫) ≡ \sum_𝐭 e^{i𝐤·𝐭} ϕ_I(𝐫-𝐭)$

where Convention 1 is the one we have been using in the theory notes and Convention 2 is another one commonly used in the literature and other packages such as [Bradlyn *et al.*](https://www.nature.com/articles/nature23268). This second convention does not includes the position of the orbital $𝐪_α$ in the phase factor of the Fourier transform.

The former is the one used in the [PythTB package](https://www.physics.rutgers.edu/pythtb/), where they claim it to be more suitable for computing topological invariants as [Berry phases](https://en.wikipedia.org/wiki/Geometric_phase) or [Wilson loops](https://en.wikipedia.org/wiki/Wilson_loop). The later is more common in the literature since it is not necessary to trace back the extra phase factor. Additionally, as we will see later, the later makes easier to compute the symmetry eigenvalues.

The arguments of which one is better than the other are out of the scope of this notes, so we are going to focus on developing both of them and pointing out their main differences. The package uses — for now — Convention 1, since it is more suitable for accounting on the 𝐤-dependence, but it also provide several tools to convert its outcome into Convention 2.

Firstly, we are going to do a similar analysis to the previous one in Convention 1, but now on Convention 2. We are going to analyze the transformation properties of the Bloch functions induced from the orbitals, the effect of this choice on the representation of the Bloch Hamiltonian and its Bloch states. Secondly, we are going to point out the main differences and similarities between both conventions. We aim to point out in which situations one more suitable than the other and when it is irrelevant. Finally, we are going to cover the conversion rules to change to one another — which are the ones we implement in this package.

### Transformation properties within Convention 2

Firstly, we are going to prove the previous statement: Convention 2 is periodic in reciprocal space, on the contrary, to Convention 1. Let us deduce how a reciprocal lattice translation $𝐆$ acts on the Bloch functions under Convention 2:

```math
φ^{(2)}_{I,𝐤+𝐆} = \sum_𝐭 e^{i(𝐤+𝐆)·𝐭} ϕ_I(𝐫-𝐭) \\
= \sum_𝐭 \cancel{e^{i𝐆·𝐭}} e^{i𝐤·𝐭} ϕ_I(𝐫-𝐭) \\
= φ^{(2)}_{I,𝐤}
```

This implies that the Bloch functions are periodic under reciprocal lattice translations within this convention. This has interesting implications in the computation of some parts of this package, such as, the representation of symmetry operations or symmetry eigenvalues. 

Secondly, let us reproduce the transformation properties of the Bloch functions under symmetry operations $g = \{ R|𝐯 \}$ of the space group under Convention 2. Following a similar approach to the one previously developed:

```math
g φ^{(2)}_{I,𝐤}(𝐫) = \frac{1}{\sqrt{N}} \sum_𝐭 e^{i𝐤·𝐭} g ϕ_{I}(𝐫-𝐭) \\
= \frac{1}{\sqrt{N}} \sum_{𝐭,J} e^{i𝐤·𝐭} [ρ(h)]_{JI} ϕ_J(𝐫-R𝐭-𝐭_{βα}) \\
= \frac{1}{\sqrt{N}} \sum_{𝐭',J} e^{i𝐤·R^{-1}(𝐭'-𝐭_{βα})} [ρ(h)]_{JI} ϕ_J(𝐫-𝐭') \\
= e^{-i(g 𝐤)·𝐭_{βα}} \sum_J [ρ(h)]_{JI} \frac{1}{\sqrt{N}} \sum_{𝐭'} e^{i(g 𝐤)·𝐭'} ϕ_J(𝐫-𝐭') \\
= e^{-i(g 𝐤)·𝐭_{βα}} \sum_J [ρ(h)]_{JI} φ^{(2)}_{J,g 𝐤}(𝐫),
```
where we made the substitution $𝐭' = R𝐭 + 𝐭_{βα}$, and used the definition stated before: $g𝐤 ≡ [R^{-1}]^T 𝐤$.

Similarly as before, we can define a representation matrix $D^{(2)}_𝐤(g)$ whose entries are $[D^{(2)}_𝐤(g)]_{IJ} = e^{-i(g𝐤)·𝐭_{βα}} [ρ(h)]_{IJ}$, where $I = (i,α)$ and $J = (j,β)$. Then, the previous expression reduces to:

```math
g φ^{(2)}_{I,𝐤}(𝐫) = \sum_j [D^{(2)}_𝐤(g)]_{JI} φ^{(2)}_{jβ,g𝐤}(𝐫)
```

Notice that the representation matrix for the space group operations differs between conventions, i.e., $D^{(1)}_𝐤(g) ≠ D^{(2)}_𝐤(g)$. The representation under Convention 1 depends on the translational part $𝐯$, as shown before, meanwhile, under Convention 2, it presents not on a global phase factor, but on a local phase factor depending on $𝐭_{βα}$.

The next step will be to build a tight-binding model using this set of functions as a basis. For that, in the following section, we will follow the same steps as in Convention 1 by introducing the creation and annihilation operators associated to such functions, and how the Bloch Hamiltonian will look like.

### Bloch Hamiltonian under Convention 2

We want to use the previously introduced Bloch functions as a basis to construct a Bloch Hamiltonian in reciprocal space. Once again, we start from the most general tight-binding Hamiltonian, which, as we described, can be written as:

```math
Ĥ = \sum_{IJ,𝐑𝐭} h_{IJ,𝐑} ĉ^†_{I,𝐭+𝐑} ĉ_{J,𝐭}
```

Since we want to use the previous Bloch functions as a basis, we have to introduce a set of creation and annihilation operators that are consistent with the convention choice — Convention 2 in this case. This is satisfied by the following relation:

```math
\hat{b}_{I,𝐤}^† = \frac{1}{\sqrt{N}} \sum_𝐭 e^{i𝐤·𝐭} ĉ_{I,𝐭}^†,
```
where we have used $\hat{b}$ as the notation for the operators under this new convention. It is interesting to notice that you can obtain one from the other by using the following relation:

```math
\hat{b}_{I,𝐤}^† = e^{-i𝐤·𝐪_α} â_{I,𝐤}^†
```

Introducing the previous transformation into the tight-binding Hamiltonian we obtain that:

```math
Ĥ = \sum_{IJ,𝐑,𝐤} h_{IJ,𝐑} e^{i𝐤·𝐑} \hat{b}_{I,𝐤}^† \hat{b}_{J,𝐤}
```

If we define $h^{(2)}_{IJ,𝐤} = \sum_𝐑 h_{IJ,𝐑} e^{i𝐤·𝐑}$, we obtain the usual expression for a tight-binding Hamiltonian in reciprocal space:

```math
Ĥ = \sum_{IJ,𝐤} h^{(2)}_{IJ,𝐤} \hat{b}_{I,𝐤}^† \hat{b}_{J,𝐤}
```

From here, we can define the Bloch Hamiltonian which will be the diagonal part of the Hamiltonian in reciprocal space, i.e.:

```math
Ĥ_𝐤 = \sum_{IJ} h^{(2)}_{IJ,𝐤} \hat{b}_{I,𝐤}^† \hat{b}_{J,𝐤}
```

Considering that now the Bloch functions are periodic it is easier to prove that the Bloch Hamiltonian is periodic also, but, nevertheless, we are going to prove it. Firstly, let us examine how the creation operator transform under a reciprocal lattice translation:

```math
\hat{b}_{I,𝐤+𝐆}^† = \frac{1}{\sqrt{N}} \sum_𝐭 e^{i(𝐤+𝐆)·𝐭} ĉ_{I,𝐭}^† \\
= \frac{1}{\sqrt{N}} \sum_𝐭 \cancel{e^{i𝐆·𝐭}} e^{i𝐤·𝐭} ĉ_{I,𝐭}^† \\
= \hat{b}_{I,𝐤}^†
```

Secondly, let us study how the coefficients transform under a reciprocal lattice translation:

```math
h^{(2)}_{IJ,𝐤+𝐆} = \sum_𝐑 h_{IJ,𝐑} e^{i(𝐤+𝐆)·𝐑} \\
= \sum_𝐑 h_{IJ,𝐑} \cancel{e^{i𝐆·𝐑}} e^{i𝐤·𝐑} \\
= h^{(2)}_{IJ,𝐤}
```

As can be seen, all the components are periodic within this convention. This is the common reason why it is heavily used in the literature. Proving now the translational invariance of the Bloch Hamiltonian:

```math
Ĥ_{𝐤+𝐆} = \sum_{IJ} h^{(2)}_{IJ,𝐤+𝐆} \hat{b}_{I,𝐤+𝐆}^† \hat{b}_{J,𝐤+𝐆} \\
= \sum_{IJ} h^{(2)}_{IJ,𝐤} \hat{b}_{I,𝐤}^† \hat{b}_{J,𝐤} = Ĥ_𝐤
```

This is an important property since the eigenvalues of this Bloch Hamiltonian, which correspond to the energies of the Bloch states of the system, must be periodic in reciprocal space, allowing us to restrict to the first Brillouin zone. Additionally, since the representation matrix $H_𝐤$ is now periodic by itself, the eigenvectors $𝐰^{(2)}_{n𝐤}$ will also be periodic. This has important implications, for example, when computing the symmetry eigenvalues.

Before studying the transformation properties of the Bloch states, we want to mention that the creation and annihilation operators and the Bloch Hamiltonian within this convention will have the same transformation properties under symmetry operations but now using the representation matrix of the operations $D^{(2)}_𝐤$ associated to Convention 2.

#### Bloch states under Convention 2

Let us now jump into the transformation properties of the Bloch states. The Bloch states are represented using the basis obtained by Convention 2 as:

```math
\ket{ψ_{n𝐤}} = \sum_I w^{(2)}_{I,n𝐤} \ket{φ^{(2)}_{I,𝐤}} = \frac{1}{\sqrt{N}} \sum_{I,𝐭} w^{(2)}_{I,n𝐤} e^{i𝐤·(𝐭+𝐪_α)} \ket{ϕ_{I,𝐭}}
```

Let us first study how the Bloch states transform under lattice translations in real and reciprocal space and, then, deduce how they transform under more complex symmetry operations.

##### Transformation under lattice translations

Firstly, we are going to consider lattice translations in real space. This can be obtained by:

```math
ψ_{n𝐤}(𝐫+𝐑) = \braket{𝐫+𝐑|ψ_{n𝐤}} \\
= \frac{1}{\sqrt{N}} \sum_{I,𝐭} w^{(2)}_{I,n𝐤} e^{i𝐤·𝐭} \braket{𝐫+𝐑|ϕ_{n,𝐭}} \\
= \frac{1}{\sqrt{N}} \sum_{I,𝐭} w^{(2)}_{I,n𝐤} e^{i𝐤·𝐭} \braket{𝐫|ϕ_{n,𝐭-𝐑}} \\
= e^{i𝐤·𝐑} \frac{1}{\sqrt{N}} \sum_{I,𝐭} w^{(2)}_{I,n𝐤} e^{i𝐤·(𝐭-𝐑)} \braket{𝐫|ϕ_{n,𝐭-𝐑}} \\
= e^{i𝐤·𝐑} \braket{𝐫|ψ_{n𝐤}} = e^{i𝐤·𝐑} ψ_{n𝐤}(𝐫)
```

The Bloch states transform as Bloch functions under translations in real space, as expected.

Secondly, let us analyze how they transform under reciprocal lattice translations. Remind that the matrix representation $H^{(2)}_𝐤$ is periodic under reciprocal lattice translations. This implies that $𝐰^{(2)}_{n,𝐤+𝐆} = 𝐰^{(2)}_{n,𝐤}$, and $E_{n,𝐤+𝐆} = E_{n,𝐤}$, i.e., they are periodic under reciprocal lattice translations.

Let us deduce how Bloch states will transform under reciprocal lattice translations:

```math
\ket{ψ_{n,𝐤+𝐆}} = \frac{1}{\sqrt{N}} \sum_{I,𝐭} w^{(2)}_{In,𝐤+𝐆} e^{i(𝐤+𝐆)·𝐭} \ket{ϕ_{I,𝐭}} \\
= \frac{1}{\sqrt{N}} \sum_{I,𝐭} w^{(2)}_{In,𝐤} \cancel{e^{i𝐆·𝐭}} e^{i𝐤·𝐭} \ket{ϕ_{I,𝐭}} \\
= \frac{1}{\sqrt{N}} \sum_{I,𝐭} w^{(2)}_{In,𝐤} e^{i𝐤·𝐭} \ket{ϕ_{I,𝐭}} = \ket{ψ_{n𝐤}}
```

The Bloch states will remain invariant under reciprocal lattice translations, as it should be.

##### Transformation properties under symmetry operations

Here we analyze how the Bloch states transform under more complex symmetry operations $g = \{R|𝐯\}$ that might involve translations $𝐯$ and site-symmetry operations $R$. The Bloch state will transform under this operations as:

```math
ĝ \ket{ψ_{n𝐤}} = \sum_I w^{(2)}_{I,n𝐤} ĝ \ket{φ^{(2)}_{I,𝐤}} \\
= \sum_{IJ} w^{(2)}_{I,n𝐤} [D^{(2)}_𝐤(g)]_{JI} \ket{φ^{(2)}_{J,g𝐤}} \\
= \sum_{IJ} [D^{(2)}_𝐤(g)]_{JI} w^{(2)}_{I,n𝐤} \ket{φ^{(2)}_{J,g𝐤}}
```

We are particularly interested in the transformation under operations $ĝ$ in the little-group $G_𝐤$ of a particular 𝐤-point. This operations will leave invariant the particular 𝐤-point up to a lattice translation, i.e., $g 𝐤 = 𝐤 + 𝐆$. How these functions transform under those operation at each high-symmetry point will allow us to assign an irrep to each of the Bloch states at that 𝐤-point. Those should coincide with the ones obtained from the original system's band structure. The irrep could be assigned by computing the symmetry eigenvalues associated to each Bloch state. Those are compute by:

```math
\braket{ψ_{n𝐤}|ĝ|ψ_{n𝐤}} = \sum_{IJ} (w^{(2)}_{I,n𝐤})^* w^{(2)}_{J,n𝐤} \braket{φ^{(2)}_{I,𝐤}|ĝ|φ^{(2)}_{J,𝐤}} \\
= \sum_{IJJ'} (w^{(2)}_{I,n𝐤})^* w^{(2)}_{J,n𝐤} [D^{(2)}_𝐤(g)]_{J'J} \braket{φ^{(2)}_{I,𝐤}|φ^{(2)}_{J',g𝐤}} \\
= \sum_{IJJ'} (w^{(2)}_{I,n𝐤})^* w^{(2)}_{J,n𝐤} [D^{(2)}_𝐤(g)]_{J'J} \braket{φ^{(2)}_{I,𝐤}|φ^{(2)}_{J',𝐤+𝐆}} \\
= \sum_{IJJ'} (w^{(2)}_{I,n𝐤})^* w^{(2)}_{J,n𝐤} [D^{(2)}_𝐤(g)]_{J'J} \braket{φ^{(2)}_{I,𝐤}|φ^{(2)}_{J',𝐤}} \\
= \sum_{IJJ'} (w^{(2)}_{I,n𝐤})^* w^{(2)}_{J,n𝐤} [D^{(2)}_𝐤(g)]_{J'J} \delta_{IJ'} \\
= \sum_{IJ} (w^{(2)}_{I,n𝐤})^* [D_𝐤(g)]_{IJ} w^{(2)}_{J,n𝐤}
```

Notice that this expression differs from the previous one due to a phase factor. This is the main reason why some authors decided to use Convention 2 instead Convention 1: it is not necessary to account for phase factors in the symmetry eigenvalues computations.

Nevertheless, it is interesting to be able to consider both conventions and that is why, in the next section, we develop conversion properties between the two conventions for several convention-dependent expressions.

### Conversions between Convention 1 and 2

Firstly, ley us start with the conversion between the Bloch functions that the different Fourier transformations inherit. Since there is just an additional phase factor, we can just convert from one convention to the other adding that extra factor as:

```math
\ket{φ^{(2)}_{I,𝐤}} = e^{-i𝐤·𝐪_α} \ket{φ^{(1)}_{I,𝐤}}
```

Obviously, the creation and annihilation operators will convert in a similar fashion, in particular, as said previously:

```math
\hat{b}_{I,𝐤}^† = e^{-i𝐤·𝐪_α} â_{I,𝐤}^†
```

Secondly, the Hamiltonian must be invariant independent of what basis we use to define it. Considering this, we can deduce how the matrix $H_𝐤$, which is representation dependent, convert from one convention to the other. Let us start with the Bloch Hamiltonian:

```math
Ĥ_𝐤 = \sum_{IJ} h^{(1)}_{IJ,𝐤} â_{I,𝐤}^† â_{J,𝐤} \\
= \sum_{IJ} e^{i𝐤·𝐪_α} h^{(1)}_{IJ,𝐤} e^{-i𝐤·𝐪_β} \hat{b}_{I,𝐤}^† \hat{b}_{J,𝐤} \\
⇒ h^{(2)}_{IJ,𝐤} = e^{i𝐤·𝐪_α} h^{(1)}_{IJ,𝐤} e^{-i𝐤·𝐪_β}
```

This allow us to convert from one representation matrix of the Hamiltonian to the other. It is more interesting to rewrite the previous relation in matrix form, which will be:

```math
H^{(2)}_𝐤 = Θ^†_𝐤 H^{(1)}_𝐤 Θ_𝐤,
```
where $Θ_𝐤$ is a diagonal matrix containing the phase factor as defined previously: $[Θ_𝐤]_{II} = e^{-i𝐤·𝐪_α}$. This relation allow us to transform from one representation of the Bloch Hamiltonian into the other.

Finally, we are interested in deducing the conversion properties of the eigenvector obtained from diagonalizing the representation matrix $H_𝐤$. Since $H^{(1)}_𝐤$ and $H^{(2)}_𝐤$ are related by a change of basis, the eigenvalue of both of them can be easily related to one each other as:

```math
𝐰^{(2)}_{n𝐤} = Θ^†_𝐤 𝐰^{(1)}_{n𝐤}
```

With this relations we are able to to go back and forward from one convention to the other, making it possible to use both depending on which one is the most suitable for each case.