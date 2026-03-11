# Developer documentation

Derivations and implementation notes for the symmetry-constrained tight-binding framework.
For the polished user-facing exposition, see [`docs/src/theory.md`](../docs/src/theory.md).
!!! note 
   Substantial parts of the contents here have been superseded by [`docs/src/theory.md`](../docs/src/theory.md) and [`devdocs/trs_notes.md`](trs_notes.md). This file is slated for eventual deletion after its contents have been checked against these files, with any non-overlapping elements transferred.

## General Hamiltonian in periodic systems

Let us start with the most general non-interacting Hamiltonian:

$$
Ĥ = \sum_{IJ,𝐑𝐑'} h_{IJ,𝐑-𝐑'} \; ĉ_{I,𝐑}^† ĉ_{J,𝐑'},
$$
where $I,J$ wrap up the internal degrees of freedom of the orbitals and the sites, i.e., 
$I = (i, α)$; and $𝐑,𝐑'$ run over the lattice translations.

> [!WARNING]
> Notice that we have assumed that hopping terms only depends on relative distances.
> We are going to denote $𝘁 ≡ 𝐑' - 𝐑$.

We can apply the same Fourier transform to go into reciprocal space:

$$
ĉ_{I,𝐑} = \frac{1}{\sqrt{N}} \sum_{𝗸} e^{i𝗸·(𝐑+𝗾_α)} â_{I,𝗸},
$$
obtaining:

$$
Ĥ = \frac{1}{N} \sum_{IJ,𝐑𝐑'} h_{IJ,𝘁} \sum_{𝗸𝗸'} e^{-i𝗸·(𝐑+𝗾_α)} e^{i𝗸'·(𝐑'+𝗾_β)}
â^†_{I,𝗸} â_{J,𝗸'} \\
= \frac{1}{N} \sum_{IJ,𝘁,𝗸𝗸'} h_{IJ,𝘁} \left[ \sum_{𝐑'} e^{-i(𝗸-𝗸')·𝐑'} \right] 
e^{-i𝗸·(𝘁+𝗾_α)} e^{i𝗸'·𝗾_β} â^†_{I,𝗸} â_{J,𝗸'} \\
= \sum_{IJ,𝘁,𝗸𝗸'} h_{IJ,𝘁} δ_{𝗸𝗸'} e^{-i𝗸·𝗾_α} e^{i𝗸'·(𝐭 + 𝗾_β)} â^†_{I,𝗸} â_{J,𝗸'} \\
= \sum_{IJ,𝘁,𝗸} h_{IJ,𝘁} e^{i𝗸·(𝘁+𝗾_β-𝗾_α)} â^†_{I,𝗸} â_{J,𝗸} \\
= \sum_{IJ,𝗸} h_{IJ,𝗸} â^†_{I,𝗸} â_{J,𝗸},
$$
where he have defined: $h_{IJ,𝗸} = \sum_𝘁 h_{IJ,𝘁} e^{i𝗸·(𝘁+𝗾_β-𝗾_α)}$.

These creation and annihilation operators will correspond to a basis set of functions such as
$\ket{φ_{I,𝗸}} = â^†_{I,𝗸} \ket{0}$.
Then if we know how this basis functions transform under the symmetries of an space group,
we can deduce how the Hamiltonian will transform.

## Representation of unitary operations using a basis

Following the deductions made by Barry in Ref. [1].

Let us start with a basis set in real space $\{ψ_{iα}(𝗿)\}$, where $i$ indicates the internal
degrees of freedom of the orbital, $α$ indicates the site $𝗾_α$ inside the Wyckoff position.
Notice that by construction we assume each function $ψ_{iα}(𝗿)$ is localized on $𝗾_α$.
Intuitively, you can think of them as Wannier functions.

We are going to focus our attention to a particular orbital $ψ_{i1}(𝗿)$ localized in the site
$𝗾₁ ≡ 𝗾$. 
This orbital will transform under the representation $ρ$ of the site-symmetry group $G_𝗾$,
associated to $𝗾$. Then, for each $h ∈ G_𝗾$:

$$
h ψ_{i1}(𝗿) = [ρ(h)]_{ji} ψ_{j1}(𝗿)
$$

> [!NOTE]
> We know this because it is tabulated and, in particular, we will obtain it by 
> `NewBandRep.siteir`.

Within the primitive unit cell, an orbital localized on each $𝗾_α$ can be defined as:

$$
ψ_{iα}(𝗿) = g_α ψ_{i1}(𝗿) = ψ_{i1}(g_α⁻¹ 𝗿),
$$
where $g_α$, with translations, generates the coset decomposition of $G_𝗾$ in $G$.
In other words, we can assign for each $𝗾_α$ a space group element $g ∈ G$, such that
$𝗾_α = g_α 𝗾$ and:

$$
G = \bigcup_{α=1}^n g_α (G_𝗾 \ltimes T).
$$

By extension, translated counterparts in other unit cells can be defined by:

$$
\{E|𝘁\} ψ_{iα}(𝗿) = ψ_{iα}(𝗿-𝘁),
$$
where $𝘁$ is a lattice translation.
The set of $n × \text{dim}(ρ) × 𝒩$ functions $ψ_{iα}(𝗿-𝘁)$, where $𝒩$ is
the number of unit cells of the system, will be the basis set on which the induced representation
$D$ will act.

Specifically, given $g = \{ R|τ \} ∈ G$, the coset decomposition implies that for each $g_α$,
there is an unique operation $g_β$ such that:

$$
g g_α = \{E|𝘁_{αβ}\} g_β h,
$$
where $h ∈ G_𝗾$ and $𝘁_{αβ} ≡ (g 𝗾_α) - 𝗾_β$.

> [!NOTE]
> Maybe I can prove this in an appendix just for completeness.

Taking all of this into consideration, we can deduce how our basis set will transform under
the action of every $g ∈ G$:

$$
g ψ_{iα}(𝗿-𝘁) = g \{E|𝘁\} ψ_{iα}(𝗿) = \\
\{E|R𝘁\} g ψ_{iα}(𝗿) = \\
\{E|R𝘁\} \{E|𝘁_{αβ}\} g_β h ψ_{i1}(𝗿) = \\
\{E|R𝘁\} \{E|𝘁_{αβ}\} g_β [ρ(h)]_{ji} ψ_{j1}(𝗿) = \\
\{E|R𝘁\} \{E|𝘁_{αβ}\} [ρ(h)]_{ji} ψ_{jβ}(𝗿) = \\
\{E|R𝘁\} [ρ(h)]_{ji} ψ_{jβ}(𝗿-𝘁_{αβ}) = \\
[ρ(h)]_{ji} ψ_{jβ}(𝐑-R𝘁-𝘁_{αβ})
$$

While it is natural to define the representation in real space, it will more useful to view
it in reciprocal space. This is more evident when $𝒩 → ∞$. 
To this end, we define the Fourier transform of our basis:

$$
φ_{iα,𝗸}(𝗿) = \sum_𝘁 e^{i(𝗸·(𝘁+𝗾_α))} ψ_{iα}(𝗿-𝘁),
$$
where the sum is over all lattice vectors $𝘁 ∈ T$.

> [!NOTE]
> Notice that this is just convention but for building a TB Hamiltonian this choice is
> better since we will eliminate all local phases.

The Fourier transform amounts to a unitary transformation that exchanges $𝒩$ unit
cells for $𝒩$ distinct $𝗸$ points.
The action of $g = \{ R|τ \} ∈ G$ in reciprocal space becomes:

$$
g φ_{iα,𝗸}(𝗿) = \sum_𝘁 e^{i(𝗸·(𝘁+𝗾_α))} g ψ_{iα}(𝗿-𝘁) = \\
\sum_𝘁 e^{i(𝗸·(𝘁+𝗾_α))} [ρ(h)]_{ji} ψ_{jβ}(𝗿-R𝘁-𝘁_{αβ}) = \\
\sum_{𝘁'} e^{i𝗸·(R⁻¹𝘁'+𝗾_β-τ)} [ρ(h)]_{ji} ψ_{jβ}(𝗿-𝘁') = \\
e^{-i([R⁻¹]ᵀ𝗸)·τ} [ρ(h)]_{ji} \sum_{𝘁'} e^{i([R⁻¹]ᵀ𝗸)·(𝘁'+𝗾_β)} ψ_{jβ}(𝗿-𝘁') = \\
e^{-i([R⁻¹]ᵀ𝗸)·τ} [ρ(h)]_{ji} φ_{jβ,[R⁻¹]ᵀ𝗸}(𝗿),
$$
where we have made the substitution: $𝘁' = R𝘁 + 𝘁_{αβ} = R𝘁 + g𝗾_α - 𝗾_β = R𝘁 + R𝗾_α + τ - 𝗾_β
= R (𝘁+𝗾_α) + τ - 𝗾_β ⇒ (𝘁+𝗾_α) = R⁻¹ (𝘁'+𝗾_β-τ)$.

> [!NOTE]
> Important note, we have used here an interesting trick that has been a source of confusion.
> Here, we made the substitution $𝗸·(R⁻¹𝗿) ≡ (g 𝗸)·𝗿$. Let me prove that:
> $$
> 𝗸·(R⁻¹𝗿) = \sum_{ij} k_i (R⁻¹_{ij} r_j) = \sum_{ij} (R⁻¹_{ij} k_i) r_j = ([R⁻¹]ᵀ 𝗸)·𝗿 ≡ (g 𝗸)·𝗿,
> $$

In reciprocal space, the matrix representation can be interpreted as a $𝒩 × 𝒩$
matrix of $n\dim(ρ) × n\dim(ρ)$ blocks, each block can be labeled by $𝗸𝗸'$.
Most of the blocks are zero: given $g = \{R|τ\} ∈ G$, there is only one non-zero block in each
row and column, corresponding to $𝗸' = R𝗸$. Mathematically, we can express this as:

$$
g φ_{iα,𝗸}(𝗿) = \sum_{jβ𝗸'} D_{jβ𝗸',iα𝗸}(g) φ_{jβ,𝗸'}(𝗿),
$$
where we have that:

$$
D_{jβ𝗸',iα𝗸}(g) = e^{-i(g𝗸)·τ} ρ_{ji}(h) δ_{g𝗸𝗸'} δ_{(g 𝗾_α)-𝗾_β \mod 𝘁},
$$
where $𝘁 ∈ T$.

We will use the following notation:

$$
Ρ_{jβ,iα}(g) = e^{-i(g𝗸)·τ} ρ_{ji}(h) δ_{g𝗾_α - 𝗾_β \mod 𝘁},
$$
where we skip the dependence on $𝗸$ which will be unnecessary.

We can vectorize the previous equation as:

$$
\boxed{g Φ_𝗸(𝗿) = Ρ^T(g) Φ_{g𝗸}(𝗿)},
$$
where $Φ_𝗸(𝗿)$ is a column vector formed by $\{φ_{iα,𝗸}(𝗿)\}$, and, $Ρ(g)$ is a $n×n$ matrix
of $\dim(ρ)×\dim(ρ)$ blocks, each of them can be labelled by $α,β$.
Most of the blocks are zero: given $g ∈ G$, there is only one non-zero block in each row and
column, corresponding to $g q_α - q_β = 0 \mod 𝘁$ with $𝘁 ∈ T$.

> [!NOTE]
> We pick the previous definition of the matrix in order to have good properties
> of composition. This is due to the fact that:
> $$
> g₁ g₂ Φ_𝗸(𝗿) = Ρ^T(g₁g₂) Φ_{g₁g₂𝗸}(𝗿) \\
> = g₁ Ρ^T(g₂) Φ_{g₂𝗸}(𝗿) = Ρ^T(g₂) Ρ^T(g₁) Φ_{g₁g₂𝗸}(𝗿) \\ 
> ⇒ \boxed{Ρ(g₁g₂) = Ρ(g₁) Ρ(g₂)}
> $$

### Action of unitary operations in a Hamiltonian

#### Quantization of the representations

The quantization of the previous (classical) theory of representations can be written using
"bra-ket" notation as:

$$
ĝ \ket{φ_{I,𝗸}} = Ρ_{JI}(g) \ket{φ_{J,g𝗸}},
$$
where $\ket{φ_{I,𝗸}} ≡ â^†_{I,𝗸} \ket{0}$. Then:

$$
ĝ â^†_{I,𝗸} ĝ⁻¹ \ket{0} = \\
ĝ â^†_{I,𝗸} \ket{0} = \\
Ρ_{JI}(g) â^†_{J,g𝗸} \ket{0} \\
⇒ \boxed{ĝ â^†_{I,𝗸} ĝ⁻¹ = Ρ_{JI}(g) 
â^†_{J,g𝗸}}
$$

> [!WARNING]
> We are using that $ĝ⁻¹ \ket{0} = \ket{0}$, but I think it is a good assumption that
> symmetries do not have any effect on vacuum.

Consequently: 
$$
ĝ â_{I,𝗸} ĝ^† = 
Ρ^*_{JI}(g) â_{J,g𝗸}
⇒ \boxed{ĝ â_{I,𝗸} ĝ⁻¹ = Ρ^*_{JI}(g) â_{J,g𝗸}}
$$

> [!WARNING]
> We are assuming that $ĝ^† = ĝ⁻¹$.


Then, if we want the Hamiltonian to be invariant under the symmetries, we must impose that:

$$
Ĥ = ĝ Ĥ ĝ⁻¹
$$

Then we obtain that:

$$
Ĥ = \sum_{IJ,𝗸} â^†_{I,𝗸} h_{IJ,𝗸} â_{J,𝗸} = \\
ĝ Ĥ ĝ⁻¹ = \sum_{IJ,𝗸} ĝ â^†_{I,𝗸} h_{IJ,𝗸} â_{J,𝗸} ĝ⁻¹ \\
= \sum_{IJ,𝗸} ĝ â^†_{I,𝗸} ĝ⁻¹ h_{IJ,𝗸} ĝ â_{J,𝗸} ĝ⁻¹ \\
= \sum_{IJ,𝗸,I'J'} â^†_{I',g𝗸} Ρ_{I'I}(g) h_{IJ,𝗸} Ρ^*_{J'J}(g) â_{J',g𝗸} \\
= \sum_{IJ,𝗸} â^†_{I,g𝗸} \left[ Ρ(g) H_{𝗸} Ρ^†(g) \right]_{IJ} â_{J,g𝗸},
$$
where we have defined $H_𝗸 ≡ h_{IJ,𝗸}$ and made the substitution $I',J' → I,J$.
Comparing the first and final rows we obtain the following relation for the Hamiltonian to
be invariant under symmetries:

$$
\boxed{H_𝗸 = Ρ(g) H_{g⁻¹𝗸} Ρ^†(g)}
$$

> [!NOTE]
> Notice that we are considering unitary operations, so we ca write the previous equation
> as:
> $$
> \boxed{H_𝗸 = Ρ(g) H_{g⁻¹𝗸} Ρ⁻¹(g)},
> $$
> which is the most usual way to write it.

In particular, if different EBRs are involved, the representation will be block diagonal in
their indices.
This is because EBRs are "closed" under the space group operations.
Then, we can write the previous expression using the EBR indices in the following way:

$$
H^{αβ}_𝗸 = [Ρ(g)]_{αα} H^{αβ}_{g⁻¹𝗸} [Ρ^†(g)]_{ββ}
$$

This implementation is more convenient code-wise, since allow us to separate the problem
into blocks which are more approachable.

## Methodology on how to write a symbolic Hamiltonian in Julia

We will explain in this section how to implement such constraints in a numerical language
such as Julia.

For that purpose, let us consider a term in a general Hamiltonian which describes the
hopping term between two EBRs. 
For the sake of simplicity, let us call them $α: (𝗾|A)$ and $β: (𝘄|B)$, where $𝗾_i$
represent a certain WP in the SG and $A$ and $B$ are two site-symmetry irreps of any
dimension.

We will denote each point in the WPs and each term in the site-symmetry irreps in a similar
fashion:

$$
𝗾: q₁, q₂, …, q_N \\
𝘄: w₁, w₂, …, w_M \\
A: A₁, A₂, …, A_J \\
B: B₁, B₂, …, B_K
$$

As we have discussed previously, in reciprocal space the Hamiltonian term involving those 
EBRs, $Ĥ_{αβ}$ can be written as a sum of matrices, where each one rows denote an orbital
from the "arriving" EBR and the column an orbital from the "departing" EBR.
Because of this, each Hamiltonian term can be described by a matrix of
$(\#𝗾,\text{dim}(A)) × (\#𝘄,\text{dim}(B))$.
Each of its components will be a complex number which depend on the vector $𝗸$ and on some
free-parameters that later on we will adjust to obtain the band spectrum.

In order to obtain such Hamiltonian term in Julia, we will need to do some previous steps 
to codify such a symbolic structure into our numeric language.

The first step we need to do is to find all the possible hopping distances that can be 
found between this two EBRs. 
Obviously that set will be infinite so we need to impose a particular cutoff.
We will impose it by defining a set of representative lattice translation to until where to
fin possible hopping parameters, and later on, adding all of their symmetry partners.
This computation is performed in the function `obtain_symmetry_related_hoppings`.
This function obtains all symmetry related hopping terms and categorize them into closed
orbits.
Those orbits will be linked to symmetry independent TB models, enabling us to separate them.

Inside of one of this orbits, we will find different hopping vectors 
$δs = [δ₁, δ₂, …, δ_n]$.
Each of them will correspond to a hopping term that will be symmetry related to the others.
In particular, due to the Fourier transform picked it will have the form: 
$δ_i = w_j + G_k - q_r$.

> [!IMPORTANT]
> The hopping vector will depend on the Fourier transform picked.


Additionally, each hopping vector can be obtained from different hopping terms.
Then, each of them will store all of the hopping terms that give rise to that hopping vector.
In particular:

$$
δ₁: q_i → w_j + G_k, q_l → w_l + G_n, … \\
δ₂: q_o → w_p + G_r, q_s → w_t + G_z, … \\
\vdots
$$
where $G$ is used to indicate a particular lattice translations.

With this information we are able to numerically codify the Hamiltonian matrix by terms as
we will show in the following.

First, as we know that all the symmetry related hopping distances, for the cutoff assumed,
we can use them to create an abstract vector $𝘃$ which will store the phases that will
appear in the Hamiltonian's term in reciprocal space. Being specific, this vector would like:

$$
𝘃^T_𝗸 = [e^{i𝗸·δ₁}, e^{i𝗸·δ₂}, …, e^{i𝗸·δ_n}]
$$

Additionally, we will need to assign a free-parameter to each orbital hopping term in the 
Hamiltonian matrix (the ones that afterwards we will tune to replicate the band structure).
Then, the vector that will store them will have a length of 
$\text{len}(δs) × \# 𝗾 × \# 𝘄 × \text{dim}(A) × \text{dim}(B)$.
In particular this vector will look like this:

$$
𝘁^T = [t(δ₁), …, t(δ₂), …, t(δ_n)],
$$
where each $t(δ_i)$ represent a collection of free-parameters, one per hopping term inside
the hopping distance $δ_i$.

Then, each term of the Hamiltonian matrix can be written as:

$$
H^{αβ}_{𝗸,ij} = 𝘃^T_𝗸 𝐌^{αβ}_{ij} 𝘁,
$$
where $𝐌^{αβ}_{ij}$ is a numerical matrix that will contain a $1$ in the position $rq$ if
$𝘃_{𝗸,r}$ and $𝘁_q$ are present on the Hamiltonian matrix term $ij$ position.
At the end what we are doing is the encoding the Hamiltonian term as a bilinear form in
Julia.

We will, then, work with a set of matrices {$𝐌^{αβ}$} that will encode the full Hamiltonian
term and will allow us to operate with it.
We will call these matrices *numerical Hamiltonian*.

In the following, we will show how symmetry operations acts on this set of matrices and how
we obtain the constraints they impose on the Hamiltonian term.

### Action of symmetries on the numerical Hamiltonian

As we stipulated at the beginning, we want to obtain a symmetry-restricted TB model.
For that reason, we want to translate the conditions imposed by unitary operations into the
numerical Hamiltonian.
We will start from the condition imposed into the Hamiltonian term:

$$
H^{αβ}_{g𝗸} = [Ρ(g)]_{αα} H^{αβ}_𝗸 [Ρ^†(g)]_{ββ}
$$

Then,

$$
𝘃^T_{g𝗸} 𝐌^{αβ}_{ij} 𝘁 = [Ρ(g)]_{αα} 𝘃^T_𝗸 𝐌^{αβ}_{ij} 𝘁 [Ρ^†(g)]_{ββ}
$$

Since the representations act on different indices as $𝘃$ and $𝘁$, we can permute them
obtaining:

$$
𝘃^T_{g𝗸} 𝐌^{αβ}_{ij} 𝘁 = 𝘃^T_𝗸 [Ρ(g)]_{αα} 𝐌^{αβ}_{ij} [Ρ^†(g)]_{ββ} 𝘁
$$

In order to compare both $𝐌$ matrices, we need to analyze what is $𝘃_{g𝗸}$.
As can be seeing above, the $𝘃$ vector is constructed as: 
$𝘃^T_𝗸 = [e^{i𝗸·δ₁}, e^{i𝗸·δ₂}, …, e^{i𝗸·δ_n}]$, where $\{ δ_i \}$ is a closed orbit.
Then, $𝘃^T_{g𝗸} = [e^{i(g𝗸)·δ₁}, e^{i(g𝗸)·δ₂}, …, e^{i(g𝗸)·δ_n}]$.
As discussed in a note above, we defined $(g𝗸)·𝗿 ≡ 𝗸·(R⁻¹𝗿)$, where $g = \{ R|τ \}$, then:
$𝘃^T_{g𝗸} = [e^{i𝗸·(R⁻¹δ₁)}, e^{i𝗸·(R⁻¹δ₂)}, …, e^{i𝗸·(R⁻¹δ_n)}]$.
Additionally, since $\{ δ_i \}$ is a closed orbit, $𝘃_{g𝗸}$ will be just a permutation of
$𝘃_𝗸$, in other words, $𝘃_{g𝗸} = σ(g) 𝘃_𝗸$, with $σ(g)$ a permutation.
This permutation is constructed in
`_permute_symmetry_related_hoppings_under_symmetry_operation`.

Then, we obtain the following relation:

$$
(σ(g) 𝘃_𝗸)^T 𝐌^{αβ}_{ij} 𝘁 = 𝘃^T_𝗸 [Ρ(g)]_{αα} 𝐌^{αβ}_{ij} [Ρ^†(g)]_{ββ} 𝘁 \\
𝘃^T_𝗸 σ(g)^T 𝐌^{αβ}_{ij} 𝘁 = 𝘃^T_𝗸 [Ρ(g)]_{αα} 𝐌^{αβ}_{ij} [Ρ^†(g)]_{ββ} 𝘁
$$

The, performing some algebra we obtain that:

$$
𝘃^T_𝗸 \left( σ(g)^T 𝐌^{αβ}_{ij} - [Ρ(g)]_{αα} 𝐌^{αβ}_{ij} [Ρ^†(g)]_{ββ} \right) 𝘁 = 0 \\
⇒ \boxed{\left( σ(g)^T 𝐌^{αβ}_{ij} - [Ρ(g)]_{αα} 𝐌^{αβ}_{ij} [Ρ^†(g)]_{ββ} \right) 𝘁 = 0}
$$

This implies that if we compute the nullspace of the previous subtraction, we will obtain
a set of free-parameter vectors that will fulfill the constrains imposed by unitary
operations.

Notice that this set of vectors will be, in general, complex vector, since the matrices
involved will have complex entries.
Then, in order to avoid compilations, we will split our free-parameter vector $𝘁$ into its
real and imaginary part, so we can work only with parameter that are reals.
This is performed in `split_complex`.
For now on, $𝘁^T = [𝘁^T_\text{real} i 𝘁^T_\text{imag}]$

## Time reversal symmetry Θ

Until now we have only focus on unitary operations.
Now we are going to introduce one anti-unitary operation being TRS $Θ$.
This operation needs a separate treatment since it cannot be described using solely matrices.
The usual approach to study TRS is to use the Wigner's time-reversal operator which writes
TRS operation as a compound operation such as:

$$
Θ = U 𝒦₀,
$$
where $U$ is an unitary operator and $𝒦₀$ is the conjugation operator.
Additionally, here we are going to assume that we are working with spinless systems so that
$Θ² = +1$.

> [!NOTE]
> We are going to focus on grey groups, but generalizing this into non-grey magnetic groups
> is possible by considering, instead of $Θ$, $𝒜 = SΘ$ with $S$ an unitary operation.
> However, we are going to leave that case for the future.

### Action of $Θ$ on a physically real basis set

Assume that we depart from a basis set $\{ψ_I(𝗿)\}$, such that it transform trivially
under TRS. In other words:

$$
Θ ψ_{iα}(𝗿) = ψ_{iα}(𝗿)
$$

This starting point is achieved using the function `physically_realify` in Crystalline.
How this function works is explained in [Appendix A](#appendix-a).

Then, following a similar procedure as the one used in this 
[section](#representation-of-unitary-operations-using-a-basis), we obtain:

$$
Θ φ_{I,𝗸}(𝗿) = \sum_𝘁 e^{-i(𝗸·(𝘁+𝗾_α))} Θ ψ_I(𝗿-𝘁) = \\
\sum_𝘁 e^{-i(𝗸·(𝘁+𝗾_α))} ψ_I(𝗿-𝘁) = \\
\sum_𝘁 e^{i((-𝗸)·(𝘁+𝗾_α))} ψ_I(𝗿-𝘁) = \\
φ_{I,-𝗸}(𝗿)
$$

#### Quantization of the action of $Θ$

We can quantize the previous expression as: $\hat{Θ} \ket{φ_{I,𝗸}} = \ket{φ_{I,-𝗸}}$.
Additionally, since $\ket{φ_{I,𝗸}} = â^†_{I,𝗸} \ket{0}$, we obtain the following relation:

$$
\hat{Θ} â^†_{I,𝗸} \hat{Θ}⁻¹ \ket{0} = \hat{Θ} \ket{φ_{I,𝗸}} = \ket{φ_{I,-𝗸}} = â^†_{I,-𝗸} 
\ket{0} \\
⇒ \boxed{\hat{Θ} â^†_{I,𝗸} \hat{Θ}⁻¹ = â^†_{I,-𝗸}}
$$

> [!CAUTION]
> Here we are cannot use that $\hat{Θ}^† = \hat{Θ}⁻¹$, since it is an anti-unitary operator,
> so we need to find another way.

We know that $â_{I,𝗸} \ket{φ_{J,𝗸'}} = δ_{𝗸𝗸'} δ_{IJ} \ket{0}$ and since 
$\hat{θ}² = +1 ⇒ \hat{Θ} = \hat{Θ}⁻¹$, then:

$$
\hat{Θ} â_{I,𝗸} \hat{Θ}⁻¹ \ket{φ_{J,𝗸'}} = \\
\hat{Θ} â_{I,𝗸} \ket{φ_{j,-𝗸'}} = \\
\hat{Θ} \left( δ_{𝗸-𝗸'} δ_{IJ} \ket{0} \right) = \\
δ_{𝗸-𝗸'} δ_{IJ} \ket{0} = â_{I,-𝗸} \ket{φ_{J,𝗸'}} \\
\boxed{\hat{Θ} â_{I,𝗸} \hat{Θ}⁻¹ = â_{I,-𝗸}}
$$

### Action of $Θ$ on the Hamiltonian

Similarly as studied in this [Section](#action-of-unitary-operations-in-a-hamiltonian), we
can deduce the constraints that TRS impose in a TB Hamiltonian. We start from the relation:

$$
Ĥ = \hat{Θ} Ĥ \hat{Θ}⁻¹
$$

Breaking down each side:

$$
\hat{Θ} Ĥ \hat{Θ}⁻¹ = \sum_{IJ,𝗸} \hat{Θ} â^†_{I,𝗸} h_{IJ,𝗸} â_{J,𝗸} \hat{Θ}⁻¹ \\
= \sum_{IJ,𝗸} \hat{Θ} â^†_{I,𝗸} \hat{Θ}⁻¹ h^*_{IJ,𝗸} \hat{Θ} â_{J,𝗸} \hat{Θ}⁻¹ \\
= \sum_{IJ,𝗸} â^†_{I,-𝗸} h^*_{IJ,𝗸} â_{J,-𝗸} = \\
Ĥ = \sum_{IJ,𝗸} â^†_{I,𝗸} h_{IJ,𝗸} â_{J,𝗸}
$$

Obtaining that:

$$
\boxed{H_𝗸 = H^*_{-𝗸}}
$$

### Action of $Θ$ on the numerical Hamiltonian

We depart form the definition of the numerical Hamiltonian made in this
[Section](#methodology-on-how-to-write-a-symbolic-hamiltonian-in-julia):

$$
H^{αβ}_{𝗸,ij} = 𝘃^T_𝗸 𝐌^{αβ}_{ij} 𝘁,
$$

Then imposing the above constraint when TRS is introduced, we find:

$$
𝘃^T_{-𝗸} 𝐌^{αβ}_{ij} 𝘁 = [𝘃^*_𝗸]^T 𝐌^{αβ}_{ij} 𝘁^*,
$$
since $𝐌$ is formed by exclusively real numbers.

Let us start with the right hand side.
As reminder, we know that $𝘃^T_𝗸 = [e^{i𝗸·δ₁}, e^{i𝗸·δ₂}, …, e^{i𝗸·δ_n}]$, where $\{ δ_i \}$
is a closed orbit under the symmetries of the space group.
Now that we are studying TRS, we face $𝘃^T_{-𝗸} = [e^{-i𝗸·δ₁}, e^{-i𝗸·δ₂}, …, e^{-i𝗸·δ_n}]$,
but now $𝘃^T_{-𝗸}$ could not be a permutation of $𝘃^T_𝗸$.
This is the case for space groups were no symmetry relates $δ$ with $-δ$.
We want to keep $𝘃^T_𝗸$ closed under the symmetries we are studying so we need to add the
$δs$ counterparts in the previous cases.
This is performed in the function `add_timereversal_related_orbits!`.

Now that $𝘃^T_𝗸$ is closed under TRS, we can find, again, the permutation that fulfil that:
$𝘃^T_{-𝗸} = σ(Θ) 𝘃^T_𝗸$.
Additionally, we have to study also $[𝘃^*_𝗸]^T = [e^{-i𝗸·δ₁}, e^{-i𝗸·δ₂}, …, e^{-i𝗸·δ_n}]$,
but as you can see the same equation holds for that vector, i.e., $[𝘃^*_𝗸]^T = σ(Θ) 𝘃^T_𝗸$.

Finally, we need to study $𝘁^*$.
As we said previously, we can split this vector into its real and imaginary part using the
function `split_complex`.
Then we obtain the following equation:

$$
H^{αβ}_{𝗸,ij} = 𝘃^T_𝗸 [𝐌^{αβ}_{ij} 𝐌^{αβ}_{ij}] 
\begin{bmatrix}
   𝘁_\text{real} \\
   i 𝘁_\text{imag}
\end{bmatrix},
$$
then:

$$
[H^{αβ}_{𝗸,ij}]^* = [𝘃^*_𝗸]^T [𝐌^{αβ}_{ij} 𝐌^{αβ}_{ij}] 
\begin{bmatrix}
   𝘁_\text{real} \\
   -i 𝘁_\text{imag}
\end{bmatrix} =
[𝘃^*_𝗸]^T [𝐌^{αβ}_{ij} -𝐌^{αβ}_{ij}] 
\begin{bmatrix}
   𝘁_\text{real} \\
   i 𝘁_\text{imag}
\end{bmatrix}
$$

Considering all off the above, we find that:

$$
H^{αβ}_{-𝗸,ij} = (σ(Θ) 𝘃_𝗸)^T [𝐌^{αβ}_{ij} 𝐌^{αβ}_{ij}] 
\begin{bmatrix}
   𝘁_\text{real} \\
   i 𝘁_\text{imag}
\end{bmatrix} = \\
 𝘃^T_𝗸 σ(Θ)^T [𝐌^{αβ}_{ij} 𝐌^{αβ}_{ij}] 
\begin{bmatrix}
   𝘁_\text{real} \\
   i 𝘁_\text{imag}
\end{bmatrix} = \\
[H^{αβ}_{𝗸,ij}]^* = [𝘃^*_𝗸]^T [𝐌^{αβ}_{ij} 𝐌^{αβ}_{ij}] 
\begin{bmatrix}
   𝘁_\text{real} \\
   -i 𝘁_\text{imag}
\end{bmatrix} = \\
𝘃^T_𝗸 σ(Θ)^T [𝐌^{αβ}_{ij} -𝐌^{αβ}_{ij}] 
\begin{bmatrix}
   𝘁_\text{real} \\
   i 𝘁_\text{imag}
\end{bmatrix}
$$

Subtracting both terms, we obtain that:

$$
𝘃^T_𝗸 σ(Θ)^T [0 \; 2𝐌^{αβ}_{ij}] 
\begin{bmatrix}
   𝘁_\text{real} \\
   i 𝘁_\text{imag}
\end{bmatrix} = 0
$$

Then, obtaining the, again, the nullspace of $σ(Θ)^T 𝐌^{αβ}_{ij}$, we will obtain
a set of free-parameter vectors that will fulfill the constrains imposed by TRS.

Then, we want to obtain the intersection of the space of solutions for TRS and the unitary
operations of the space group.
In order to do so, we implement the Zassenhaus' algorithm in the function
`zassenhaus_intersection`, which allow us to obtain a basis set of free-parameters that 
fulfils at the same time the constraints imposed by TRS and all operations of the space
group under study.

## Hermiticity vs. anti-hermiticity

## Appendix A

Assume that we depart from a basis set $Φ^T_𝗸 = [φ_{𝗸,1}, φ_{𝗸,2}, …, φ_{𝗸,n}]$ which 
engenders a representation $Ρ$, namely:

$$
g Φ_𝗸 = Ρ^T(g) Φ_{g𝗸}
$$

Assuming that $Θ g = g Θ$, we obtain that

$$
g (Θ Φ_𝗸) = Θ (g Φ_{g𝗸}) = Θ (Ρ^T(g) Φ_{g𝗸}) = [P^T(g)]^* (Θ Φ_{g𝗸})
$$

Then the time reversed representation is *identical* to the complex representation $Ρ^*$.

If we now construct the representation Γ engendered from the combined basis
$ℱ^T_𝗸 = [Φ^T_𝗸 (ΘΦ_𝗸)^T]$,
we obtain:

$$
g ℱ_𝗸 = Γ^T(g) ℱ_{g𝗸} = 
\begin{pmatrix}
   Ρ^T(g) & 0 \\
   0 & [P^*(g)]^T
\end{pmatrix}
\begin{bmatrix}
   Φ_{g𝗸} \\
   Θ Φ_{g𝗸}
\end{bmatrix} \\
⇒ \boxed{Γ(g) = 
\begin{pmatrix}
   Ρ(g) & 0 \\
   0 & P^*(g)
\end{pmatrix}}
$$

Next, we apply the same operation to $Θ$, obtaining:

$$
Θ (Θ Φ_𝗸) = Θ² Φ_𝗸 = Φ_𝗸,
$$
then:

$$
Θ ℱ_𝗸 = Γ^T(Θ) ℱ_𝗸 = 
\begin{pmatrix}
   0 & 1 \\
   1 & 0
\end{pmatrix}
\begin{bmatrix}
   Φ_𝗸 \\
   Θ Φ_𝗸
\end{bmatrix} \\
⇒ \boxed{Γ(Θ) = 
\begin{pmatrix}
   0 & 1 \\
   1 & 0
\end{pmatrix}}
$$

Having this result, three different scenarios arises:

1. $Θ Φ_𝗸$ reproduces the same set as $Φ_𝗸$. Then, the co-representation $Γ$ correspond to
   a single representation $Ρ$. In this case no new degeneracy is introduced. This case is
   usually called *real* representation.
2. $Θ Φ_𝗸$ doesn't reproduce the same set as $Φ_𝗸$, but also form a basis for Ρ. Then, the
   co-representation correspond again to a single representation $Ρ$, but with twice the
   dimension. In this case the dimension of $Ρ$ is doubled. This case is usually called
   *pseudo-real* representation. This case will not happen in site-symmetry groups so we
   will not experience them - it needs even dimension.
3. $Θ Φ_𝗸$ doesn't reproduce the same set as $Φ_𝗸$ and forms a basis for a representation
   $Ρ'$ which is not equivalent to $Ρ$. In this case the co-representation will correspond
   to two different representations. This implies that the anti-unitary cause the two
   representations to become degenerate. This case is usually called *complex*
   representation.

In order to always be in the first case, which is the most convenient one, we previously
realify the representations.
This operation is implemented in Crystalline in the function `realify`.

As we discussed above, there are three possible scenarios.
We are going to analyze each one of them here separately.

### Scenario 1: the representation is real

In this case, nothing is necessary to do since the co-representation will coincide with the
representation itself.

### Scenario 2: the representation is pseudo-real

In this case the co-representation will correspond to just doubling the representation itself.
Namely:

$$
Γ(g) =
\begin{pmatrix}
   Ρ(g) & 0 \\
   0 & Ρ(g)
\end{pmatrix}
$$