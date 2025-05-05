# Documentation for the package `TETB.jl`

This code is aimed at finding a decomposition in terms of EBRs of the band representation of
the first two photonic bands in a 3D PhC and building a symmetry-constrained TB using the
EBRs to achieve it.
This code could be also interesting to construct TB models for spinless systems with or
without time-reversal symmetry and with hermiticity or anti-hermiticity.

## Notation

Let me first specify the notation we are going to use through out the code in other to avoid
any misunderstandings.
We are going to use Unicode characters and face as best as we can with their limitations.

First, we need to differentiate between symmetry vectors obtained from MPB (or given by the 
user) and the ones obtained from the code as possible candidates for a TETB model.
For making this differentiation we are going to use the following notation:

1. Symmetry vectors obtained from MPB will be denoted by $𝗺$.
2. However symmetry vectors coming as outputs of the code will be denoted by $𝗻$. There will
   be different vectors coming as outputs from the code.
   We will denoted them by a super-index indicating their "polarization".
   Being more specific, the symmetry vector representing the whole TETB model (transverse +
   longitudinal modes) and will be denoted by $𝗻^{t+l}$; and the one representing only the
   longitudinal modes $𝗻^l$. You can obtain the symmetry vector of transversal modes from
   those two by $𝗻^t = 𝗻^{t+l} - 𝗻^l$.

Additionally, it is interesting to separate those symmetry vectors into the irreps belonging 
to Γ and other HSPs.
This is because we need to test if the zero frequency modes are well represented by the model
as stipulated by this [paper](https://link.aps.org/doi/10.1103/PhysRevX.12.021066), and that
any negative multiplicity coming from the ill-definition of that content is present on the
longitudinal modes $𝗻^l$. For this we do the following:

1. Separate all the symmetry vectors into a part belonging to Γ and the other HSPs.
   For example for the symmetry vector coming from MPB: $𝗺 = 𝗺_Γ + 𝗺_{-Γ}$, but this could
   be dome similarly for the other symmetry vectors.
2. Additionally, the part belonging to Γ can be separated into two: 
   $𝗺_Γ = 𝗺_Γ^{=0} + 𝗺_Γ^{>0}$, one part belonging to the modes at $ω = 0$ and the ones at
   $ω > 0$.

Now having the notation clear we can explain how our code works.

## The problem under study

Basically the problem we want to solve can be written as: $A 𝗰^t = 𝗺$, where $A$ is the
matrix of BRs in the SG and $𝗰^t$ is a vector of coefficients in the BRs.
It's important to have in mind that $𝗻^t ≡ A 𝗰^t$.
In the code $A$ will be usually called `brs`.
The notation explained for the symmetry vectors can also be applied to this matrix.
We can divide this problem into two enabling us to study the Γ-point separately:

$$
\begin{bmatrix} A_Γ \\ A_{-Γ} \end{bmatrix} 𝗰^t = 
\begin{bmatrix} 𝗺_Γ \\ 𝗺_{-Γ} \end{bmatrix}
$$

### Problem 1

$$
A_{-Γ} 𝗰^t = 
𝗺_{-Γ} ⇒ A_{-Γ} 𝗰^{t+l} = 
𝗺_{-Γ} + A_{-Γ} 𝗰^l
$$

In general $𝗰^t$ could have some irreps with negative multiplicities, that why we separated
the problem into $𝗰^{t+l}$ and $𝗰^l$ so the problem will be strictly positive.
For more details check this [paper](https://doi.org/10.48550/arXiv.2305.18257).

First, we find all possible longitudinal modes $𝗻^l = A 𝗰^l$ in the SG for a given (usually,
smallest possible) occupation `t`. 
This is performed with the function `find_auxiliary_modes`.

Then, using the identified set of of fixed-occupation $\{ 𝗻^l \}$, we search if there are
transverse+longitudinal modes $𝗻^{t+l}$ that satisfy the previous equation.
This is done in the function `find_apolar_modes`.

### Problem 2

$$
A_{Γ} 𝗰^t = 
𝗺_{Γ} = 
𝗺_Γ^{=0} + 𝗺_Γ^{>0}
$$

`MPBUtils.jl` forces the content at zero frequency to be exactly what is shown in Table 
(S6-8) of this [paper](https://link.aps.org/doi/10.1103/PhysRevX.12.021066) which we are
going to call `n_fixed`, but we have some freedom coming from $Q 𝗽$ which can appear
in our model.
Then, the previous equality should be thought of as `n_fixed`$\mod Q$.
More explicitly, a compatible solution must solve the following equation with 
$𝗽 ∈ ℤ$:

$$
A_{Γ} 𝗰^T - 𝗺_{Γ} =
Q 𝗽
$$

See the details in [Thomas' paper](https://link.aps.org/doi/10.1103/PhysRevX.12.021066).

## Physicality

Description of what a *physical* solutions means.

Assume we have a solution provided by our code, which consist on a longitudinal part $𝗻^l$ 
and a transverse+longitudinal one $𝗻^{t+l}$.
This is obtained by solving [Problem 1](#problem-1). 
Given a symmetry vector $𝗺$ from MPB, we are considering if our solution properly represents 
the content at Γ, which has been in the computation of the solution. For this we 
need to check two things:

1. Whether our solution sub-duces properly to the $O(3)$ representation at Γ and zero 
   frequency. This can be checked easily using `PhotonicBandConnectivity.jl`. As explained 
   before in [Problem 2](#problem-2), this is fulfilled if there exists a $𝗽 ∈ ℤ$ solving 
   $A_{Γ} 𝗰^T - 𝗺_{Γ} = Q 𝗽$.
2. Whether our solution doesn't make use of the higher frequency irreps ($ω > 0$) present in
   $𝗺_Γ^{>0}$ to regularize the symmetry content at $ω = 0$, and that instead those negative
   multiplicities in the irreps are cancelled out by the longitudinal modes $𝗻^l$. 
   We ensure this by the following check:𝗻
    Define the candidate-solution's zero-frequency content at Γ by:

    $$
    𝗻_Γ^{t,=0} = 𝗻_{Γ}^{t} - 𝗻_{Γ}^{>0} = 
    𝗻_{Γ}^{t+l} - 𝗻_{Γ}^l - 𝗻_{Γ}^{>0} = 
    𝗺_{Γ}^{=0} + Q 𝗽
    $$
  
    Consider the following two cases:
    - If $n_{Γ,i}^{t,=0} < 0$ for some $i$, then $n_{Γ,i}^l \geq |n_{Γ,i}^{t,=0}|$ for that
    $i$; equivalently, in this case $n_{Γ,i}^l \geq -n_{Γ,i}^{t,=0}$.
    - Conversely, if  $n_{Γ,i}^{t,=0} ≥ 0$ for some $i$, we still have $n_{Γ,i}^l ≥ 0$ and
    consequently also $n_{Γ,i}^l ≥ -n_{Γ,i}^{t,=0}$.

    Thus, regardless of the sign of $n_{Γ,i}^{t,=0}$, we may require that:

    $$
    n_{Γ}^l \geq -n_Γ^{t,=0}
    $$

These constraints are directly imposed in the function `find_apolar_modes` thanks to the
functionalities of `find_all_admissible_expansions`.

Now that we have an EBR decomposition of our band set, we can use them as a basis to build a
TB Hamiltonian.

## General Hamiltonian in periodic systems

Let us start with the most general non-interacting Hamiltonian:

$$
Ĥ = \sum_{IJ,𝐑𝐑'} h_{IJ,𝐑-𝐑'} \; ĉ_{I,𝐑}^† ĉ_{J,𝐑'},
$$
where $I,J$ wrap up the internal degrees of freedom of the orbitals and the sites, i.e., 
$I = (i, α)$; and $𝐑,𝐑'$ run over the lattice translations.

> [!WARNING]
> Notice that we have assumed that hopping terms only depends on relative distances.
> We are going to denote $𝘁 ≡ 𝐑 - 𝐑'$.

We can apply the same Fourier transform to go into reciprocal space:

$$
ĉ_{I,𝐑} = \frac{1}{\sqrt{N}} \sum_{𝗸} e^{-i𝗸·(𝐑+𝗾_α)} â_{I,𝗸},
$$
obtaining:

$$
Ĥ = \frac{1}{N} \sum_{IJ,𝐑𝐑'} h_{IJ,𝘁} \sum_{𝗸𝗸'} e^{i𝗸·(𝐑+𝗾_α)} e^{-i𝗸'·(𝐑'+𝗾_β)}
â^†_{I,𝗸} â_{J,𝗸'} \\
= \frac{1}{N} \sum_{IJ,𝘁,𝗸𝗸'} h_{IJ,𝘁} \left[ \sum_{𝐑'} e^{i(𝗸-𝗸')·𝐑'} \right] 
e^{i𝗸·(𝘁+𝗾_α)} e^{-i𝗸'·𝗾_β} â^†_{I,𝗸} â_{J,𝗸'} \\
= \sum_{IJ,𝘁,𝗸𝗸'} h_{IJ,𝘁} δ_{𝗸𝗸'} e^{i𝗸·(𝘁+𝗾_α)} e^{-i𝗸'·𝗾_β} â^†_{I,𝗸} â_{J,𝗸'} \\
= \sum_{IJ,𝘁,𝗸} h_{IJ,𝘁} e^{i𝗸·𝗸·(𝘁+𝗾_α-𝗾_β)} â^†_{I,𝗸} â_{J,𝗸} \\
= \sum_{IJ,𝗸} h_{IJ,𝗸} â^†_{I,𝗸} â_{J,𝗸},
$$
where he have defined: $h_{IJ,𝗸} = \sum_𝘁 h_{IJ,𝘁} e^{i𝗸·𝗸·(𝘁+𝗾_α-𝗾_β)}$.

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
$𝗾_1 ≡ 𝗾$. 
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
ψ_{iα}(𝗿) = g_α ψ_{i1}(𝗿) = ψ_{i1}(g_α^{-1} 𝗿),
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
The set of $n × \text{dim}(ρ) × \mathcal{N}$ functions $ψ_{iα}(𝗿-𝘁)$, where $\mathcal{N}$ is
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
it in reciprocal space. This is more evident when $\mathcal{N} → ∞$. 
To this end, we define the Fourier transform of our basis:

$$
φ_{iα,𝗸}(𝗿) = \sum_𝘁 e^{i(𝗸·(𝘁+𝗾_α))} ψ_{iα}(𝗿-𝘁),
$$
where the sum is over all lattice vectors $𝘁 ∈ T$.

> [!NOTE]
> Notice that this is just convention but for building a TB Hamiltonian this choice is
> better since we will eliminate all local phases.

The Fourier transform amounts to a unitary transformation that exchanges $\mathcal{N}$ unit
cells for $\mathcal{N}$ distinct $𝗸$ points.
The action of $g = \{ R|τ \} ∈ G$ in reciprocal space becomes:

$$
g φ_{iα,𝗸}(𝗿) = \sum_𝘁 e^{i(𝗸·(𝘁+𝗾_α))} g ψ_{iα}(𝗿-𝘁) = \\
\sum_𝘁 e^{i(𝗸·(𝘁+𝗾_α))} [ρ(h)]_{ji} ψ_{jβ}(𝗿-R𝘁-𝘁_{αβ}) = \\
\sum_{𝘁'} e^{i𝗸·(R^{-1}𝘁'+𝗾_β-τ)} [ρ(h)]_{ji} ψ_{jβ}(𝗿-𝘁') = \\
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

In reciprocal space, the matrix representation can be interpreted as a $\mathcal{N}×\mathcal{N}$
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
ĝ â^†_{I,𝗸} ĝ^{-1} \ket{0} = \\
ĝ â^†_{I,𝗸} \ket{0} = \\
Ρ_{JI}(g) â^†_{J,g𝗸} \ket{0} \\
⇒ \boxed{ĝ â^†_{I,𝗸} ĝ^{-1} = Ρ_{JI}(g) 
â^†_{J,g𝗸}}
$$

> [!WARNING]
> We are using that $ĝ^{-1} \ket{0} = \ket{0}$, but I think it is a good assumption that
> symmetries do not have any effect on vacuum.

Consequently: 
$$
ĝ â_{I,𝗸} ĝ^† = 
Ρ^*_{JI}(g) â_{J,g𝗸}
⇒ \boxed{ĝ â_{I,𝗸} ĝ^{-1} = Ρ^*_{JI}(g) â_{J,g𝗸}}
$$

> [!WARNING]
> We are assuming that $ĝ^† = ĝ^{-1}$.


Then, if we want the Hamiltonian to be invariant under the symmetries, we must impose that:

$$
Ĥ = ĝ Ĥ ĝ^{-1}
$$

Then we obtain that:

$$
Ĥ = \sum_{IJ,𝗸} â^†_{I,𝗸} h_{IJ,𝗸} â_{J,𝗸} = \\
ĝ Ĥ ĝ^{-1} = \sum_{IJ,𝗸} ĝ â^†_{I,𝗸} h_{IJ,𝗸} â_{J,𝗸} ĝ^{-1} \\
= \sum_{IJ,𝗸} ĝ â^†_{I,𝗸} ĝ^{-1} h_{IJ,𝗸} ĝ â_{J,𝗸} ĝ^{-1} \\
= \sum_{IJ,𝗸,I'J'} â^†_{I',g𝗸} Ρ_{I'I}(g) h_{IJ,𝗸} Ρ^*_{J'J}(g) â_{J',g𝗸} \\
= \sum_{IJ,𝗸} â^†_{I,g𝗸} \left[ Ρ(g) H_{𝗸} Ρ^†(g) \right]_{IJ} â_{J,g𝗸},
$$
where we have defined $H_𝗸 ≡ h_{IJ,𝗸}$ and made the substitution $I',J' → I,J$.
Comparing the first and final rows we obtain the following relation for the Hamiltonian to
be invariant under symmetries:

$$
\boxed{H_𝗸 = Ρ(g) H_{g^{-1}𝗸} Ρ^†(g)}
$$

> [!NOTE]
> Notice that we are considering unitary operations, so we ca write the previous equation
> as:
> $$
> \boxed{H_𝗸 = Ρ(g) H_{g^{-1}𝗸} Ρ⁻¹(g)},
> $$
> which is the most usual way to write it.

In particular, if different EBRs are involved, the representation will be block diagonal in
their indices.
This is because EBRs are "closed" under the space group operations.
Then, we can write the previous expression using the EBR indices in the following way:

$$
H^{αβ}_𝗸 = [Ρ(g)]_{αα} H^{αβ}_{g^{-1}𝗸} [Ρ^†(g)]_{ββ}
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
𝗾: q_1, q_2, …, q_N \\
𝘄: w_1, w_2, …, w_M \\
A: A_1, A_2, …, A_J \\
B: B_1, B_2, …, B_K
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
$δs = [δ_1, δ_2, …, δ_n]$.
Each of them will correspond to a hopping term that will be symmetry related to the others.
In particular, due to the Fourier transform picked it will have the form: 
$δ_i = w_j + G_k - q_r$.

> [!CAUTION]
> The hopping vector will depend on the Fourier transform picked.


Additionally, each hopping vector can be obtained from different hopping terms.
Then, each of them will store all of the hopping terms that give rise to that hopping vector.
In particular:

$$
δ_1: q_i → w_j + G_k, q_l → w_l + G_n, … \\
δ_2: q_o → w_p + G_r, q_s → w_t + G_z, … \\
\vdots
$$
where $G$ is used to indicate a particular lattice translations.

With this information we are able to numerically codify the Hamiltonian matrix by terms as
we will show in the following.

First, as we know that all the symmetry related hopping distances, for the cutoff assumed,
we can use them to create an abstract vector $𝘃$ which will store the phases that will
appear in the Hamiltonian's term in reciprocal space. Being specific, this vector would like:

$$
𝘃^T_𝗸 = [e^{i𝗸·δ_1}, e^{i𝗸·δ_2}, …, e^{i𝗸·δ_n}]
$$

Additionally, we will need to assign a free-parameter to each orbital hopping term in the 
Hamiltonian matrix (the ones that afterwards we will tune to replicate the band structure).
Then, the vector that will store them will have a length of 
$\text{len}(δs) × \# 𝗾 × \# 𝘄 × \text{dim}(A) × \text{dim}(B)$.
In particular this vector will look like this:

$$
𝘁^T = [t(δ_1), …, t(δ_2), …, t(δ_n)],
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
$𝘃^T_𝗸 = [e^{i𝗸·δ_1}, e^{i𝗸·δ_2}, …, e^{i𝗸·δ_n}]$, where $\{ δ_i \}$ is a closed orbit.
Then, $𝘃^T_{g𝗸} = [e^{i(g𝗸)·δ_1}, e^{i(g𝗸)·δ_2}, …, e^{i(g𝗸)·δ_n}]$.
As discussed in a note above, we defined $(g𝗸)·𝗿 ≡ 𝗸·(R^{-1}𝗿)$, where $g = \{ R|τ \}$, then: $𝘃^T_{g𝗸} = [e^{i𝗸·(R^{-1}δ_1)}, e^{i𝗸·(R^{-1}δ_2)}, …, e^{i𝗸·(R^{-1}δ_n)}]$.
Additionally, since $\{ δ_i \}$ is a closed orbit, $𝘃_{g𝗸}$ will be just a permutation of
$𝘃_𝗸$, in other words, $𝘃_{g𝗸} = σ 𝘃_𝗸$, with $σ$ a permutation.
This permutation is constructed in
`_permute_symmetry_related_hoppings_under_symmetry_operation`.

Then, we obtain the following relation:

$$
(σ 𝘃_𝗸)^T 𝐌^{αβ}_{ij} 𝘁 = 𝘃^T_𝗸 [Ρ(g)]_{αα} 𝐌^{αβ}_{ij} [Ρ^†(g)]_{ββ} 𝘁 \\
𝘃^T_𝗸 σ^T 𝐌^{αβ}_{ij} 𝘁 = 𝘃^T_𝗸 [Ρ(g)]_{αα} 𝐌^{αβ}_{ij} [Ρ^†(g)]_{ββ} 𝘁
$$

The, performing some algebra we obtain that:

$$
𝘃^T_𝗸 \left( σ^T 𝐌^{αβ}_{ij} - [Ρ(g)]_{αα} 𝐌^{αβ}_{ij} [Ρ^†(g)]_{ββ} \right) 𝘁 = 0 \\
⇒ \boxed{\left( σ^T 𝐌^{αβ}_{ij} - [Ρ(g)]_{αα} 𝐌^{αβ}_{ij} [Ρ^†(g)]_{ββ} \right) 𝘁 = 0}
$$

This implies that if we compute the nullspace of the previous subtraction, we will obtain
a set of free-parameter vectors that will fulfill the constrains imposed by unitary
operations.

Notice that this set of vectors will be, in general, complex vector, since the matrices
involved will have complex entries.
Then, in order to avoid compilations, we will split our free-parameter vector $𝘁$ into its
real and imaginary part.
This is performed in `split_complex`.

## Time reversal symmetry

## Hermiticity vs. anti-hermiticity