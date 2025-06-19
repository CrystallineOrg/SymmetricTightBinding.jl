# Documentation for the package `TETB.jl`

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

## The problem we want to solve

Basically the problem we want to solve can be written as: $A 𝗰^t = 𝗺$, where $A$ is the
matrix of BRs in the SG and $𝗰^t$ is a vector of coefficients in the BRs.
It's important to have in mind that $A 𝗰^t = 𝗻^t$.
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
𝗺_{-Γ} \Rightarrow A_{-Γ} 𝗰^{t+l} = 
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

## Representation of the SG operations in 𝗸-space

As discussed in Section 3 of 
[Barry's article](https://doi.org/10.1146/annurev-conmatphys-041720-124134), we know Wannier
functions transform in the following way:

$$
\rho_G(g)a_{i\alpha}(\mathbf{r}-\mathbf{t}) 
=
 \rho_G(g) \{E|\mathbf{t}\} a_{i\alpha}(\mathbf{r}) 
\\
= \{E|R\mathbf{t}\} \{E|\mathbf{t}_{\beta\alpha}\} g_\beta \rho(h) g_\alpha^{-1} 
    a_{i\alpha}(\mathbf{r}) 
\\
= \{E|R\mathbf{t}\} \{E|\mathbf{t}_{\beta\alpha}\} g_\beta 
    \rho(h) a_{i1}(\mathbf{r}) 
\\
= \{E|R\mathbf{t}\} \{E|\mathbf{t}_{\beta\alpha}\}
    g_\beta [\rho(h)]_{ji} a_{j1}(\mathbf{r}) 
\\
= [\rho(h)]_{ji} \{E|R\mathbf{t}\} 
    a_{j\beta}(\mathbf{r}-\mathbf{t}_{\beta\alpha}) 
\\
= [\rho(h)]_{ji} 
    a_{j\beta}(\mathbf{r}-R\mathbf{t}-\mathbf{t}_{\beta\alpha})
$$

being $g \in G$ and $h \in G_{\mathbf{q}}$, such that given a coset decomposition $G = 
\bigcup_\alpha g_\alpha G_\mathbf{q}$ : $g = \{E|\mathbf{t}_{\beta\alpha}\} g_\beta h 
g_\alpha^{-1}$ and $\mathbf{t}_{\beta\alpha} = g\mathbf{q}_\alpha - \mathbf{q}_\beta$. 
We assume we know a site-symmetry group representation: $h a_{i1}(\mathbf{r}) = [\rho(h)]_{ji} 
a_{j1}(\mathbf{r})$, where we are assuming Einstein's summation rule.

Considering the following Fourier transform we want to study how Bloch functions will transform:

$$
a_{i\alpha}(\mathbf{k},\mathbf{r}) = 
\frac{1}{\sqrt{N}} \sum_{\mathbf{t}} e^{i\mathbf{k} \cdot (\mathbf{t}+\mathbf{q}_\alpha)}
a_{i\alpha}(\mathbf{r}-\mathbf{t})
$$

Then:

$$
\rho_G(g) a_{i\alpha}(\mathbf{k},\mathbf{r}) = 
\frac{1}{\sqrt{N}} \sum_{\mathbf{t}} e^{i\mathbf{k} \cdot (\mathbf{t}+\mathbf{q}_\alpha)} 
\rho_G(g) a_{i\alpha}(\mathbf{r}-\mathbf{t}) = 
\\ 
\frac{1}{\sqrt{N}} \sum_{\mathbf{t}} e^{i\mathbf{k} \cdot (\mathbf{t}+\mathbf{q}_\alpha)}
[\rho(h)]_{ji} a_{j\beta}(\mathbf{r}-R\mathbf{t}-\mathbf{t}_{\beta\alpha})
$$

If we now make the substitution: $\mathbf{t}^\prime = R\mathbf{t} + \mathbf{t}_{
\beta\alpha} \Rightarrow \mathbf{t} = R^{-1}\mathbf{t}^\prime - \mathbf{q}_\alpha + R^{-1}
\mathbf{q}_\beta - R^{-1}\mathbf{v}$, then:

$$
\rho_G(g) a_{i\alpha}(\mathbf{k},\mathbf{r})
= \frac{1}{\sqrt{N}} \sum_{\mathbf{t}^\prime} e^{i\mathbf{k} \cdot 
R^{-1}(\mathbf{t}^\prime+\mathbf{q}_\beta - \mathbf{v})} [\rho(h)]_{ji} 
a_{j\beta}(\mathbf{r}-\mathbf{t}^\prime) = 
\\ 
e^{-i (R\mathbf{k})\cdot\mathbf{v}} [\rho(h)]_{ji} a_{j\beta}(R\mathbf{k},\mathbf{r})
$$

This functionality is implemented in the function `sgrep_induced_by_siteir_excl_phase`.
Note, however, that the global phase $e^{-i (R\mathbf{k})\cdot\mathbf{v}}$ is ignored in
`sgrep_induced_by_siteir_excl_phase`; to include the phase, e.g., for use in assessing
the symmetry transformation of a single Bloch state, see `sgrep_induced_by_siteir`.

## Hamiltonian formalism with translational symmetries

Suppose we have an imaginary Bravais lattice with a set of translations denoted by $\mathcal{T}$.
The translations inside such set will be denoted in capital letters $\mathbf{R}$. Additionally,
we may consider that we have different atomic positions inside the unit cell. Those positions
will be denoted by Greek letter $\alpha$ and their positions in the unit cell will be 
$\mathbf{q}_\alpha$. At first, we aren't going to differentiate between atomic positions 
belonging to different WPs or any integral degrees of freedom. Then the most general
Hamiltonian of such Bravais lattice can be written as:

$$
\mathcal{H} = \sum_{\mathbf{r}\mathbf{R},\alpha\beta} t_{\alpha\beta,\mathbf{R}} 
a_{\alpha,\mathbf{r}}^\dagger a_{\beta,\mathbf{r-R}} + \text{c.c.}
$$

We can consider on applying a Fourier transformation to this Hamiltonian, considering the 
translational invariance of it. Let me consider the following Fourier transform of the fields:

$$
a_{\alpha,\mathbf{k}} = \frac{1}{\sqrt{N}} \sum_{\mathbf{t}} e^{i\mathbf{k}\cdot
(\mathbf{t}+\mathbf{q}_\alpha)} a_{\alpha,\mathbf{r-t}}
$$

Then the Hamiltonian in $k$-space will look like:

$$
\mathcal{H} = \sum_{\mathbf{kR},\alpha\beta} t_{\alpha\beta,\mathbf{R}} e^{i\mathbf{k}\cdot
(\mathbf{q}_\alpha-\mathbf{q}_\beta-\mathbf{R})} a_{\alpha,\mathbf{k}}^\dagger 
a_{\beta,\mathbf{k}}
$$
where we have used that $\sum_\mathbf{r} e^{i(\mathbf{k-k'}) \cdot \mathbf{r}} = N 
\delta_{\mathbf{k}\mathbf{k'}}$.

Obviously this will be an infinite problem due to $\mathbf{R}$, so in order to obtain it
numerically, we will need to impose a clever cutoff.

### Our approach: Provide a list of translations $\mathbf{R}$ to consider

Our approach will be to give a set of lattice translations that we want to consider and try
to compute all symmetry related hopping that involve such translation. To illustrate this 
approach, let me consider on lattice translation that I will denote by $\mathbf{R}$. Then 
given to points $\alpha$ and $\beta$, the hopping distance between those two points will be
$\Delta_{\alpha\to \beta} = \mathbf{q}_\beta + \mathbf{R} - \mathbf{q}_\alpha$. Notice that 
this distance is exactly the term in the exponential of the Hamiltonian written in $k$-space.

Then if we want to obtain a symmetry preserving Hamiltonian, we need to include in our study 
at least all other hopping terms related by symmetry. This could be easily obtained by 
considering the symmetries of our system. Let us assume that $g=\{R|\tau\}$ is a symmetry of 
the system. Then the symmetry related hopping terms will be given by $\Delta^\prime
_{\alpha\to \beta} = g(\mathbf{q}_\beta + \mathbf{R}) - g\mathbf{q}_\alpha = R(
\mathbf{q}_\beta + \mathbf{R}) - R\mathbf{q}_\alpha = R\Delta_{\alpha\to \beta}$. Then if we
want to obtain a symmetry constrained Hamiltonian we will need to consider classes of hopping
distances $\{\Delta_{\alpha\to \beta}\}$, which are invariant under the space group symmetries.

Considering all this hopping terms our Hamiltonian will be closed to respect to the symmetry
operations so we can be certain that after imposing the symmetry constraints using the 
representations of the operations as:

$$ \rho_G(g) H(\mathbf{k}) \rho^{-1}_G(g) = H(g\mathbf{k}) $$
our Hamiltonian will be the most general symmetry constrained Hamiltonian up to hopping terms
considering a set of lattice translations.

## Example on 1D bipartite lattice with inversion

Assume we have two sites in a one dimensional lattice of parameter $a=1$ where we place an
inversion-even orbital at the origin denoted by (1a|A); and an inversion-odd orbital at $x=1/2$
denoted by (1b|B).

**Notation:** for simplicity we will denote orbitals (1a|A) as $a$ and (1b|B) as $b$.
Additionally, we will indicate with a subindex the unit cell it belongs to. For example, 
$a_0$ will be placed at $x=0$, while $b_1$ will be placed at $x=3/2$ or $a_{-1}$ at $x=-1$.

### Deduction by inspection

This orbitals will transform under inversion symmetry in the following way:

$$
\mathcal{I} a_n = a_{-n}; \quad \mathcal{I} b_n = -b_{-n-1}
$$ 

Then the most general inversion-symmetric Hamiltonian for 1-st nearest neighbors will be:

$$
\mathcal{H} = \sum_n t (a_n^\dagger b_n - a_n^\dagger b_{n-1}) + \text{c.c.}
$$

It is easy to check that this Hamiltonian is inversion symmetric.

If we consider the following Fourier transform: $a_n = \frac{1}{\sqrt{N}} \sum_n e^{-ikn} a_k$,
and $b_n = \frac{1}{\sqrt{N}} \sum_n e^{-ik(n+1/2)} b_k$
then the Hamiltonian in $k$-space will look like:

$$
\mathcal{H} = \sum_k 2it \sin(k/2) a_k^\dagger b_k + \text{c.c}
$$

### Deduction from our method

First, let me consider the translation $t=0$. Then, remember that $\Delta_{\alpha\to\beta+R} = 
\mathbf{q}_\beta + \mathbf{R} - \mathbf{q}_\alpha$.

- Possible $\Delta$'s: $\Delta_{a\to a} = \Delta_{b\to b} = 0; \quad \Delta_{a\to b} = 
    -\Delta_{b\to a} = -1/2$.
- Orbits of $\Delta$'s: $\{\Delta_{a\to a}\} = \{\Delta_{b\to b}\} = \{0\}; \quad 
    \{\Delta_{a\to b}\} = \{\Delta_{b\to a}\} = \{1/2,-1/2\}$.

Then the non-symmetrized Hamiltonian for this translation will be:

$$
H_{k,0} = \begin{pmatrix} t_{a\to a} & t_{b\to a}^1e^{ik/2}+t_{b\to a}^2e^{-ik/2} \\ 
t_{a\to b}^1e^{ik/2}+t_{a\to b}^2e^{-ik/2} & t_{b\to b} \end{pmatrix}
$$

Let's proceed now to symmetrize this Hamiltonian. Since we only have inversion we only need 
to check that: $H_{k,0} = \rho(\mathcal{I})H_{k,0}\rho^{-1}(\mathcal{I})$, where:

$$
\rho(\mathcal{I}) = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}
$$

Then, that constraint impose that: $\left\{ \begin{matrix} t_{b\to a}^1 = -t_{b\to a}^2 \\
t_{a\to b}^1 = -t_{a\to b}^2 \end{matrix} \right.$, so the Hamiltonian will look like:

$$
H_{k,0} = \begin{pmatrix} t_{a\to a} & 2it_{b\to a} \sin(k/2) \\ 
2it_{a\to b} \sin(k/2) & t_{b\to b} \end{pmatrix}
$$

Which is exactly the Hamiltonian deduced with the previous method but with onsite terms and 
without hermiticity imposed.

### Implementation of the free-parameter Hamiltonian into Julia code

For a non-symbolic programming language as Julia we cannot implement this Hamiltonian directly
with some free parameters. In order to do so, we will rewrite it in a clever way so we can 
perform linear algebra with it.

Let me start with our previous example. For this let me focus on only one hopping term. In 
this case it will be $b \to a$. The Hamiltonian for this term will look like:

$$
H_{k,0}^{b \to a} = t_{b\to a}^{1/2} e^{ik/2} + t_{b\to a}^{-1/2} e^{-ik/2}
$$

This term can be obtained by just knowing the class $\{\Delta_{a\to b}\}$. Knowing this we 
can define two vectors such that:

$$
\mathbf{v}^T = \begin{pmatrix}
    e^{ik/2} & e^{-ik/2} \\
\end{pmatrix}; \quad \mathbf{t} = \begin{pmatrix}
    t_{b\to a}^{1/2} \\
    t_{b\to a}^{-1/2}
\end{pmatrix}
$$

Then $H_{k,0}^{b \to a}$ can be written as:

$$
H_{k,0}^{b \to a} = \mathbf{v}^T \underbrace{\begin{pmatrix}
    1 & 0 \\
    0 & 1
\end{pmatrix}}_{M^{b \to a}_0} \mathbf{t}
$$

In general, I think, this can be performed by creating an identity matrix of dimension 
$\#\{\Delta_{b\to a}\} \times (\#\{\Delta_{b\to a}\} \times \#\Delta_{b\to a} \times \mu_a 
\times \mu_b)$.

How a symmetry transformation will look like if I apply the symmetries to $M$ instead of $H$.
Then they will look as:

$$ \rho_G(g)_{im} H^{mn}_{k} \rho_G^{-1}(g)_{mj} = v_\alpha \rho_G(g)_{im} M^{mn}_{\alpha\beta} 
    \rho_G^{-1}(g)_{nj} t_\beta $$

So, we can represent the action of $\rho_G(g)H_k\rho_G(g)^{-1}$ as:

$$
    (\rho_G(g)H_k\rho_G(g)^{-1})_{ij}
    =
    v_\alpha R_{ij,\alpha\beta}(g) t_\beta
    \text{ with }
    R_{ij,\alpha\beta}(g) = \rho_G(g)_{im}M_{\alpha\beta}^{mn}\rho_G^{-1}(g)_{nj}
$$

or, in matrix notation,

$$
    \rho_G(g)H_k\rho_G(g)^{-1}
    =
    v_\alpha \mathbf{R}_{\alpha\beta}(g) t_\beta
    \text{ with }
    \mathbf{R}_{\alpha\beta}(g) = \bm{\rho}_G(g)𝗺_{\alpha\beta}\bm{\rho}_G^{-1}(g)
$$

Or, defining $\mathbf{R}$ as a block-matrix with elements $\mathbf{R}_{\alpha\beta}$, we have:

$$
    \rho_G(g)H_k\rho_G(g)^{-1} = \mathbf{v}^T \mathbf{R} \mathbf{t}
$$

The transformation $H(gk)$, on the other hand, becomes:

$$ H(gk)_{ij} = \text{permuted}(\mathbf{v})^T 𝗺^{ij} \mathbf{t} = \mathbf{v}^T 
\text{permuted}(𝗺)^{ij} \mathbf{t} $$

The permutation of $\mathbf{v}$ can be realized by a permutation matrix $𝗽(g)$, s.t.
$\text{permuted}(\mathbf{v}^T) = 𝗽(g)\mathbf{v}$, s.t.:

$$
    \text{permuted}(M)^{ij} = 𝗽(g)^T 𝗺^{ij}
$$

#### How do we obtain such permutation?

Since we have $\{\Delta_{b \to a}\}$, we can just do $R^{-1}\{\Delta_{b\to a}\}$, with $g=\{R|
\mathbf{v}\}$, and find the permutations made in $\{\Delta_{b \to a}\}$. Those permutations 
will be exactly the ones we are looking for. If we apply those permutations to the rows of 
$M$ then the operation will be performed.

**NOTE**: this is done with $R^{-1}$ since $(g\mathbf{k})\cdot \delta = [(R^{-1})^T\mathbf{k}]
\cdot\delta = \mathbf{k}\cdot(R^{-1}\delta)$, where $g=\{R|\mathbf{v}\}$. Then, we need to use 
the transpose of the operation not the actual operation for this trick to ork.


## Methodology on how to write a symbolic Hamiltonian in Julia

Let us know explain in a formal way the strategy followed before.

For that purpose, let us consider a term in a general Hamiltonian which describes the hopping 
term between two EBRs. For the sake of simplicity let us call them $\alpha: (\mathbf{q}|A)$ 
and $\beta: (\mathbf{w}|B)$, where $\mathbf{q}_i$ represent a certain WP in the SG and $A$ 
and $B$ are two site-symmetry irreps of any dimension.

For the sake of notation, we will denote each point in the WPs and each term in the 
site-symmetry irreps in a similar fashion:

$$
\mathbf{q}: q_1, q_2, \dots, q_N \\
\mathbf{w}: w_1, w_2, \dots, w_M \\
A: A_1, A_2, \dots, A_J \\
B: B_1, B_2, \dots, B_K
$$

As we have discussed previously, in reciprocal space the Hamiltonian term involving those 
EBRs, $H_{\alpha\beta}$ can be written as a matrix where each row denote an orbital from the 
"arriving" EBR and the column an orbital from the "departing" EBR. Because of this the 
Hamiltonian term can be described by a $(\#\mathbf{q},\text{dim}(A)) \times 
(\#\mathbf{w},\text{dim}(B))$. Each of its components will be a complex number which depend 
on the vector $\mathbf{k}$ and on some free-parameters that later on we will adjust to obtain 
the band structure.

In order to obtain such Hamiltonian term in `Julia`, we will need to do some previous steps 
so we can code such a symbolic structure into our numeric language.

The first step we need to do is to find all the possible hopping distances that can be 
found between this two EBRs. Obviously that set will be infinite so we need to impose a 
particular cutoff. As explained above, we will impose it by constraining the hopping terms
to a particular set of lattice translations (and obviously theirs symmetry partners). This 
complex structure is computed in the function `obtain_symmetry_related_hoppings`, where we
provide a set of representatives of hopping distances which which is associated to a set of
hopping terms that are symmetry related.

Inside of one of this representatives we will find different hopping distances $\delta s = 
[\delta_1, \delta_2, \dots, \delta_n]$, which will be associated to different hopping terms:

$$
\delta_1: q_i \to w_j + G_k, q_l \to w_l + G_n, \dots \\
\delta_2: q_o \to w_p + G_r, q_s \to w_t + G_z, \dots \\
\vdots
$$

where $G_k$ are some particular lattice translations.

With this information we are able to numerically codify the Hamiltonian matrix by terms as 
we will show in the following.

First, as we know that all the symmetry related hopping distances, for the cutoff assumed, 
are stored in the output of `obtain_symmetry_related_hoppings`, we can use them to create an 
abstract vector $\mathbf{v}$ which will store the phases that will appear in the 
Hamiltonian's term in reciprocal space. Being specific, this vector would like:

$$
\mathbf{v}^T = [e^{i\mathbf{k}\cdot\delta_1}, e^{i\mathbf{k}\cdot\delta_2}, \dots, 
e^{i\mathbf{k}\cdot\delta_n}]
$$

Note that we are going to use here the order provided by the function 
`obtain_symmetry_related_hoppings` to store this phases.

Additionally, we will need to assign a free-parameter to each orbital hopping term in the 
Hamiltonian matrix (the ones that afterwards we will tune to replicate the band structure).
This vector then will have a length of $\text{len}(\delta s) \times \# \mathbf{q} \times \# 
\mathbf{w} \times \text{dim}(A) \times \text{dim}(B)$. In particular this vector will look 
like this:

$$
\mathbf{t}^T = [t(\delta_1), \dots, t(\delta_2), \dots, t(\delta_n)]
$$

where each $t(\delta_i)$ represent a collection of free-parameter. One per hopping term 
inside the hopping distance $\delta_i$.

Then each term of the Hamiltonian matrix can be written as:

$$
H_{\alpha\beta,ij} = \mathbf{v}^T 𝗺_{\alpha\beta,ij} \mathbf{t}
$$

where $𝗺_{\alpha\beta,ij}$ is a numerical matrix that will relate a phase with a 
free-parameter present on the Hamiltonian matrix term. At the end what we are doing is the 
encoding the bilinear form of the Hamiltonian matrix term so we can operate with it in 
Julia.

We will, then, work with a set of matrices {$𝗺_{\alpha\beta}$} that will encode the 
full Hamiltonian term and will allow us to operate with it.

Allow us now to show how symmetry operations acts on this set of matrices and how to obtain 
the constraints they impose on the Hamiltonian term.

### Action of symmetries on the $𝗺$ matrices

