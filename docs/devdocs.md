# Documentation for the package `TETB.jl`

## Notation

Let me first specify the notation we are going to use through out the code in other to avoid
any misunderstandings. We are going to use Unicode characters and face as best as we can with 
their limitations.

First, we need to differentiate between symmetry vectors obtained from MPB (or given by the 
user) and the ones obtained from the code as possible candidates for a TETB model. For making 
this differentiation we are going to use the following notation:

1. Symmetry vectors obtained from MPB will be denoted by $m$.
2. However symmetry vectors coming as outputs of the code will be denoted by $n$. There will 
   be different vectors coming as outputs from the code. We will denoted them by a super-index 
   indicating their "polarization". Being mere specific, the symmetry vector representing the 
   whole TETB model (transverse + longitudinal modes) and will be denoted by $n^{T+L}$; and 
   the one representing only the longitudinal modes $n^L$. You can obtain the symmetry vector 
   of transversal modes from those two by $n^T=n^{T+L}-n^L$.

Additionally, it is interesting to separate those symmetry vectors into the irreps belonging 
to $\Gamma$ and the others. This is because we need to test if the zero frequency modes are 
well represented by the model as stipulated by 
[Thomas' paper](https://link.aps.org/doi/10.1103/PhysRevX.12.021066), and that any negative 
multiplicity coming from the ill-definition of that content is present on the longitudinal 
modes $n^L$. For this we do the following:

1. Separate all the symmetry vectors into a part belonging to $\Gamma$ and the other high 
   symmetry points. For example for the symmetry vector coming from MPB: 
   $m = m_\Gamma + m_{-\Gamma}$, but this could be dome similarly for the other symmetry vectors.
2. Additionally, the part belonging to $\Gamma$ can be separated into two: 
   $m_\Gamma=m_\Gamma^{=0}+m_\Gamma^{>0}$, one part belonging to the modes at $\omega=0$ and 
   the ones at $\omega>0$.

Now having the notation clear we can explain how our code works.

## The problem we want to solve

Basically the problem we want to solve can be written as: $Ac^T=m$, where $A$ is the matrix 
of BRs in the SG and $c$ is a vector of coefficients in the BRs. It's important to have in mind
that $Ac^T=n^T$. In the code $A$ will be called `brs`. The notation explained for the symmetry 
vectors can also be applied to this matrix. We can divide this problem into two enabling us 
to study the $\Gamma$-point separately:

$$ \begin{bmatrix} A_\Gamma \\ A_{-\Gamma} \end{bmatrix} c^T = 
    \begin{bmatrix} m_\Gamma \\ m_{-\Gamma} \end{bmatrix} $$

### Problem 1

$$ A_{-\Gamma} c^T = m_{-\Gamma} \Rightarrow A_{-\Gamma} c^{T+L} = m_{-\Gamma} + A_{-\Gamma}
    c^L $$

In general $n^T$ could have some irreps with negative multiplicities, that why we separated 
the problem into $n^{T+L}$ and $n^L$ so the problem will be strictly positive. For more details 
check [Antonio's paper](https://doi.org/10.48550/arXiv.2305.18257).

First we find all possible longitudinal modes $n^L$ in the SG for a given (usually, smallest 
possible) occupation `t`. 
This is performed with the function `find_auxiliary_modes`.

Then, using the identified set of of fixed-occupation $\{n^L\}$, we search if there are 
transverse+longitudinal modes $n^{T+L}$ that satisfy the previous equation. This is done in 
the function `find_apolar_modes`.

### Problem 2

$$ A_{\Gamma} c^T = m_{\Gamma} = m_\Gamma^{=0} + m_\Gamma^{>0} $$

`MPBUtils.jl` forces the content at zero frequency to be exactly that shown in Table 
(S6-8) of [Thomas' paper](https://link.aps.org/doi/10.1103/PhysRevX.12.021066) which we are 
going to call `n_fixed`, but we have some freedom coming from $Q\mathbf{p}$ which can appear 
in our model. Then, the previous equality should be thought of as `n_fixed`$\mod Q$.
More explicitly, a compatible solution must solve the following with some $\mathbf{p}\in\mathbb{Z}$:

$$ A_{\Gamma} c^T - m_{\Gamma} = Q\mathbf{p}. $$

See the details in [Thomas' paper](https://link.aps.org/doi/10.1103/PhysRevX.12.021066).

## Physicality

Description of what a *physical* solutions means.

Assume we have a solution provided by our code, which consist on a longitudinal part $n^L$ 
and a trasverse+longitudinal one $n^{T+L}$. This is obtained by solving [Problem 1](#problem-1). 
Given a symmetry vector $m$ from MPB, we are considering if our solution properly represents 
the content at $\Gamma$, which has been in the computation of the solution. For this we 
need to check two things:

1. Whether our solution subduces properly to the $O(3)$ representation at $\Gamma$ and zero 
   frequency. This can be checked easily using `PhotonicBandConnectivity.jl`. As explained 
   before in [Problem 2](#problem-2), this is fulfilled if there exists a $\mathbf{p}\in\mathbb{Z}$ 
   solving $A_{\Gamma} c^T - m_{\Gamma} = Q\mathbf{p}$.
2. Whether our solution doesn't make use of the higher frequency irreps present in 
   $m_\Gamma^{>0}$ to regularize the symmetry content at zero frequency, and that instead 
   those negative multiplicities in the irreps are cancelled out by the longitudinal modes $n^L$. 
   We ensure this by the following check:

    Define the candidate-solution's zero-frequency content at $\Gamma$ by:

    $$n_\Gamma^{T,=0} = n_{\Gamma}^{T} - m_{\Gamma}^{>0} = n_{\Gamma}^{T+L} - n_{\Gamma}^L 
    - m_{\Gamma}^{>0} = m_{\Gamma}^{=0} + Q\mathbf{p}.$$
  
    Consider the following two cases:
    - If $n_{\Gamma,i}^{T,=0} < 0$ for some $i$, then $n_{\Gamma,i}^L \geq |n_{\Gamma,i}^{T,=0}|$ 
        for that $i$; equivalently, in this case $n_{\Gamma,i}^L \geq -n_{\Gamma,i}^{T,=0}$.
    - Conversely, if  $n_{\Gamma,i}^{T,=0} ≥ 0$ for some $i$, we still have $n_{\Gamma,i}^L ≥ 0$
         and consequently also $n_{\Gamma,i}^L ≥ -n_{\Gamma,i}^{T,=0}$.

    Thus, regardless of the sign of $n_{\Gamma,i}^{T,=0}$, we may require that:

    $$ n_{\Gamma}^L \geq -n_\Gamma^{T,=0}$$

These constraints are directly imposed in the function `find_apolar_modes` thanks to the
functionalities of `find_all_admissible_expansions`.

## Representation of the SG operations in $\mathbf{k}$-space

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

$$ a_{i\alpha}(\mathbf{k},\mathbf{r}) = \frac{1}{\sqrt{N}} \sum_{\mathbf{t}} e^{i\mathbf{k} 
    \cdot (\mathbf{t}+\mathbf{q}_\alpha)} a_{i\alpha}(\mathbf{r}-\mathbf{t}) $$

Then:

$$ \rho_G(g) a_{i\alpha}(\mathbf{k},\mathbf{r}) = \frac{1}{\sqrt{N}} \sum_{\mathbf{t}} 
    e^{i\mathbf{k} \cdot (\mathbf{t}+\mathbf{q}_\alpha)} \rho_G(g) 
    a_{i\alpha}(\mathbf{r}-\mathbf{t}) = \frac{1}{\sqrt{N}} \sum_{\mathbf{t}} e^{i\mathbf{k}
    \cdot (\mathbf{t}+\mathbf{q}_\alpha)} [\rho(h)]_{ji} a_{j\beta}(\mathbf{r}-R\mathbf{t}-
    \mathbf{t}_{\beta\alpha}) $$

If we now make the substitution: $\mathbf{t}^\prime = R\mathbf{t} + \mathbf{t}_{
\beta\alpha} \Rightarrow \mathbf{t} = R^{-1}\mathbf{t}^\prime - \mathbf{q}_\alpha + R^{-1}
\mathbf{q}_\beta - R^{-1}\mathbf{v}$, then:

$$ \rho_G(g) a_{i\alpha}(\mathbf{k},\mathbf{r})
    = \frac{1}{\sqrt{N}} \sum_{\mathbf{t}^\prime} e^{i\mathbf{k} \cdot 
    R^{-1}(\mathbf{t}^\prime+\mathbf{q}_\beta - \mathbf{v})} [\rho(h)]_{ji} 
    a_{j\beta}(\mathbf{r}-\mathbf{t}^\prime) = e^{-i (R\mathbf{k})\cdot\mathbf{v}} 
    [\rho(h)]_{ji} a_{j\beta}(R\mathbf{k},\mathbf{r}) $$

This functionality is implemented under the function `sgrep_induced_by_siteir_generators`.
Take into consideration that the global phase $e^{-i (R\mathbf{k})\cdot\mathbf{v}}$ is 
ignored in `sgrep_induced_by_siteir_generators`.

## Hamiltonian formalism with translational symmetries

Suppose we have an imaginary Bravais lattice with a set of translations denoted by $\mathcal{T}$.
The translations inside such set will be denoted in capital letters $\mathbf{R}$. Additionally,
we may consider that we have different atomic positions inside the unit cell. Those positions
will be denoted by Greek letter $\alpha$ and their positions in the unit cell will be 
$\mathbf{q}_\alpha$. At first, we aren't going to differentiate between atomic positions 
belonging to different WPs or any integral degrees of freedom. Then the most general
Hamiltonian of such Bravais lattice can be written as:

$$ \mathcal{H} = \sum_{\mathbf{r}\mathbf{R},\alpha\beta} t_{\alpha\beta,\mathbf{R}} 
    a_{\alpha,\mathbf{r}}^\dagger a_{\beta,\mathbf{r-R}} + \text{c.c.} $$

We can consider on applying a Fourier transformation to this Hamiltonian, considering the 
translational invariance of it. Let me consider the following Fourier transform of the fields:

$$ a_{\alpha,\mathbf{k}} = \frac{1}{\sqrt{N}} \sum_{\mathbf{t}} e^{i\mathbf{k}\cdot
    (\mathbf{t}+\mathbf{q}_\alpha)} a_{\alpha,\mathbf{r-t}} $$

Then the Hamiltonian in $k$-space will look like:

$$ \mathcal{H} = \sum_{\mathbf{kR},\alpha\beta} t_{\alpha\beta,\mathbf{R}} e^{i\mathbf{k}\cdot
    (\mathbf{q}_\alpha-\mathbf{q}_\beta-\mathbf{R})} a_{\alpha,\mathbf{k}}^\dagger 
    a_{\beta,\mathbf{k}} $$
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

$$ \mathcal{I} a_n = a_{-n}; \quad \mathcal{I} b_n = -b_{-n-1} $$ 

Then the most general inversion-symmetric Hamiltonian for 1-st nearest neighbors will be:

$$ \mathcal{H} = \sum_n t (a_n^\dagger b_n - a_n^\dagger b_{n-1}) + \text{c.c.} $$

It is easy to check that this Hamiltonian is inversion symmetric.

If we consider the following Fourier transform: $a_n = \frac{1}{\sqrt{N}} \sum_n e^{-ikn} a_k$,
and $b_n = \frac{1}{\sqrt{N}} \sum_n e^{-ik(n+1/2)} b_k$
then the Hamiltonian in $k$-space will look like:

$$ \mathcal{H} = \sum_k 2it \sin(k/2) a_k^\dagger b_k + \text{c.c} $$

### Deduction from our method

First, let me consider the translation $t=0$. Then, remember that $\Delta_{\alpha\to\beta+R} = 
\mathbf{q}_\beta + \mathbf{R} - \mathbf{q}_\alpha$.

- Possible $\Delta$'s: $\Delta_{a\to a} = \Delta_{b\to b} = 0; \quad \Delta_{a\to b} = 
    -\Delta_{b\to a} = -1/2$.
- Orbits of $\Delta$'s: $\{\Delta_{a\to a}\} = \{\Delta_{b\to b}\} = \{0\}; \quad 
    \{\Delta_{a\to b}\} = \{\Delta_{b\to a}\} = \{1/2,-1/2\}$.

<font color=red>**NOTE:**</font> be careful that we need to differentiate between clases of 
hopping not only the class itself.

Then the non-symmetrized Hamiltonian for this translation will be:

$$ H_{k,0} = \begin{pmatrix} t_{a\to a} & t_{b\to a}^1e^{ik/2}+t_{b\to a}^2e^{-ik/2} \\ 
    t_{a\to b}^1e^{ik/2}+t_{a\to b}^2e^{-ik/2} & t_{b\to b} \end{pmatrix} $$

Let's proceed now to symmetrize this Hamiltonian. Since we only have inversion we only need 
to check that: $H_{k,0} = \rho(\mathcal{I})H_{k,0}\rho^{-1}(\mathcal{I})$, where:

$$\rho(\mathcal{I}) = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}$$

Then, that constraint impose that: $\left\{ \begin{matrix} t_{b\to a}^1 = -t_{b\to a}^2 \\
t_{a\to b}^1 = -t_{a\to b}^2 \end{matrix} \right.$, so the Hamiltonian will look like:

$$ H_{k,0} = \begin{pmatrix} t_{a\to a} & 2it_{b\to a} \sin(k/2) \\ 
    2it_{a\to b} \sin(k/2) & t_{b\to b} \end{pmatrix} $$

Which is exactly the Hamiltonian deduced with the previous method but with onsite terms and 
without hermiticity imposed.

### Implementation of the free-parameter Hamiltonian into Julia code

For a non-symbolic programming language as Julia we cannot implement this Hamiltonian directly
with some free parameters. In order to do so, we will rewrite it in a clever way so we can 
perform linear algebra with it.

Let me start with our previous example. For this let me focus on only one hopping term. In 
this case it will be $b \to a$. The Hamiltonian for this term will look like:

$$ H_{k,0}^{b \to a} = t_{b\to a}^1e^{ik/2}+t_{b\to a}^2e^{-ik/2}  $$
<font color=red>**NOTE:**</font> I know we discussed to take also the complex conjugate of it
but I think it is not necessary. Let me explained it in the following.

This term can be obtained by just knowing the class $\{\Delta_{a\to b}\}$. Knowing this we 
can define two vectors such that:

$$ \mathbf{v}^T = \begin{pmatrix}
    e^{ik/2} & e^{-ik/2} \\
\end{pmatrix}; \quad \mathbf{t} = \begin{pmatrix}
    t_{b\to a}^1 \\
    t_{b\to a}^2
\end{pmatrix} $$

Then $H_{k,0}^{b \to a}$ can be written as:

$$ H_{k,0}^{b \to a} = \mathbf{v}^T \underbrace{\begin{pmatrix}
    1 & 0 \\
    0 & 1
\end{pmatrix}}_{M^{b \to a}_0} \mathbf{t} $$

In general, I think, this can be performed by creating an identity matrix of dimension 
$\#\{\Delta_{b\to a}\} \times (\mu_a + \mu_b)$.

<font color=red> **NOTE:** </font> I think this dimensions will depend on the internal degrees
of freedom of the orbitals '$\mu$' under study in the way I have specified. Maybe wrong.

How a symmetry transformation will look like if I apply the symmetries to $M$ instead of $H$.
Then they will look as:

$$ \rho_G(g)_{im} H^{mn}_{k} \rho_G^{-1}(g)_{mj} = v_\alpha \rho_G(g)_{im} M^{mn}_{\alpha\beta} 
    \rho_G^{-1}(g)_{nj} t_\beta $$

and

$$ H_{gk}^{ij} = permuted(\mathbf{v}^T) M^{ij} \mathbf{t} = \mathbf{v}^T permuted(M^{ij}) 
    \mathbf{t} $$

#### How do we obtain such permutation?

Since we have $\{\Delta_{b \to a}\}$, we can just do $R\{\Delta_{b\to a}\}$ with $g=\{R|
\mathbf{v}\}$ and find the permutations. Those permutations will be exactly the ones we are 
looking for.