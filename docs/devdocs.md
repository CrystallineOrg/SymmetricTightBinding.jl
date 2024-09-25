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
    $$n_\Gamma^{T,=0} = n_{\Gamma}^{T} - m_{\Gamma}^{>0} = n_{\Gamma}^{T+L} - n_{\Gamma}^L - m_{\Gamma}^{>0} = m_{\Gamma}^{=0} + Q\mathbf{p}.$$
    Consider the following two cases:
    - If $n_{\Gamma,i}^{T,=0} < 0$ for some $i$, then $n_{\Gamma,i}^L \geq |n_{\Gamma,i}^{T,=0}|$ for
    that $i$; equivalently, in this case $n_{\Gamma,i}^L \geq -n_{\Gamma,i}^{T,=0}$.
    - Conversely, if  $n_{\Gamma,i}^{T,=0} ≥ 0$ for some $i$, we still have $n_{\Gamma,i}^L ≥ 0$ and consequently also $n_{\Gamma,i}^L ≥ -n_{\Gamma,i}^{T,=0}$.

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

being $g \in G$ and $h \in G_{\mathbf{q}}$, such that given a coset decomposition $G = \bigcup_\alpha
g_\alpha G_\mathbf{q}$ : $g = \{E|\mathbf{t}_{\beta\alpha}\} g_\beta h g_\alpha^{-1}$ and
$\mathbf{t}_{\beta\alpha} = g\mathbf{q}_\alpha - \mathbf{q}_\beta$. We assume we know a 
site-symmetry group representation: $h a_{i1}(\mathbf{r}) = [\rho(h)]_{ji} a_{j1}(\mathbf{r})$,
where we are assuming Einstein's summation rule.

Considering the following Fourier transform we want to study how Bloch functions will transform:

$$ a_{i\alpha}(\mathbf{k},\mathbf{r}) = \sum_{\mathbf{t}} e^{i\mathbf{k} \cdot 
    (\mathbf{t}+\mathbf{q}_\alpha)} a_{i\alpha}(\mathbf{r}-\mathbf{t}) $$

Then:

$$ \rho_G(g) a_{i\alpha}(\mathbf{k},\mathbf{r}) = \sum_{\mathbf{t}} e^{i\mathbf{k} \cdot 
    (\mathbf{t}+\mathbf{q}_\alpha)} \rho_G(g) a_{i\alpha}(\mathbf{r}-\mathbf{t}) = 
    \sum_{\mathbf{t}} e^{i\mathbf{k} \cdot 
    (\mathbf{t}+\mathbf{q}_\alpha)} [\rho(h)]_{ji} a_{j\beta}(\mathbf{r}-R\mathbf{t}-
    \mathbf{t}_{\beta\alpha}) $$

If we now make the substitution: $\mathbf{t}^\prime = R\mathbf{t} + \mathbf{t}_{
\beta\alpha} \Rightarrow \mathbf{t} = R^{-1}\mathbf{t}^\prime - \mathbf{q}_\alpha + R^{-1}
\mathbf{q}_\beta - R^{-1}\mathbf{v}$, then:

$$ \rho_G(g) a_{i\alpha}(\mathbf{k},\mathbf{r})
 = \sum_{\mathbf{t}^\prime} 
 e^{i\mathbf{k} \cdot R^{-1}(\mathbf{t}^\prime+\mathbf{q}_\beta - \mathbf{v})} [\rho(h)]_{ji} 
    a_{j\beta}(\mathbf{r}-\mathbf{t}^\prime)
= e^{-i (R\mathbf{k})\cdot\mathbf{v}} [\rho(h)]_{ji} a_{j\beta}(R\mathbf{k},\mathbf{r}) $$

This functionality is implemented under the function `obtain_sg_representation`.