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

Additionally it is interesting to separate those symmetry vectors into the irreps belonging 
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

Basically the problem we want to solve can be written as: $An^T=m$, where $A$ is the matrix 
of BRs in the SG. In the code $A$ will be called `brs`. The notation explained for the 
symmetry vectors can also be applied to this matrix. We can divide this problem into two 
enabling us to study the $\Gamma$-point separately:

$$ \left[ \begin{matrix} A_\Gamma \\ A_{-\Gamma} \end{matrix} \right] n^T = \left[ 
    \begin{matrix} m_\Gamma \\ m_{-\Gamma} \end{matrix} \right] $$

### Problem 1

$$ A_{-\Gamma} n^T = m_{-\Gamma} \Rightarrow A_{-\Gamma} n^{T+L} = m_{-\Gamma} + A_{-\Gamma}
    n^L $$

In general $n^T$ could have some irreps with negative multiplicities, that why we separated 
the problem into $n^{T+L}$ and $n^L$ so the problem will be strictly positive. For more details 
check [Antonio's paper](https://doi.org/10.48550/arXiv.2305.18257).

First we find all possible longitudinal modes $n^L$ in the SG for a given dimension `t`. 
This is performed with the function `find_auxiliary_modes`.

Then, using such $n^Ls$ founded, we search if there are transverse+longitudinal modes $n^{T+L}$ 
that satisfies the previous equation. This is done inside the function 
`find_all_band_representations`.

### Problem 2

$$ A_{\Gamma} n^T = m_{\Gamma} = m_\Gamma^{=0} + m_\Gamma^{>0} $$

`MPBUtils.jl` forced the content at zero frequency to me exactly at the ones shown in Table 
S(6-8) in [Thomas' paper](https://link.aps.org/doi/10.1103/PhysRevX.12.021066) which we are 
going to called `n_fixed`, but we have some freedom coming from $Q\mathbf{p}$ which can appear 
in our model. Then, the previous equality should be thought in term of $\mod Q$. If we write 
it without such ambivalence:

$$ A_{\Gamma} n^T - m_{\Gamma} = Q\mathbf{p} $$

where $\mathbf{p}\in\mathbb{Z}$.

See the details in [Thomas' paper](https://link.aps.org/doi/10.1103/PhysRevX.12.021066).

This problem is solved in the function `physical`.

## Physicality

Description on how the function `physical` works.

Assume we have a solution provided by our code, which consist on a longitudinal part $n^L$ 
and a trasverse+longitudinal one $n^{T+L}$. This is obtain solving [Problem 1](#problem-1). 
Given a symmetry vector from MPB $m$ we are considering if our solution represents properly 
the the content at $\Gamma$, which has been in the computation of the solution. For this we 
need to check two things:

1. Whether if our solution subduce properly the $O(3)$ representation at $\Gamma$ and zero 
   frequency. This can be check easily using `PhotonicBandConnectivity.jl`. As explained 
   before in [Problem 2](#problem-2), this is fulfilled if $\mathbf{p}\in\mathbb{Z}$.
2. Whether if our solution doesn't make use of the higher frequency irreps present in 
   $m_\Gamma^{>0}$ to regularize the symmetry content at zero frequency, and that instead 
   those negative multiplicities in the irreps are cancel out by the longitudinal modes $n^L$. 
   This is checked in the following way:

    - If $n_{\Gamma,i}^{T+L} > m^{>0}_\Gamma \quad \forall i$, all higher frequency irreps 
        are present on the model so the model perfectly represent the system bans.
    - If $n_{\Gamma,i}^{T+L} < m^{>0}_\Gamma$ for some $i$ that higher frequency irrep is 
        lacking in the model so the solution do not proper represent the higher frequency 
        bands so it must be discarded.

    Those conditions can be easily checked by the following condition:

    $$ n_{\Gamma,i}^{T+L} - m^{>0}_\Gamma = abs(n_{\Gamma,i}^{T+L} - m^{>0}_\Gamma) $$

Both of this conditions are checked for every solution on the function `physical`.