# Time-reversal symmetry in tight-binding models

This document follows — mainly — Chapter 12 of
[Wooten's book](https://www.cambridge.org/core/books/symmetry-and-condensed-matter-physics/218B3D7B149076E63A618D4584E3379B).
First, we present a general way to introduce a non-unitary
transformation into the formalism that can be identified with TRS. Then, we
specialize to TRS and, particularly, to bosonic TRS:
$\Theta² = +\mathbb{I}$.

## Time-reversal representation: theory of corepresentations

We start by considering the combined action of $\Theta$ with another linear
or non-linear operator $\mathcal{O}$.

!!! note
    Here $\Theta$ could be any anti-unitary operator. For our purposes it will be the TRS operator.

```math
\Theta \mathcal{O} \psi_\mu = \Theta \sum \psi_\nu \Gamma_{\nu\mu}(\mathcal{O})
= \sum (\Theta \psi_\nu) \Gamma^*_{\nu\mu}(\mathcal{O}) = \sum_{\nu\mu} \psi_\lambda
\Sigma_\lambda\nu(\Theta) \Gamma^*_{\nu\mu}(\mathcal{O})
```

This demonstrates that the product of the two operators does not lead to just a 
product of the corresponding matrix representatives, but leads, in addition, to 
a c-conjugation of the matrix representative of $\mathcal{O}$.

### Construction of corepresentations (corep)

We consider a magnetic space group $\mathcal{M}$ which we write as

```math
\mathcal{M} = \mathcal{N} + \mathcal{AN},
```

where $\mathcal{N}$ is a unitary subgroup of index 2 (normal subgroup), and 
$\mathcal{A} \notin \mathcal{N}$ an anti-unitary element of $\mathcal{M}$.

!!! note
    This formalism is quite general and can be applied to all kind of magnetic space groups. However, since we are interested in space groups, $\mathcal{N}$ can be identified as the space group and $\mathcal{A}$ as TRS.

!!! note "Notation"
    We denote elements of $\mathcal{N}$ by $R$, $S$, $T$, etc., and those of $\mathcal{AN}$ by $\mathcal{A}$, $\mathcal{B}$, etc..

We start with applying $\Theta$ to a basis set $\{\psi_\mu\} \equiv \Psi$ 
which engenders an irrep $\Delta$ of $\mathcal{N}$, namely,

```math
R \psi_\mu = \sum_\nu \psi_\nu \Delta_{\nu\mu}(R), \\

R \Psi = \Psi \Delta(R).
```

The effect of $\Theta$ on this basis was shown above, and is summarized by

```math
\Theta R \Psi = \Theta (\Psi \Delta(R)) = (\Theta \Psi) \Delta^*(R)
```

Next, we define the generalized time reversed set $\Phi \equiv \{\phi_\mu\} = 
\{\mathcal{A} \psi_\mu\}$ such that

```math
R \Phi = \Phi\quad {}^\mathcal{A}\Delta(R),
```

but, since $\Theta R = R \Theta$, we have

```math
R (\mathcal{A} \Psi) = (\mathcal{A} \Psi) ^\mathcal{A}\Delta(R) = R S \Theta \Psi 
= S (S⁻¹ R S) \Theta \Psi \\
= S \Theta (S⁻¹ R S) \Psi = (S \Theta \Psi) \Delta^*(S⁻¹ R S) = (\mathcal{A} \Psi)
\Delta^*(S⁻¹ R S). \\
\boxed{^\mathcal{A}\Delta(R) = \Delta^*(S⁻¹ R S),}
```

where $\Delta(S⁻¹ R S)$ is an irrep conjugate to $\Delta(R)$, since $\mathcal{N}$ 
is a normal subgroup.

We now construct the rep $\Gamma$ engendered by the combined basis $\mathbf{F} = [\Psi, \mathcal{A}\Psi]$, namely,

```math
\boxed{R \mathbf{F} = \mathbf{F} \Gamma(R) = [\Psi \mathcal{A} \Psi] 
\begin{pmatrix}
    \Delta(R) & \mathbf{0} \\
    \mathbf{0} & \Delta^*(S⁻¹ R S)
\end{pmatrix}, \qquad \forall R \in \mathcal{N}.}
```

Next we apply an operation $\mathcal{B} = \mathcal{A} T \in \mathcal{AN}$, and obtain

```math
\mathcal{B} \Psi = \mathcal{A} T \Psi = \mathcal{A} \Psi \Delta(T) = (\mathcal{A} 
\Psi) \Delta^*(T) = (\mathcal{A} \Psi) \Delta^*(\mathcal{A}⁻¹ \mathcal{B}), \\

\mathcal{B} (\mathcal{A} \Psi) = \mathcal{B} \mathcal{A} \Psi = \Psi 
\Delta(\mathcal{B} \mathcal{A}), \qquad \mathcal{BA} \in \mathcal{N}.
```

We then find

```math
\boxed{\mathcal{B} \mathbf{F} = \mathbf{F} \Gamma(\mathcal{B}) = [\Psi 
\mathcal{A} \Psi] 
\begin{pmatrix}
    \mathbf{0} & \Delta(\mathcal{BA}) \\
    \Delta^*(\mathcal{A⁻¹B}) & \mathbf{0}
\end{pmatrix}, \qquad \forall \mathcal{B} \in \mathcal{AN}.}
```

!!! danger
    The matrix representatives $\Gamma$ do not obey the ordinary multiplication 
    relations associated with unitary groups, but rather obey

    ```math
    \Gamma(R) \Gamma(S) = \Gamma(RS), \qquad \Gamma(R) \Gamma(\mathcal{B}) = 
    \Gamma(R\mathcal{B}), \\
    \Gamma(\mathcal{B}) \Gamma^*(R) = \Gamma(\mathcal{B}R), \qquad \Gamma(\mathcal{B})
    \Gamma^*(\mathcal{C}) = \Gamma(\mathcal{BC}).
    ```

The set of unitary matrices obtained forms a *corepresentation* (corep) of 
$\mathcal{M}$, derived from the unitary irrep $\Delta$ of its normal subgroup 
$\mathcal{N}$.

#### Specialization to grey groups

We now specialize to grey groups (type II), the case relevant to us.
Here $\mathcal{A} = \Theta$ and $\mathcal{N} = \mathcal{G}$, where $\mathcal{G}$
is the space group, so that

```math
\mathcal{M} = \mathcal{G} \oplus \Theta \mathcal{G} = \mathcal{G} \otimes \{E,
\Theta\}.
```

Constructing the rep $\Gamma$ for a transformation $R \in \mathcal{G}$, we
obtain

```math
\boxed{R \mathbf{F} = \mathbf{F} \Gamma(R) = [\Psi \quad \Theta \Psi] 
\begin{pmatrix}
    \Delta(R) & \mathbf{0} \\
    \mathbf{0} & \Delta^*(R)
\end{pmatrix}.}
```

!!! note
    The time-reversed representation $^\mathcal{A}\Delta$ is **identical** to the 
    complex conjugate representation $\Delta^*$.

!!! warning
    What is the difference here with the statement at the beginning? Why are we not
    able to impose the conditions where $R$ and $\Theta$ commute?

    The distinction is due to the fact that $\Psi \not = \Phi = \Theta \Psi$. If
    that were the case, the previous statement could be applied.

Next we do it for any element $\mathcal{B} = \Theta T \in \Theta\mathcal{G}$, and 
obtain

```math
\boxed{\mathcal{B} \mathbf{F} = \mathbf{F} \Gamma(\mathcal{B}) = [\Psi \Theta\Psi] 
\begin{pmatrix}
    \mathbf{0} & \Delta(T\Theta²) \\
    \Delta^*(T) & \mathbf{0}
\end{pmatrix} = [\Psi \Theta\Psi] 
\begin{pmatrix}
    \mathbf{0} & \Delta(T) \\
    \Delta^*(T) & \mathbf{0}
\end{pmatrix}}
```

!!! danger
    Notice that we are using the bosonic TRS, i.e., $\Theta^2 = +\mathbb{I}$.
    To generalize this to fermionic TRS, a minus sign is needed.

!!! note
    In the implementation, we only consider $\Theta$ to include new constraints.
    This can be justified by the fact that grey groups can be decomposed as $\mathcal{M} = \mathcal{G} \otimes \{E, \Theta\}$.

#### Three scenarios for the co-representation

Given the co-representation $\Gamma$ constructed above, three scenarios arise
depending on the relationship between $\Psi$ and $\Theta\Psi$:

1. **Real representation:** $\Theta\Psi$ reproduces the same set as $\Psi$.
   The co-representation coincides with the original representation $\Delta$
   and no new degeneracy is introduced.
2. **Pseudo-real representation:** $\Theta\Psi \neq \Psi$ but $\Theta\Psi$
   also forms a basis for $\Delta$. The co-representation corresponds to
   $\Delta$ with doubled dimension. This case does not occur for
   site-symmetry groups (it requires even-dimensional representations).
3. **Complex representation:** $\Theta\Psi$ forms a basis for a
   representation $\Delta'$ inequivalent to $\Delta$. The co-representation
   pairs $\Delta$ and $\Delta'$, forcing them to become degenerate.

To work uniformly in scenario 1 (the most convenient case), we realify
the representations beforehand using `realify` in Crystalline.jl, which
constructs explicit co-representations in cases 2 and 3 so that the
combined basis transforms under a single real representation.

## Quantization of TRS action on creation and annihilation operators

Assuming a physically real basis (achieved via `physically_realify` in
Crystalline.jl), $\Theta$ acts trivially on the real-space orbitals and
therefore simply conjugates the Bloch phases. The action on creation
operators follows directly:

```math
\hat{\Theta} \hat{a}^\dagger_{I,\mathbf{k}} \hat{\Theta}^{-1} =
\hat{a}^\dagger_{I,-\mathbf{k}}
```

For annihilation operators, we cannot use $\hat{\Theta}^\dagger =
\hat{\Theta}^{-1}$ since $\Theta$ is anti-unitary. Instead, we use
$\hat{\Theta}^2 = +1 \Rightarrow \hat{\Theta} = \hat{\Theta}^{-1}$ and act on a
general single-particle state:

```math
\hat{\Theta} \hat{a}_{I,\mathbf{k}} \hat{\Theta}^{-1}
\ket{\varphi_{J,\mathbf{k}'}} =
\hat{\Theta} \hat{a}_{I,\mathbf{k}} \ket{\varphi_{J,-\mathbf{k}'}} =
\hat{\Theta} \left( \delta_{\mathbf{k},-\mathbf{k}'} \delta_{IJ} \ket{0}
\right) = \delta_{\mathbf{k},-\mathbf{k}'} \delta_{IJ} \ket{0} =
\hat{a}_{I,-\mathbf{k}} \ket{\varphi_{J,\mathbf{k}'}} \\
\Rightarrow \boxed{\hat{\Theta} \hat{a}_{I,\mathbf{k}} \hat{\Theta}^{-1} =
\hat{a}_{I,-\mathbf{k}}}
```

## Finding an explicitly real form of irrep matrices

An explicit real, or physically real, form of a set of irrep matrices is one 
where the associated matrices $D(g)$ have the following property:

```math
D(g) = D^*(g),
```

for all operations $g$ in the considered group $G$.

The standard listings of irreps are not explicitly real. However, if an irrep is 
either intrinsically real - or has been made into a corep in the complex or 
pseudoreal case - it is always equivalent to an intrinsically real form. That is,
there exists a unitary transform $S$ such that:

```math
S D(g) S^{-1} = S D(g) S^\dagger = D^*(g).
```

Suppose we can find this unitary transformation $S$ by some means. What we are 
interested in, is finding a related transform $W$, defining an explicitly real 
form of the irrep:

```math
\tilde{D}(g) = W D(g) W^{-1} = W D(g) W^\dagger,
```

where $W$ is some other unitary transformation and where $\tilde{D}(g)$ is an 
intrinsically real form of $D(g)$, i.e., where

```math
\tilde{D}(g) = \tilde{D}^*(g) \quad \forall g\in G.
```

Our aim is to find $W$, assuming we know $S$. For brevity, we will often write 
$D_g$ in place of $D(g)$.
First, note that $S$ is not merely a unitary matrix: rather, since, by assumption 
$D(g)$ is a "real" matrix, what we really mean is that $S$ is also a _symmetric_ 
unitary matrix, i.e., $S = S^{\mathrm{T}}$ and $S^{-1} = S^\dagger$ (implying, 
jointly, $S^* = S^\dagger = S^{-1}$); this is e.g., derived in Inui p. 74 (bottom) 
to 75 (top). Accordingly $S$ is also normal, i.e., $S S^* = S^* S$.

This property in turn implies that we can express $S$ as the square of another 
symmetric unitary matrix, say $W$, in the sense that $S = W^2$. This follows from 
the following manipulations (Inui, p. 75 bottom), involving the eigendecomposition 
$S = V Λ V^{-1}$ where $\Lambda$ is a diagonal matrix with unit-modulus values 
and $V$ are a set of real eigenvectors (real because $S$ is symmetric unitary) 
and $V^{-1} = V^\dagger = V^{\mathrm{T}}$ (since $S$ is normal).

```math
S = VΛV^{-1} = VΛ^{1/2}Λ^{1/2}V^{\mathrm{T}} = (VΛ^{1/2}V^{\mathrm{T}})
(VΛ^{1/2}V^{\mathrm{T}}),
```

so we can pick $W = VΛ^{1/2}V^{\mathrm{T}}$ (note also that the square root of 
$\Lambda$ must exist and is well-defined since $S$ is invertible, i.e., has full 
rank).
Hence $W^* = V(Λ^{1/2})^*V^{\mathrm{T}} = VΛ^{-1/2}V^{\mathrm{T}} = W^{-1}$ and 
$W^{\mathrm{T}} = (VΛ^{1/2}V^{\mathrm{T}})^{\mathrm{T}} = 
(V^{\mathrm{T}})^{\mathrm{T}}(Λ^{1/2})^{\mathrm{T}} V^{\mathrm{T}} = 
VΛ^{1/2}V^{\mathrm{T}} = W$. I.e., $W$ is also unitary symmetric and normal.

Now, let us rewrite $S D(g) S^{-1} = D^*(g)$ in terms of $W$:

```math
WW D_g W^{-1}W^{-1} = D_g^* \\
```

Multiply from LHS by $W^*$ and from RHS by $W$:

```math
W^*WW D_g W^{-1}W^{-1} W = W^*D_g^* W \\
\Leftrightarrow W D_g W^{-1} = W^*D_g^* W \\
\Leftrightarrow W D_g W^{-1} = W^* D_g^* (W^{-1})^*
```

where we have used the properties of $W$ to reduce the expressions. 

Identifying $\tilde{D}(g) = W D(g) W^{-1}$ we clearly obtain the desired 
invariance under complex conjugation since $\tilde{D}^*(g) = (W D_g W^{-1})^* = 
W^* D_g^* (W^{-1})^* = W D_g W^{-1} = \tilde{D}(g)$.

## Theory of representations in crystalline systems

This section describes how to express a general Hamiltonian in a symmetry-adapted
basis and derives the resulting symmetry constraints. We first study the action
of a symmetry transformation on a basis set in $k$-space, and then derive the
constraints that the Hamiltonian matrix must satisfy.

### Representation of symmetry operators using a basis

Following the deductions made by Bradlyn *et al.* in Ref. [1].

Let us start with a basis set in real space $\{ψ_{iα}
(\mathbf{r})\}$, where $i$ indicates the internal degrees of freedom of 
the orbital, $α$ indicates the site $\mathbf{q}_α$ inside the Wyckoff
position. Notice that by construction we assume each function $ψ_{iα}
(\mathbf{r})$ is localized on $\mathbf{q}_α$. Intuitively, these can be thought of as Wannier functions.

We focus on a particular orbital $ψ_{i1}(\mathbf{r})$
localized in the site $\mathbf{q}_1 \equiv \mathbf{q}$. This orbital will 
transform under the representation $ρ$ of the site-symmetry group $G_\mathbf{q}$,
associated with $\mathbf{q}$. Then, for each $h \in G_\mathbf{q}$:

```math
h ψ_{i1}(\mathbf{r}) = [ρ(h)]_{ji} ψ_{j1}(\mathbf{r})
```

Within the primitive unit cell, an orbital localized on each $\mathbf{q}_α$
can be defined as:

```math
ψ_{iα}(\mathbf{r}) = g_α ψ_{i1}(\mathbf{r}) = ψ_{i1}(g_α^{-1} \mathbf{r}),
```

where $g_α$, with translations, generates the coset decomposition of 
$G_\mathbf{q}$ in $G$. In other words, we can assign for each $\mathbf{q}_α$
a space group element $g \in G$, such that $\mathbf{q}_α = g_α\mathbf{q}$
and:

```math
G = \bigcup_{α=1}^n g_α (G_\mathbf{q} \ltimes T).
```

By extension, translated counterparts in other unit cells can be defined by:

```math
\{E|\mathbf{t}\} ψ_{iα}(\mathbf{r}) = ψ_{iα}(\mathbf{r-t}),
```

where $\mathbf{t}$ is a lattice translation. The set of $n \times \text{dim}(ρ)
\times \mathcal{N}$ functions $ψ_{iα}(\mathbf{r-t})$, where $\mathcal{N}$ is
the number of unit cells of the system, will be the basis set on which the induced 
representation $D$ will act.

Specifically, given $g = \{R|\mathbf{v}\} \in G$, the coset decomposition implies 
that for each $g_α$, there is a unique operation $g_β$ such that:

```math
g g_α = \{E|\mathbf{t}_{αβ}\} g_β h,
```

where $h \in G_\mathbf{q}$ and $t_{αβ} \equiv g\mathbf{q}_α - 
\mathbf{q}_β$.

!!! note
    This follows from the coset decomposition; see, e.g., Bradlyn *et al.* [1].
    Note that this reference does not include an explicit proof of the above equation.

!!! todo
    Consider adding a proof as an appendix in the future.

Taking all of this into consideration, we can deduce how our basis set will 
transform under the action of every $g \in G$:

```math
g ψ_{iα}(\mathbf{r-t}) = g \{E|\mathbf{t}\} ψ_{iα}(\mathbf{r}) = \\
\{E|R\mathbf{t}\} g ψ_{iα}(\mathbf{r}) = \\
\{E|R\mathbf{t}\} \{E|\mathbf{t}_{αβ}\} g_β h ψ_{i1}(\mathbf{r}) = \\
\{E|R\mathbf{t}\} \{E|\mathbf{t}_{αβ}\} g_β [ρ(h)]_{ji} 
ψ_{j1}(\mathbf{r}) = \\
\{E|R\mathbf{t}\} \{E|\mathbf{t}_{αβ}\} [ρ(h)]_{ji} 
ψ_{jβ}(\mathbf{r}) = \\
\{E|R\mathbf{t}\} [ρ(h)]_{ji} ψ_{jβ}(\mathbf{r-\mathbf{t}_{αβ}}) = \\
[ρ(h)]_{ji} ψ_{jβ}(\mathbf{r}-R\mathbf{\mathbf{t}-\mathbf{t}_{αβ}})
```

While it is natural to define the representation in real space, it will be more
useful to view it in reciprocal space. This is more evident when $\mathcal{N} 
\to \infty$. To this end, we define the Fourier transform of our basis:

```math
φ_{iα,\mathbf{k}}(\mathbf{r}) = \sum_\mathbf{t} e^{i\mathbf{k\cdot(t+q_α)}} 
ψ_{iα}(\mathbf{r-t}),
```

where the sum is over all lattice vectors $\mathbf{t} \in T$.

!!! note
    Notice that this is just convention but for building a tight-binding Hamiltonian this choice is better since we will eliminate all local phases.

The Fourier transform amounts to a unitary transformation that exchanges 
$\mathcal{N}$ unit cells for $\mathcal{N}$ distinct $\mathbf{k}$ points. The 
action of $g \in G$ in reciprocal space becomes:

```math
g φ_{iα,\mathbf{k}}(\mathbf{r}) = \sum_\mathbf{t} e^{i\mathbf{k}\cdot(\mathbf{t}+\mathbf{q}_α)} g
ψ_{iα}(\mathbf{r-t}) = \\
\sum_\mathbf{t} e^{i\mathbf{k\cdot(t+q_α)}} [ρ(h)]_{ji} 
ψ_{jβ}(\mathbf{r}-R\mathbf{\mathbf{t}-\mathbf{t}_{αβ}}) = \\
\sum_\mathbf{t}' e^{i\mathbf{k}\cdot R^{-1}(\mathbf{t}'+\mathbf{q}_β-\mathbf{v})}} [ρ(h)]_{ji} 
ψ_{jβ}(\mathbf{r}-\mathbf{t}') = \\
e^{-i([R⁻¹]ᵀ \mathbf{k}) \cdot \mathbf{v}} [ρ(h)]_{ji} \sum_\mathbf{t'} 
e^{i([R⁻¹]ᵀ \mathbf{k}) \cdot (\mathbf{t}'+\mathbf{q}_β)} ψ_{jβ}(\mathbf{r-t'}) = \\
e^{-i([R⁻¹]ᵀ \mathbf{k}) \cdot \mathbf{v}} [ρ(h)]_{ji} φ_{jβ,[R⁻¹]ᵀ\mathbf{k}}(\mathbf{r}),
```

where we have made the substitution: $\mathbf{t}' = R\mathbf{t} + \mathbf{t}_{αβ}
= R\mathbf{t} + g\mathbf{q_α - q_β} = R\mathbf{t} + R\mathbf{q_α + v - q_β} =
R(\mathbf{t+q_α}) + \mathbf{v - q_β} \Rightarrow (\mathbf{t+q_α}) = R^{-1}
(\mathbf{t'+q_β-v})$.

!!! note
    We used the identity $\mathbf{k}·(R⁻¹\mathbf{r}) \equiv (g \mathbf{k})·\mathbf{r}$, which follows from:

    ```math
    \mathbf{k}·(R⁻¹\mathbf{r}) = \sum_{ij} k_i (R⁻¹_{ij} r_j) = \sum_{ij} (R⁻¹_{ij} k_i)
    r_j = ([R⁻¹]ᵀ \mathbf{k}) · \mathbf{r} \equiv (g \mathbf{k}) · \mathbf{r},
    ```

In reciprocal space, the matrix representation can be interpreted as a $\mathcal{N}
\times \mathcal{N}$ matrix of $n\dim(ρ) \times n\dim(ρ)$ blocks, each block can 
be labeled by $\mathbf{k,k'}$. Most of the blocks are zero: given $g = \{R|
\mathbf{v}\} \in G$, there is only one non-zero block in each row and column, 
corresponding to $\mathbf{k'} = R\mathbf{k}$. Mathematically, we can express this
as:

```math
g φ_{iα,\mathbf{k}}(\mathbf{r}) = \sum_{jβ\mathbf{k'}} D_{jβ\mathbf{k'},iα\mathbf{k}}(g)
φ_{jβ,\mathbf{k}'}(\mathbf{r}),
```

where we have that:

```math
D_{jβ\mathbf{k'},iα\mathbf{k}}(g) = e^{-i(g\mathbf{k) \cdot v}} ρ_{ji}(h)
\delta_{g\mathbf{k,k'}} \delta_{g\mathbf{q_α - q_β} \mod \tau},
```

where $\tau \in T$.

We will use the following notation:

```math
Ρ_{jβ,iα}(g) = e^{-i(g\mathbf{k) \cdot v}} ρ_{ji}(h) 
\delta_{g\mathbf{q_α - q_β} \mod \tau},
```

where the dependence on $\mathbf{k}$ is left implicit.

We can vectorize the previous equation as:

```math
\boxed{g Φ_\mathbf{k}(\mathbf{r}) = Ρ^T(g) Φ_{g\mathbf{k}}(\mathbf{r})},
```

where $Φ_\mathbf{k}(\mathbf{r})$ is a column vector formed by 
$\{φ_{iα,\mathbf{k}}(\mathbf{r})\}$, and, $Ρ(g)$ is a $n \times n$ matrix of 
$\dim(ρ) \times \dim(ρ)$ blocks, each of them can be labelled by $α,β$. Most of 
the blocks are zero: given $g \in G$, there is only one non-zero block in each 
row and column, corresponding to $g q_α - q_β = 0 \mod \tau$ with $\tau \in T$, 
and is equal to:

```math
Ρ_{jβ,iα}(g)= e^{-i(g\mathbf{k}) \cdot \mathbf{v}} [ρ(h)]_{ji} 
\delta_{g\mathbf{q}_α - \mathbf{q}_β} \mod \tau}
```

!!! note
    We pick the previous definition of the matrix in order to have good properties of composition.
    This is due to the fact that:

    ```math
    g₁ g₂ Φ_\mathbf{k}(\mathbf{r}) = Ρ^T(g₁g₂) Φ_{g₁g₂\mathbf{k}}(\mathbf{r}) \\
    = g₁ Ρ^T(g₂) Φ_{g₂\mathbf{k}}(\mathbf{r}) = Ρ^T(g₂) Ρ^T(g₁) 
    Φ_{g₁g₂\mathbf{k}}(\mathbf{r}) \\ 
    \Rightarrow \boxed{Ρ(g₁g₂) = Ρ(g₁) Ρ(g₂)}
    ```

### Action of symmetry operators on a Hamiltonian

Let us start with the most general non-interacting Hamiltonian:

```math
\hat{H} = \sum_{IJ,\mathbf{R}\mathbf{R}'} h_{IJ,\mathbf{R}-\mathbf{R}'} \; 
\hat{c}_{I,\mathbf{R}}^\dagger \hat{c}_{J,\mathbf{R}'},
```

where $I,J$ collect the internal degrees of freedom of the orbitals and the sites
of the WP, i.e., $I = (i, α)$; and $\mathbf{R,R}'$ run over the lattice translations.

!!! note
    We have here assumed that hopping terms depend only on relative distances.
    We denote $\mathbf{t} \equiv \mathbf{R}' - \mathbf{R}$.

To be consistent with the Fourier transform convention above, the creation operator transforms as:

```math
ĉ_{I,𝐑}^† = \frac{1}{\sqrt{N}} \sum_{𝐤} e^{-i𝐤·(𝐑+𝐪_α)} â_{I,𝐤}^†,
```

obtaining:

```math
\hat{H} = \frac{1}{N} \sum_{IJ,\mathbf{RR}'} h_{IJ,\mathbf{t}} \sum_{\mathbf{kk}'}
e^{-i\mathbf{k·(R+q_α)}} e^{i\mathbf{k'·(R'+q_β)}} \hat{a}^\dagger_{I,\mathbf{k}} 
\hat{a}_{J,\mathbf{k}'} \\
= \frac{1}{N} \sum_{IJ,\mathbf{t},\mathbf{kk}'} h_{IJ,\mathbf{t}} 
\left[ \sum_{\mathbf{R}'} e^{i\mathbf{(k'-k)·R}} \right] e^{i\mathbf{k·(t-q_α)}} 
e^{i\mathbf{k'·q_β}} \hat{a}^\dagger_{I,\mathbf{k}} \hat{a}_{J,\mathbf{k}'} \\
= \sum_{IJ,\mathbf{t},\mathbf{kk}'} h_{IJ,\mathbf{t}} 
\delta_{\mathbf{k,k'}} e^{i\mathbf{k·(t-q_α)}} e^{i\mathbf{k'·q_β}} 
\hat{a}^\dagger_{I,\mathbf{k}} \hat{a}_{J,\mathbf{k}'} \\
= \sum_{IJ,\mathbf{t},\mathbf{k}} h_{IJ,\mathbf{t}} 
e^{i\mathbf{k·(t+q_β-q_α)}} \hat{a}^\dagger_{I,\mathbf{k}} \hat{a}_{J,\mathbf{k}} \\
= \sum_{IJ,\mathbf{k}} h_{IJ,\mathbf{k}} \hat{a}^\dagger_{I,\mathbf{k}} 
\hat{a}_{J,\mathbf{k}},
```

where we have defined: $h_{IJ,\mathbf{k}} = \sum_\mathbf{t} h_{IJ,\mathbf{t}}
e^{i\mathbf{k}·(\mathbf{t}+\mathbf{q}_β-\mathbf{q}_α)}$.

#### Quantization of the representations

The quantization of the previous (classical) theory of representations can be 
written using "braket" notation as:

```math
\hat{g} \ket{φ_{I,\mathbf{k}}} = Ρ_{JI}(g) \ket{φ_{J,g\mathbf{k}}},
```

where $\ket{φ_{I,\mathbf{k}}} \equiv a^\dagger_{I,\mathbf{k}} \ket{0}$. Then:

```math
\hat{g} \hat{a}^\dagger_{I,\mathbf{k}} \hat{g}^{-1} \ket{0} = \hat{g} 
\hat{a}^\dagger_{I,\mathbf{k}} \ket{0} = Ρ_{JI}(g) 
\hat{a}^\dagger_{J,g\mathbf{k}} \ket{0} \\
\Rightarrow \boxed{\hat{g} \hat{a}^\dagger_{I,\mathbf{k}} \hat{g}^{-1} = Ρ_{JI}(g) 
\hat{a}^\dagger_{J,g\mathbf{k}}}
```

!!! note
    In the above, we used that $\hat{g}^{-1} \ket{0} = \ket{0}$, i.e., symmetries act trivially on the vacuum.

Consequently: 

```math
\hat{g} \hat{a}_{I,\mathbf{k}} \hat{g}^\dagger = Ρ^*_{JI}(g)
\hat{a}_{J,g\mathbf{k}}
\Rightarrow \boxed{\hat{g} \hat{a}_{I,\mathbf{k}} \hat{g}^{-1} = Ρ^*_{JI}(g) 
\hat{a}_{J,g\mathbf{k}}}
```

!!! note
    This step assumes that $\hat{g}$ is unitary, i.e., $\hat{g}^\dagger = \hat{g}^{-1}$.


Then, if we want the Hamiltonian to be invariant under the symmetries, we must 
impose that:

```math
\hat{H} = \hat{g} \hat{H} \hat{g}^{-1}
```

Then we obtain that:

```math
\hat{H} = \sum_{IJ,\mathbf{k}} \hat{a}^\dagger_{I,\mathbf{k}} h_{IJ,\mathbf{k}}
\hat{a}_{J,\mathbf{k}} = \\
\hat{g} \hat{H} \hat{g}^{-1} = \sum_{IJ,\mathbf{k}} \hat{g}
\hat{a}^\dagger_{I,\mathbf{k}} h_{IJ,\mathbf{k}} \hat{a}_{J,\mathbf{k}} \hat{g}^{-1} \\
= \sum_{IJ,\mathbf{k}} \hat{g} \hat{a}^\dagger_{I,\mathbf{k}} \hat{g}^{-1} 
h_{IJ,\mathbf{k}} \hat{g} \hat{a}_{J,\mathbf{k}} \hat{g}^{-1} \\
= \sum_{IJ,\mathbf{k},I'J'} \hat{a}^\dagger_{I',g\mathbf{k}} Ρ_{I'I}(g) h_{IJ,\mathbf{k}}
Ρ^*_{J'J}(g) \hat{a}_{J',g\mathbf{k}} \\
= \sum_{IJ,\mathbf{k}} \hat{a}^\dagger_{I,g\mathbf{k}} \left[ Ρ(g) H_{\mathbf{k}}
Ρ^\dagger(g) \right]_{IJ} \hat{a}_{J,g\mathbf{k}},
```

where we have defined $H_\mathbf{k} \equiv h_{IJ,\mathbf{k}}$ and made the substitution
$I',J' \to I,J$. Comparing the first and final rows we obtain the following 
relation for the Hamiltonian to be invariant under symmetries:

```math
\boxed{H_\mathbf{k} = Ρ(g) H_{g^{-1}\mathbf{k}} Ρ^\dagger(g)}
```

!!! note
    Notice that the representations of spatial operations are unitary, so we end up with:

    ```math
    \boxed{H_\mathbf{k} = Ρ(g) H_{g^{-1}\mathbf{k}} Ρ⁻¹(g)},
    ```

    which is the more common form.

### Time reversal symmetry

For TRS, a similar computation can be performed. Let us assume that the action of
TRS over our basis is the following:

```math
Θ \ket{φ_{I,\mathbf{k}}} = \ket{φ_{I,\mathbf{-k}}},
```

then, we obtain the following relations:

```math
\boxed{\hat{Θ} \hat{a}^\dagger_{I,\mathbf{k}} \hat{Θ}^{-1} = 
\hat{a}^\dagger_{I,-\mathbf{k}}; \quad \hat{Θ} \hat{a}_{I,\mathbf{k}} 
\hat{Θ}^{-1} = \hat{a}_{I,-\mathbf{k}}}
```

Then, the invariance under TRS of the Hamiltonian is simply reduced to:

```math
\hat{Θ} \hat{H} \hat{Θ}^{-1} = \sum_{IJ,\mathbf{k}} \hat{Θ} 
\hat{a}^\dagger_{I,\mathbf{k}} h_{IJ,\mathbf{k}} \hat{a}_{J,\mathbf{k}} \hat{Θ}^{-1} \\
= \sum_{IJ,\mathbf{k}} \hat{Θ} \hat{a}^\dagger_{I,\mathbf{k}} \hat{Θ}^{-1} 
h^*_{IJ,\mathbf{k}} \hat{Θ} \hat{a}_{J,\mathbf{k}} \hat{Θ}^{-1} \\
= \sum_{IJ,\mathbf{k}} \hat{a}^\dagger_{I,-\mathbf{k}} h^*_{IJ,\mathbf{k}} 
\hat{a}_{J,-\mathbf{k}} = \\
\hat{H} = \sum_{IJ,\mathbf{k}} \hat{a}^\dagger_{I,\mathbf{k}} h_{IJ,\mathbf{k}} 
\hat{a}_{J,\mathbf{k}}
```

Obtaining the following relation:

```math
\boxed{H_\mathbf{k} = H^*_{-\mathbf{k}}}
```

## Proof that physically real representations admit a real basis 

In the presence of time-reversal symmetry, i.e., $Θ g = g Θ$, and for an explicitly
real representation $Ρ(g) = Ρ^*(g)$, we now show that this implies the existence of
a real basis for the representation.
First, we recall the definition of the representation as it acts on a basis
element $φ_{I,\mathbf{k}}$: $g φ_{I,\mathbf{k}} = Ρ_{JI}(g) φ_{J,g\mathbf{k}}$.
On the other hand, in the presence of time-reversal symmetry, we also have
$Θ g = g Θ$. Accordingly:

```math
g Θ φ_{I,\mathbf{k}} = Θ (g φ_{I,\mathbf{k}}) = Θ (Ρ_{JI}(g) φ_{J,g\mathbf{k}}) 
= Ρ^*_{JI}(g) (Θ φ_{J,g \mathbf{k}}) = Ρ_{JI}(g) (Θ φ_{J,g\mathbf{k}}).
```

Then $Θ φ_{I,\mathbf{k}}$ is another basis element for the _same_ (i.e., identical)
representation $Ρ_{JI}(g)$ as $φ_{I,\mathbf{k}}$. As a result, they can be chosen such
that $φ_{I,\mathbf{k}} = Θ φ_{I,\mathbf{k}}$.

!!! note
    The previous argument holds for irreducible representations. Specifically, since
    both $\{φ_{I,𝐤}\}$ and $\{Θ φ_{I,𝐤}\}$ are orthonormal bases for the same irrep
    space $V$, there exists a unitary mapping $U : V \to V$ with $Uφ_{I,𝐤} = Θφ_{I,𝐤}$.
    To apply Schur's lemma, we need to show that $U$ commutes with the group elements $g$.
    Therefore, we act with $g$ on both sides of $Uφ_{I,𝐤} = Θφ_{I,𝐤}$:

    ```math
    gUφ_{I,𝐤} = gΘφ_{I,𝐤} = P_{JI}(g)(Θφ_{J,g𝐤}) = P_{JI}(g)(Uφ_{J,g𝐤}) = U(gφ_{I,𝐤}),
    ```

    so $[U, g] = 0$ for all $g$. By Schur's lemma (applied to the linear operator 
    $U$ on an irreducible space), $U = e^{iα}𝟏$, giving:

    ```math
    Θφ_{I,𝐤} = e^{iα}φ_{I,𝐤}.
    ```

    I.e., $φ_{I,𝐤}$ and $Θφ_{I,𝐤}$ differ only by a phase. This phase is removable, since
    we can introduce a new, real basis element by defining $\tilde{φ}_{I,𝐤} = e^{iα/2} φ_{I,𝐤}$.
    Then, using the antilinearity of $Θ$,

    ```math
    Θ \tilde{φ}_{I,𝐤} = e^{-iα/2} Θφ_{I,𝐤} = e^{-iα/2} e^{iα} φ_{I,𝐤} = e^{iα/2} φ_{I,𝐤} 
    = \tilde{φ}_{I,𝐤},
    ```

    showing that $\tilde{φ}_{I,𝐤}$ is a real function.
    Although this argument assumed irreducible representations (via its invocation of Schur's
    lemma), it can be extended to general representations by simply considering them brought
    to block-diagonal form, where the argument then applies to each block individually.

## Appendix A

For completeness, we derive the commutators of a symmetry $g ∈ G$ with
the creation and annihilation operators, as an alternative to the conjugation
relations used above.

```math
[\hat{g}, \hat{a}^\dagger_{I,\mathbf{k}}] \ket{0} = \hat{g} 
\hat{a}^\dagger_{I,\mathbf{k}} \ket{0} - \hat{a}^\dagger_{I,\mathbf{k}} \hat{g}
\ket{0} = \hat{g} \ket{φ_{I,\mathbf{k}}} - \hat{a}^\dagger_{I,\mathbf{k}} \ket{0} \\
= Ρ_{JI}(g) \ket{φ_{J,g\mathbf{k}}} - \hat{a}^\dagger_{I,\mathbf{k}} \ket{0} = 
\left( Ρ_{JI}(g) \hat{a}^\dagger_{J,g\mathbf{k}} - \hat{a}^\dagger_{I,\mathbf{k}} 
\right) \ket{0} \\
\Rightarrow \boxed{[\hat{g}, \hat{a}^\dagger_{I,\mathbf{k}}] = 
Ρ_{JI}(g) \hat{a}^\dagger_{J,g\mathbf{k}} - \hat{a}^\dagger_{I,\mathbf{k}}}
```

Now we want to do a similar computation for the annihilation operator. However, 
since $[a,a^\dagger] ≠ 0$, we cannot do the previous trick. We will use a more 
general single-particle state $\ket{φ_{I',\mathbf{k}'}}$:

```math
\left[\hat{g}, \hat{a}_{I,\mathbf{k}} \right] \ket{φ_{I',\mathbf{k}'}} = \hat{g} 
\hat{a}_{I,\mathbf{k}} \ket{φ_{I',\mathbf{k}'}} - \hat{a}_{I,\mathbf{k}} \hat{g} 
\ket{φ_{I',\mathbf{k}'}} = δ_{II'} δ_{\mathbf{kk}'} \hat{g} \ket{0} - P_{JI'}(g)
\hat{a}_{I,\mathbf{k}} \ket{φ_{J,\mathbf{gk}}} \\
\Rightarrow \boxed{\left[\hat{g}, \hat{a}_{I,\mathbf{k}} \right] = 0}.
```

Then:

```math
[\hat{g}, \hat{H}] = \sum_{IJ,\mathbf{k}} h_{IJ,\mathbf{k}} [\hat{g}, 
\hat{a}^\dagger_{I,\mathbf{k}} \hat{a}_{J,\mathbf{k}}] = \sum_{IJ,\mathbf{k}} 
h_{IJ,\mathbf{k}} \left( [\hat{g}, \hat{a}^\dagger_{I,\mathbf{k}}] 
\hat{a}_{J,\mathbf{k}} + \hat{a}^\dagger_{I,\mathbf{k}} [\hat{g}, 
\hat{a}_{J,\mathbf{k}}] \right) \\
= \sum_{IJ,\mathbf{k}, I'} h_{IJ,\mathbf{k}} \left[ Ρ_{I'I}(g) 
\hat{a}^\dagger_{I',g\mathbf{k}} - \hat{a}^\dagger_{I,\mathbf{k}} \right] 
\hat{a}_{J,\mathbf{k}}
```

!!! warning
    The resulting expression is not straightforwardly useful; one would need to
    relate $\hat{a}_{J,\mathbf{k}}$ to $\hat{a}_{J,g\mathbf{k}}$ to simplify further.

## References

[1] Band Representations and Topological Quantum Chemistry by Bradlyn *et al.*
(2021) https://doi.org/10.1146/annurev-conmatphys-041720-124134 