# Notes on How to Implement TRS in TB Models

I am going to follow - mainly - Chapter 12 of 
[Wooten's book](https://www.cambridge.org/core/books/symmetry-and-condensed-matter-physics/218B3D7B149076E63A618D4584E3379B). 
Firstly, I am going to present a general way to introduce a non-unitary 
transformation into the formalism which could be thought as TRS. Secondly, I am 
going to particularize it into TRS symmetry and, particularly, to bosonic TRS: 
$\Theta² = +\mathbb{I}$.

## How TRS Acts in the Presence of Another Spatial Symmetry?

One interesting question to consider and the key point of all of this discussion 
is to understand how TRS interacts with other spatial symmetries. In general, 
and in particular for this section, TRS operator $\Theta$ will commute with 
any element $R$ of a point group or a space group, i.e.,

$$
R \Theta = \Theta R 
$$

and therefore, have *real representations*. This is really important for our 
case and it is something I didn't fully understand yet. Here is a sketch of the 
proof:

---

We right the system wave function as 

$$
\psi(\mathbf{r}) = \psi_\mathcal{R}(\mathbf{r}) + i \psi_\mathcal{I}(\mathbf{r})
$$

where $\psi_\mathcal{R}(\mathbf{r})$ and $\psi_\mathcal{I}(\mathbf{r})$ are 
real functions. Operating with $\Theta R$ and later with $R \Theta$, we 
obtain

$$
\Theta R \psi(\mathbf{r}) = \Theta \psi(R⁻¹\mathbf{r}) = 
(\psi_\mathcal{R}(R⁻¹\mathbf{r}) + i \psi_\mathcal{I}(R⁻¹\mathbf{r}))^* \\

R \Theta \psi(\mathbf{r}) = \Theta \psi(R⁻¹\mathbf{r}) = 
\psi_\mathcal{R}(R⁻¹\mathbf{r}) - i \psi_\mathcal{I}(R⁻¹\mathbf{r})
$$

The two results are equal since $R$ is orthogonal. This, somehow, explains 
that the representations are also real...

---

> [!NOTE] 
> Only valid for spinless systems. I do not understand fully the implications of 
> this...

However, this is a really simple case for state kets. What happens when we 
consider the action of TRS on a more abstract representation basis set?

## Time-Reversed Representation: Theory of Corepresentations

Let me start considering how the combined action of $\Theta$ with other linear 
or non-linear operator $\mathcal{O}$.

> [!NOTE] 
> Here $\Theta$ could be any anti-unitary operator. For our purpose it will be 
> just TRS operator

$$
\Theta \mathcal{O} \psi_\mu = \Theta \sum \psi_\nu \Gamma_{\nu\mu}(\mathcal{O})
= \sum (\Theta \psi_\nu) \Gamma^*_{\nu\mu}(\mathcal{O}) = \sum_{\nu\mu} \psi_\lambda
\Sigma_\lambda\nu(\Theta) \Gamma^*_{\nu\mu}(\mathcal{O})
$$

This demonstrate that the product of the two operators does not lead to just a 
product of the corresponding matrix representatives, but leads, in addition, to 
a c-conjugation of the matrix representative of $\mathcal{O}$.

### Construction of Corepresentations (CoRep)

We consider a magnetic space group $\mathcal{M}$ which we write as

$$
\mathcal{M} = \mathcal{N} + \mathcal{AN},
$$
where $\mathcal{N}$ is a unitary subgroup of index 2 (normal subgroup), and 
$\mathcal{A} \notin \mathcal{N}$ an anti-unitary element of $\mathcal{M}$.

> [!NOTE] 
> For our purpose $\mathcal{N}$ is the space group and $\mathcal{A}$ is TRS.
> This is because we are interested in grey groups (for now).

> [!IMPORTANT] 
> **Notation:** We denote elements of $\mathcal{N}$ by $R$, $S$, $T$, etc., and 
> those of $\mathcal{AN}$ by $\mathcal{A}$, $\mathcal{B}$, etc..

We start with applying $\Theta$ to a basis set $\{\psi_\mu\} \equiv \Psi$ 
which engenders an Irrep $\Delta$ of $\mathcal{N}$, namely,

$$
R \psi_\mu = \sum_\nu \psi_\nu \Delta_{\nu\mu}(R), \\

R \Psi = \Psi \Delta(R).
$$

The effect of $\Theta$ on this basis was showed above, and it is summarized by

$$
\Theta R \Psi = \Theta (\Psi \Delta(R)) = (\Theta \Psi) \Delta^*(R)
$$

Next, we define the generalized time reversed set $\Phi \equiv \{\phi_\mu\} = 
\{\mathcal{A} \psi_\mu\}$ such that

$$
R \Phi = \Phi\quad {}^\mathcal{A}\Delta(R),
$$
but, since $\Theta R = R \Theta$, we have

$$
R (\mathcal{A} \Psi) = (\mathcal{A} \Psi) ^\mathcal{A}\Delta(R) = R S \Theta \Psi 
= S (S⁻¹ R S) \Theta \Psi \\
= S \Theta (S⁻¹ R S) \Psi = (S \Theta \Psi) \Delta^*(S⁻¹ R S) = (\mathcal{A} \Psi)
\Delta^*(S⁻¹ R S). \\
\boxed{^\mathcal{A}\Delta(R) = \Delta^*(S⁻¹ R S),}
$$
where $\Delta(S⁻¹ R S)$ is an Irrep conjugate to $\Delta(R)$, since $\mathcal{N}$ 
is a normal subgroup.

We now construct the Rep $\Gamma$ engendered by the combined basis $\mathbf{F} = 
[\Psi, \mathcal{A}\Psi]$, namely,

$$
\boxed{R \mathbf{F} = \mathbf{F} \Gamma(R) = [\Psi \mathcal{A} \Psi] 
\begin{pmatrix}
    \Delta(R) & \mathbf{0} \\
    \mathbf{0} & \Delta^*(S⁻¹ R S)
\end{pmatrix}, \qquad \forall R \in \mathcal{N}.}
$$

We now construct the Rep $\Gamma$ engendered byt the combined basis $\mathbf{F} 
= [\Psi, \mathcal{A}\Psi]$, namely,

Next we apply an operation $\mathcal{B} = \mathcal{A} T \in \mathcal{AN}$, and 
obtain

$$
\mathcal{B} \Psi = \mathcal{A} T \Psi = \mathcal{A} \Psi \Delta(T) = (\mathcal{A} 
\Psi) \Delta^*(T) = (\mathcal{A} \Psi) \Delta^*(\mathcal{A}⁻¹ \mathcal{B}), \\

\mathcal{B} (\mathcal{A} \Psi) = \mathcal{B} \mathcal{A} \Psi = \Psi 
\Delta(\mathcal{B} \mathcal{A}), \qquad \mathcal{BA} \in \mathcal{N}.
$$

We then find

$$
\boxed{\mathcal{B} \mathbf{F} = \mathbf{F} \Gamma(\mathcal{B}) = [\Psi 
\mathcal{A} \Psi] 
\begin{pmatrix}
    \mathbf{0} & \Delta(\mathcal{BA}) \\
    \Delta^*(\mathcal{A⁻¹B}) & \mathbf{0}
\end{pmatrix}, \qquad \forall \mathcal{B} \in \mathcal{AN}.}
$$

> [!WARNING]
> The matrix representatives $\Gamma$ do not obey the ordinary multiplication 
> relations associated with unitary groups, but rather obey
> $$
> \Gamma(R) \Gamma(S) = \Gamma(RS), \qquad \Gamma(R) \Gamma(\mathcal{B}) = 
> \Gamma(R\mathcal{B}), \\
> \Gamma(\mathcal{B}) \Gamma^*(R) = \Gamma(\mathcal{B}R), \qquad \Gamma(\mathcal{B})
> \Gamma^*(\mathcal{C}) = \Gamma(\mathcal{BC}).
> $$

The set of unitary matrices obtained form a *corepresentation* (CoRep) of 
$\mathcal{M}$, derived from the unitary Irrep $\Delta$ of its normal subgroup 
$\mathcal{N}$.

#### Specialization into Grey Groups

Now we are going to focus our attention into the case that interests us. That 
case is the case of grey groups (type II). In this case we just have that 
$\mathcal{A} = \Theta$ and that $\mathcal{N} = \mathcal{G}$, where $\mathcal{G}$ 
is a certain space group, namely,

$$
\mathcal{M} = \mathcal{G} \oplus \Theta \mathcal{G} = \mathcal{G} \otimes \{E,
\Theta\}.
$$

If you construct the Rep $\Gamma$ for a transformation $R \in \mathcal{G}$, you 
obtain

$$
\boxed{R \mathbf{F} = \mathbf{F} \Gamma(R) = [\Psi \quad \Theta \Psi] 
\begin{pmatrix}
    \Delta(R) & \mathbf{0} \\
    \mathbf{0} & \Delta^*(R)
\end{pmatrix}.}
$$

> [!IMPORTANT]
> The time-reversed representation $^\mathcal{A}\Delta$ is **identical** to the 
> complex conjugate representation $\Delta^*$.

> [!CAUTION]
> What is the difference here with the statement at the beginning? Why are we not
> able to impose the conditions where $R$ and $\Theta$ commute?
>
> In my opinion it all due to th fact that $\Psi \not = \Phi = \Theta \Psi$. If 
> that was the case the previous statement could be applied.

Next we do it for any element $\mathcal{B} = \Theta T \in \Theta\mathcal{G}$, and 
obtain

$$
\boxed{\mathcal{B} \mathbf{F} = \mathbf{F} \Gamma(\mathcal{B}) = [\Psi \Theta\Psi] 
\begin{pmatrix}
    \mathbf{0} & \Delta(T\Theta²) \\
    \Delta^*(T) & \mathbf{0}
\end{pmatrix} = [\Psi \Theta\Psi] 
\begin{pmatrix}
    \mathbf{0} & \Delta(T) \\
    \Delta^*(T) & \mathbf{0}
\end{pmatrix}}
$$

> [!IMPORTANT]
> Notice that we are using the bosonic TRS, i.e., $\Theta² = +\mathbb{I}$. If 
> we want to include fermionic TRS we will need to add a minus sing, consequently.

> [!NOTE]
> Do we need to compute this for a general $\mathcal{B} \in \Theta\mathcal{G}$ or 
> are we just fine with $\Theta$?
>
> Since we are interested only on a generator set of $\mathcal{M}$, I think we are 
> good with only considering $\Theta$.

## Finding an explicitly real form of irrep matrices

An explicit real, or physically real, form of a set of irrep matrices is one 
where the associated matrices $D(g)$ have the following property:

$$
D(g) = D^*(g),
$$

for all operations $g$ in the considered group $G$.

The standard listings of irreps is not explicitly real. However, if an irrep is 
either intrinsically real - or has been made into a corep in the complex or 
pseudoreal case - it is always equivalent to an intrinsically real form. That is,
there exists a unitary transform $S$ such that:

$$
S D(g) S^{-1} = S D(g) S^\dagger = D^*(g).
$$

Suppose we can find this unitary transformation $S$ by some means. What we are 
interested in, is finding a related transform $W$, defining an explicitly real 
form of the irrep:

$$
\tilde{D}(g) = W D(g) W^{-1} = W D(g) W^\dagger,
$$

where $W$ is some other unitary transformation and where $\tilde{D}(g)$ is an 
intrinsically real form of $D(g)$, i.e., where

$$
\tilde{D}(g) = \tilde{D}^*(g) \quad \forall g\in G.
$$

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

$$
S = VΛV^{-1} = VΛ^{1/2}Λ^{1/2}V^{\mathrm{T}} = (VΛ^{1/2}V^{\mathrm{T}})
(VΛ^{1/2}V^{\mathrm{T}}),
$$

so we can pick $W = VΛ^{1/2}V^{\mathrm{T}}$ (note also that the square root of 
$\Lambda$ must exist and is well-defined since $S$ is invertible, i.e., has full 
rank).
Hence $W^* = V(Λ^{1/2})^*V^{\mathrm{T}} = VΛ^{-1/2}V^{\mathrm{T}} = W^{-1}$ and 
$W^{\mathrm{T}} = (VΛ^{1/2}V^{\mathrm{T}})^{\mathrm{T}} = 
(V^{\mathrm{T}})^{\mathrm{T}}(Λ^{1/2})^{\mathrm{T}} V^{\mathrm{T}} = 
VΛ^{1/2}V^{\mathrm{T}} = W$. I.e., $W$ is also unitary symmetric and normal.

Now, let us rewrite $S D(g) S^{-1} = D^*(g)$ in terms of $W$:

$$
WW D_g W^{-1}W^{-1} = D_g^* \\
$$

Multiply from LHS by $W^*$ and from RHS by $W$:

$$
W^*WW D_g W^{-1}W^{-1} W = W^*D_g^* W \\
\Leftrightarrow W D_g W^{-1} = W^*D_g^* W \\
\Leftrightarrow W D_g W^{-1} = W^* D_g^* (W^{-1})^*
$$
where we have used the properties of $W$ to reduce the expressions. 

Identifying $\tilde{D}(g) = W D(g) W^{-1}$ we clearly obtain the desired 
invariance under complex conjugation since $\tilde{D}^*(g) = (W D_g W^{-1})^* = 
W^* D_g^* (W^{-1})^* = W D_g W^{-1} = \tilde{D}(g)$.

## Theory of representations in crystalline systems

Here I am going to explain how to present a general Hamiltonian using a basis
and what is the behavior of this representation under symmetries. First of all, 
we will study the action of a symmetry transformation over a basis set in 
$k$-space, and, secondly, the constraints the Hamiltonian matrix, written using
such basis set, will fullfil due to symmetry.

### Representation of symmetry operators using a basis

Following the deductions made by Barry in Ref. [1].

Let us start with a basis set in real space $\{ψ_{iα}
(\mathbf{r})\}$, where $i$ indicates the internal degrees of freedom of 
the orbital, $α$ indicates the site $\mathbf{q}_α$ inside the Wyckoff
position. Notice that by construction we assume each function $ψ_{iα}
(\mathbf{r})$ is localized on $\mathbf{q}_α$. Intuitively, you can think of 
them as Wannier functions.

We are going to focus our attention to a particular orbital $ψ_{i1}(\mathbf{r})$
localized in the site $\mathbf{q}_1 \equiv \mathbf{q}$. This orbital will 
transform under the representation $ρ$ of the site-symmetry group $G_\mathbf{q}$,
associated to $\mathbf{q}$. Then, for each $h \in G_\mathbf{q}$:

$$
h ψ_{i1}(\mathbf{r}) = [ρ(h)]_{ji} ψ_{j1}(\mathbf{r})
$$

Within the primitive unit cell, an orbital localized on each $\mathbf{q}_α$
can be defined as:

$$
ψ_{iα}(\mathbf{r}) = g_α ψ_{i1}(\mathbf{r}) = ψ_{i1}(g_α^{-1} \mathbf{r}),
$$
where $g_α$, with translations, generates the coset decomposition of 
$G_\mathbf{q}$ in $G$. In other words, we can assign for each $\mathbf{q}_α$
a space group element $g \in G$, such that $\mathbf{q}_α = g_α\mathbf{q}$
and:

$$
G = \bigcup_{α=1}^n g_α (G_\mathbf{q} \ltimes T).
$$

By extension, translated counterparts in other unit cells can be defined by:

$$
\{E|\mathbf{t}\} ψ_{iα}(\mathbf{r}) = ψ_{iα}(\mathbf{r-t}),
$$
where $\mathbf{t}$ is a lattice translation. The set of $n \times \text{dim}(ρ)
\times \mathcal{N}$ functions $ψ_{iα}(\mathbf{r-t})$, where $\mathcal{N}$ is
the number of unit cells of the system, will be the basis set on which the induced 
representation $D$ will act.

Specifically, given $g = \{R|\mathbf{v}\} \in G$, the coset decomposition implies 
that for each $g_α$, there is an unique operation $g_β$ such that:

$$
g g_α = \{E|\mathbf{t}_{αβ}\} g_β h,
$$
where $h \in G_\mathbf{q}$ and $t_{αβ} \equiv g\mathbf{q}_α - 
\mathbf{q}_β$.

> [!NOTE]
> Maybe I can prove this in an appendix just for completeness.

Taking all of this into consideration, we can deduce how our basis set will 
transform under the action of every $g \in G$:

$$
g ψ_{iα}(\mathbf{r-t}) = g \{E|\mathbf{t}\} ψ_{iα}(\mathbf{r}) = \\
\{E|R\mathbf{t}\} g ψ_{iα}(\mathbf{r}) = \\
\{E|R\mathbf{t}\} \{E|\mathbf{t}_{αβ}\} g_β h ψ_{i1}(\mathbf{r}) = \\
\{E|R\mathbf{t}\} \{E|\mathbf{t}_{αβ}\} g_β [ρ(h)]_{ji} 
ψ_{j1}(\mathbf{r}) = \\
\{E|R\mathbf{t}\} \{E|\mathbf{t}_{αβ}\} [ρ(h)]_{ji} 
ψ_{jβ}(\mathbf{r}) = \\
\{E|R\mathbf{t}\} [ρ(h)]_{ji} ψ_{jβ}(\mathbf{r-\mathbf{t}_{αβ}}) = \\
[ρ(h)]_{ji} ψ_{jβ}(\mathbf{r}-R\mathbf{\mathbf{t}-\mathbf{t}_{αβ}})
$$

While it is natural to define the representation in real space, it will more 
useful to view it in reciprocal space. This is more evident when $\mathcal{N} 
\to \infty$. To this end, we define the Fourier transform of our basis:

$$
φ_{iα,\mathbf{k}}(\mathbf{r}) = \sum_\mathbf{t} e^{i\mathbf{k\cdot(t+q_α)}} 
ψ_{iα}(\mathbf{r-t}),
$$
where the sum is over all lattice vectors $\mathbf{t} \in T$.

> [!NOTE]
> Notice that this is just convention but for building a tight-binding Hamiltonian
> this choice is better since we will eliminate all local phases.

The Fourier transform amounts to a unitary transformation that exchanges 
$\mathcal{N}$ unit cells for $\mathcal{N}$ distinct $\mathbf{k}$ points. The 
action of $g \in G$ in reciprocal space becomes:

$$
g φ_{iα,\mathbf{k}}(\mathbf{r}) = \sum_\mathbf{t} e^{i\mathbf{k\cdot(t+q_α)}} g
ψ_{iα}(\mathbf{r-t}) = \\
\sum_\mathbf{t} e^{i\mathbf{k\cdot(t+q_α)}} [ρ(h)]_{ji} 
ψ_{jβ}(\mathbf{r}-R\mathbf{\mathbf{t}-\mathbf{t}_{αβ}}) = \\
\sum_\mathbf{t'} e^{i\mathbf{k}\cdot R^{-1}\mathbf{(t'+q_β-v)}} [ρ(h)]_{ji} 
ψ_{jβ}(\mathbf{r-t'}) = \\
e^{-i([R⁻¹]ᵀ \mathbf{k) \cdot v}} [ρ(h)]_{ji} \sum_\mathbf{t'} 
e^{i([R⁻¹]ᵀ \mathbf{k) \cdot (t'+q_β)}} ψ_{jβ}(\mathbf{r-t'}) = \\
e^{-i([R⁻¹]ᵀ \mathbf{k) \cdot v}} [ρ(h)]_{ji} φ_{jβ,[R⁻¹]ᵀ\mathbf{k}}(\mathbf{r}),
$$
where we have made the substitution: $\mathbf{t}' = R\mathbf{t} + \mathbf{t}_{αβ}
= R\mathbf{t} + g\mathbf{q_α - q_β} = R\mathbf{t} + R\mathbf{q_α + v - q_β} =
R(\mathbf{t+q_α}) + \mathbf{v - q_β} \Rightarrow (\mathbf{t+q_α}) = R^{-1}
(\mathbf{t'+q_β-v})$.

> [!NOTE]
> Important note, we have used here an interesting trick that has been a source 
> of confusion. Here, we made the substitution $\mathbf{k}·(R⁻¹\mathbf{r}) \equiv
> (g \mathbf{k)·r}$. Let me prove that:
> $$
> \mathbf{k}·(R⁻¹\mathbf{r}) = \sum_{ij} k_i (R⁻¹_{ij} r_j) = \sum_{ij} (R⁻¹_{ij} k_i)
> r_j = ([R⁻¹]ᵀ \mathbf{k}) · \mathbf{r} \equiv (g \mathbf{k}) · \mathbf{r},
> $$

In reciprocal space, the matrix representation can be interpreted as a $\mathcal{N}
\times \mathcal{N}$ matrix of $n\dim(ρ) \times n\dim(ρ)$ blocks, each block can 
be labeled by $\mathbf{k,k'}$. Most of the blocks are zero: given $g = \{R|
\mathbf{v}\} \in G$, there is only one non-zero block in each row and column, 
corresponding to $\mathbf{k'} = R\mathbf{k}$. Mathematically, we can express this
as:

$$
g φ_{iα,\mathbf{k}}(\mathbf{r}) = \sum_{jβ\mathbf{k'}} D_{jβ\mathbf{k'},iα\mathbf{k}}(g)
φ_{jβ,\mathbf{k}'}(\mathbf{r}),
$$
where we have that:

$$
D_{jβ\mathbf{k'},iα\mathbf{k}}(g) = e^{-i(g\mathbf{k) \cdot v}} ρ_{ji}(h)
\delta_{g\mathbf{k,k'}} \delta_{g\mathbf{q_α - q_β} \mod \tau},
$$
where $\tau \in T$.

We will use the following notation:

$$
Ρ_{jβ,iα}(g) = e^{-i(g\mathbf{k) \cdot v}} ρ_{ji}(h) 
\delta_{g\mathbf{q_α - q_β} \mod \tau},
$$
where we skip the dependence on $\mathbf{k}$ which will be unnecessary.

We can vectorize the previous equation as:

$$
\boxed{g Φ_\mathbf{k}(\mathbf{r}) = Ρ^T(g) Φ_{g\mathbf{k}}(\mathbf{r})},
$$
where $Φ_\mathbf{k}(\mathbf{r})$ is a column vector formed by 
$\{φ_{iα,\mathbf{k}}(\mathbf{r})\}$, and, $Ρ(g)$ is a $n \times n$ matrix of 
$\dim(ρ) \times \dim(ρ)$ blocks, each of them can be labelled by $α,β$. Most of 
the blocks are zero: given $g \in G$, there is only one non-zero block in each 
row and column, corresponding to $g q_α - q_β = 0 \mod \tau$ with $\tau \in T$, 
and it is going to be equal to:

$$
Ρ_{jβ,iα}(g)= e^{-i(g\mathbf{k) \cdot v}} [ρ(h)]_{ji} 
\delta_{g\mathbf{q_α - q_β} \mod \tau}
$$

> [!NOTE]
> We pick the previous definition of the matrix in order to have good properties
> of composition. This is due to the fact that:
> $$
> g₁ g₂ Φ_\mathbf{k}(\mathbf{r}) = Ρ^T(g₁g₂) Φ_{g₁g₂\mathbf{k}}(\mathbf{r}) \\
> = g₁ Ρ^T(g₂) Φ_{g₂\mathbf{k}}(\mathbf{r}) = Ρ^T(g₂) Ρ^T(g₁) 
> Φ_{g₁g₂\mathbf{k}}(\mathbf{r}) \\ 
> \Rightarrow \boxed{Ρ(g₁g₂) = Ρ(g₁) Ρ(g₂)}
> $$

### Action of symmetry operators in a Hamiltonian

Let us start with the most general non-interacting Hamiltonian:

$$
\hat{H} = \sum_{IJ,\mathbf{R}\mathbf{R}'} h_{IJ,\mathbf{R-R}'} \; 
\hat{c}_{I,\mathbf{R}}^\dagger \hat{c}_{J,\mathbf{R}'},
$$
where $I,J$ wrap up the internal degrees of freedom of the orbitals and the sites
of the WP, i.e., $I = (i, α)$; and $\mathbf{R,R}'$ run over the lattice translations.

> [!WARNING]
> Notice that we have assumed that hopping terms only depends on relative distances.
> We are going to denote $\mathbf{t \equiv R - R}'$.

We can apply the same Fourier transform to go into reciprocal space:

$$
\hat{c}_{I,\mathbf{R}} = \frac{1}{\sqrt{N}} \sum_{\mathbf{k}} e^{-i\mathbf{k·(R+q_α)}}
\hat{a}_{I,\mathbf{k}},
$$

obtaining:

$$
\hat{H} = \frac{1}{N} \sum_{IJ,\mathbf{RR}'} h_{IJ,\mathbf{t}} \sum_{\mathbf{kk}'}
e^{i\mathbf{k·(R+q_α)}} e^{-i\mathbf{k'·(R'+q_β)}} \hat{a}^\dagger_{I,\mathbf{k}} 
\hat{a}_{J,\mathbf{k}'} \\
= \frac{1}{N} \sum_{IJ,\mathbf{t},\mathbf{kk}'} h_{IJ,\mathbf{t}} 
\left[ \sum_{\mathbf{R}'} e^{i\mathbf{(k-k')·R+'}} \right] e^{i\mathbf{k·(t+q_α)}} 
e^{-i\mathbf{k'·q_β}} \hat{a}^\dagger_{I,\mathbf{k}} \hat{a}_{J,\mathbf{k}'} \\
= \sum_{IJ,\mathbf{t},\mathbf{kk}'} h_{IJ,\mathbf{t}} 
\delta_{\mathbf{k,k'}} e^{i\mathbf{k·(t+q_α)}} e^{-i\mathbf{k'·q_β}} 
\hat{a}^\dagger_{I,\mathbf{k}} \hat{a}_{J,\mathbf{k}'} \\
= \sum_{IJ,\mathbf{t},\mathbf{k}} h_{IJ,\mathbf{t}} 
e^{i\mathbf{k·(t+q_α-q_β)}} \hat{a}^\dagger_{I,\mathbf{k}} \hat{a}_{J,\mathbf{k}} \\
= \sum_{IJ,\mathbf{k}} h_{IJ,\mathbf{k}} \hat{a}^\dagger_{I,\mathbf{k}} 
\hat{a}_{J,\mathbf{k}},
$$
where he have defined: $h_{IJ,\mathbf{k}} = \sum_\mathbf{t} h_{IJ,\mathbf{t}}
e^{i\mathbf{k·(t+q_α-q_β)}}$.

#### Quantization of the representations

The quantization of the previous (classical) theory of representations can be 
written using "braket" notation as:

$$
\hat{g} \ket{φ_{I,\mathbf{k}}} = Ρ_{JI}(g) \ket{φ_{J,g\mathbf{k}}},
$$
where $\ket{φ_{I,\mathbf{k}}} \equiv a^\dagger_{I,\mathbf{k}} \ket{0}$. Then:

$$
\hat{g} \hat{a}^\dagger_{I,\mathbf{k}} \hat{g}^{-1} \ket{0} = \hat{g} 
\hat{a}^\dagger_{I,\mathbf{k}} \ket{0} = Ρ_{JI}(g) 
\hat{a}^\dagger_{J,g\mathbf{k}} \ket{0} \\
\Rightarrow \boxed{\hat{g} \hat{a}^\dagger_{I,\mathbf{k}} \hat{g}^{-1} = Ρ_{JI}(g) 
\hat{a}^\dagger_{J,g\mathbf{k}}}
$$

> [!WARNING]
> We are using that $\hat{g}^{-1} \ket{0} = \ket{0}$, but I think it is a good 
> assumption that symmetries do not act on vacuum.

Consequently: 
$$
\hat{g} \hat{a}_{I,\mathbf{k}} \hat{g}^\dagger = Ρ^*_{JI}(g)
\hat{a}_{J,g\mathbf{k}} \\
\Rightarrow \boxed{\hat{g} \hat{a}_{I,\mathbf{k}} \hat{g}^{-1} = Ρ^*_{JI}(g) 
\hat{a}_{J,g\mathbf{k}}}
$$

> [!WARNING]
> Here I am passing the operator through but I need to check it if this could be
> done formally. Additionally, we are assuming that $\hat{g}^\dagger = \hat{g}^{-1}$.


Then, if we want the Hamiltonian to be invariant under the symmetries, we must 
impose that:

$$
\hat{H} = \hat{g} \hat{H} \hat{g}^{-1}
$$

Then we obtain that:

$$
\hat{H} = \sum_{IJ,\mathbf{k}} \hat{a}^\dagger_{I,\mathbf{k}} h_{IJ,\mathbf{k}}
\hat{a}_{J,\mathbf{k}} = \\
\hat{g} \hat{H} \hat{g}^{-1} = \sum_{IJ,\mathbf{k}} \hat{g}
\hat{a}^\dagger_{I,\mathbf{k}} h_{IJ,\mathbf{k}} \hat{a}_{J,\mathbf{k}} \hat{g}^{-1} \\
= \sum_{IJ,\mathbf{k}} \hat{g} \hat{a}^\dagger_{I,\mathbf{k}} \hat{g}^{-1} 
h_{IJ,\mathbf{k}} \hat{g} \hat{a}_{J,\mathbf{k}} \hat{g}^{-1} \\
= \sum_{IJ,\mathbf{k},I'J'} \hat{a}^\dagger_{I',g\mathbf{k}} Ρ_{I'I}(g) h_{IJ,\mathbf{k}}
Ρ^*_{J'J}(g) \hat{a}_{J',g\mathbf{k}} \\
= \sum_{IJ,\mathbf{k}} \hat{a}^\dagger_{I,g\mathbf{k}} Ρ(g) H_{\mathbf{k}}
Ρ^\dagger(g) \hat{a}_{J,g\mathbf{k}},
$$
where we have defined $H_\mathbf{k} \equiv h_{IJ,\mathbf{k}}$ and made the substitution
$I',J' \to I,J$. Comparing the first and final rows we obtain the following 
relation for the Hamiltonian to be invariant under symmetries:

$$
\boxed{H_\mathbf{k} = Ρ(g) H_{g^{-1}\mathbf{k}} Ρ^\dagger(g)}
$$

> [!NOTE]
> Notice that the representations of spatial operations are unitary, so we end up 
> with:
> 
> $$
> \boxed{H_\mathbf{k} = Ρ(g) H_{g^{-1}\mathbf{k}} Ρ⁻¹(g)},
> $$
> which is the most usual way to write it.

### Time reversal symmetry

For TRS a similar computation can be performed. Let us assume that the action of
TRS over our basis is the following:

$$
Θ \ket{φ_{I,\mathbf{k}}} = \ket{φ_{I,\mathbf{-k}}},
$$
then we obtain the following relations:

$$
\boxed{\hat{Θ} \hat{a}^\dagger_{I,\mathbf{k}} \hat{Θ}^{-1} = 
\hat{a}^\dagger_{I,-\mathbf{k}} \quad \hat{Θ} \hat{a}_{I,\mathbf{k}} 
\hat{Θ}^{-1} = \hat{a}_{I,-\mathbf{k}}}
$$

Then, the invariance under TRS of the Hamiltonian is simply reduced to:

$$
\hat{Θ} \hat{H} \hat{Θ}^{-1} = \sum_{IJ,\mathbf{k}} \hat{Θ} 
\hat{a}^\dagger_{I,\mathbf{k}} h_{IJ,\mathbf{k}} \hat{a}_{J,\mathbf{k}} \hat{Θ}^{-1} \\
= \sum_{IJ,\mathbf{k}} \hat{Θ} \hat{a}^\dagger_{I,\mathbf{k}} \hat{Θ}^{-1} 
h^*_{IJ,\mathbf{k}} \hat{Θ} \hat{a}_{J,\mathbf{k}} \hat{Θ}^{-1} \\
= \sum_{IJ,\mathbf{k}} \hat{a}^\dagger_{I,-\mathbf{k}} h^*_{IJ,\mathbf{k}} 
\hat{a}_{J,-\mathbf{k}} = \\
\hat{H} = \sum_{IJ,\mathbf{k}} \hat{a}^\dagger_{I,\mathbf{k}} h_{IJ,\mathbf{k}} 
\hat{a}_{J,\mathbf{k}}
$$

Obtaining the following relation:

$$
\boxed{H_\mathbf{k} = H^*_{-\mathbf{k}}}
$$

## Proof that physically real representations need a real basis

If the representation is real then $Ρ(g) = Ρ^*(g)$. On one side we have that: 
$g φ_{I,\mathbf{k}} = Ρ_{JI}(g) φ_{J,g\mathbf{k}}$. On the other hand, since 
$Θ g = g Θ$:

$$
g Θ φ_{I,\mathbf{k}} = Θ (g φ_{I,\mathbf{k}}) = Θ (Ρ_{JI}(g) φ_{J,g\mathbf{k}}) 
= Ρ^*_{JI}(g) (Θ φ_{I,\mathbf{k}}) = Ρ_{JI}(g) (Θ φ_{I,\mathbf{k}}).
$$

Then $Θ φ_{I,\mathbf{k}}$ yields the same representation as $φ_{I,\mathbf{k}} $ 
so they can picked such that $φ_{I,\mathbf{k}} = Θ φ_{I,\mathbf{k}}$

> [!CAUTION]
> Maybe this is too much of an assumption…

## References

[1] Band Representations and Topological Quantum Chemistry by Bradlyn *et al.*
(2021) https://doi.org/10.1146/annurev-conmatphys-041720-124134 