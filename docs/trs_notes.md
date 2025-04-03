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
unitary matrix, i.e., $S = S^{\mathrm{T}}$ and $S^{-1} = S^\dagger$ (implifying, 
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

## Representation of the Hamiltonian using a real basis

Here I am going to explain how to present a general Hamiltonian using a real basis
and what is the behavior of this representation under TRS.

First of all, let me start with a general Hamiltonian $H$. If we assume that the 
system is periodic, we can apply Bloch theorem and decompose the Hamiltonian as:

$$
H = \sum_\mathbf{k} Ψ_\mathbf{k}^\dagger H_\mathbf{k} Ψ_\mathbf{k},
$$
where $Ψ_\mathbf{k}$ is a column operator whose components are single particle
operators.

In our language, this operators will correspond to a set of orbitals that we know
how they transform under the space group operations.

> [!NOTE]
> In the case of interest, $Ψ_\mathbf{k}$ are explicitly real operators. In 
> other words, $Θ Ψ_\mathbf{k} = Ψ_\mathbf{k}$. This also implies that 
> $Ψ_\mathbf{k}^\dagger = Ψ_\mathbf{k}^T$.

> [!WARNING]
> I am assuming that an explicitly real representation $\Leftrightarrow$ an 
> explicitly real basis. I think this is the only possible way but we might
> dedicate some time to prove it.

Since we have a real basis (and a real representation), we can easily study how 
TRS symmetry will interact with another crystalline symmetry $g \in G$. 

### Real space symmetries

Let me first start with the action of an spatial symmetry on the Hamiltonian:

$$
g H =  g \sum_\mathbf{k} Ψ_\mathbf{k}^\dagger H_\mathbf{k} Ψ_\mathbf{k} = 
\sum_\mathbf{k} (g Ψ_\mathbf{k}^\dagger) H_{g\mathbf{k}} (g Ψ_\mathbf{k}) g = \\
\sum_\mathbf{k} Ψ_{g\mathbf{k}}^\dagger D_\mathbf{k}^\dagger(g) H_{g\mathbf{k}} 
D_\mathbf{k}(g) Ψ_{g\mathbf{k}} g
$$

If we want our Hamiltonian to be invariant, the we must impose that:

$$
H g = \sum_\mathbf{k} Ψ_\mathbf{k}^\dagger H_\mathbf{k} Ψ_\mathbf{k} g = g H =
\sum_\mathbf{k} Ψ_{g\mathbf{k}}^\dagger D_\mathbf{k}^\dagger(g) H_{g\mathbf{k}} 
D_\mathbf{k}(g) Ψ_{g\mathbf{k}} g \\
\Rightarrow H_\mathbf{k} = D_\mathbf{k}^\dagger(g) H_{g\mathbf{k}} 
D_\mathbf{k}(g) \Rightarrow H_{g\mathbf{k}} = D_\mathbf{k}(g) H_\mathbf{k} 
D_\mathbf{k}^\dagger(g)
$$

> [!NOTE]
> Notice that the representations of spatial operations are unitary so 
> $D_\mathbf{k}^\dagger = D_\mathbf{k}⁻¹$. Additionally, in this case they are 
> real so $D_\mathbf{k}^\dagger = D_\mathbf{k}⁻¹ = D_\mathbf{k}^T$.

### Time reversal symmetry

For TRS a similar computation can be performed:

$$
Θ H =  Θ \sum_\mathbf{k} Ψ_\mathbf{k}^\dagger H_\mathbf{k} Ψ_\mathbf{k} = 
\sum_\mathbf{k} (Θ Ψ_\mathbf{k}^\dagger) H_{-\mathbf{k}}^* (Θ Ψ_\mathbf{k}) Θ
$$

Assuming that $Ψ_\mathbf{k}$ is a real basis, we should have that $Θ Ψ_\mathbf{k} 
= Ψ_\mathbf{k}$. Then, the invariance with TRS is simply reduced to:

$$
H Θ = \sum_\mathbf{k} Ψ_\mathbf{k}^\dagger H_\mathbf{k} Ψ_\mathbf{k} Θ = Θ H = 
\sum_\mathbf{k} Ψ_\mathbf{k}^\dagger H_{-\mathbf{k}}^* Ψ_\mathbf{k} Θ \\
\Rightarrow H_\mathbf{k} = H_\mathbf{-k}^*
$$

> [!CAUTION]
> The pull request made with the update of this file assumes everything here is 
> right, so the code will proceed under the previous assumptions.

## Proof that physically real representations need a real basis

If the representation is real then $D = D^*$. On one side we have that: $g Ψ = 
D(g) Ψ$. On the other hand, since $Θ g = g Θ$:

$$
g Θ Ψ = Θ (g Ψ) = Θ (D(g) Ψ) = D^*(g) (Θ Ψ) = D(g) (Θ Ψ).
$$

Then $Θ Ψ$ yields the same representation as $Ψ$ so they can be consider equal (?)

> [!WARNING]
> Maybe this is too much of an assumption…