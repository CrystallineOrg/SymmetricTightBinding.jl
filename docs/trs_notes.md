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
R \Phi = \Phi ^\mathcal{A}\Delta(R),
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
\boxed{R \mathbf{F} = \mathbf{F} \Gamma(R) = [\Psi \Theta \Psi] 
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

## "Realification" of Irrep and basis

What happens if we "realify" the irreps? In fact, in Crystalline.jl there is a 
function implemented called `realify` which seems as what we are looking for - 
better discuss this with @thchr.

The idea focus on the fact that if $\Psi = \Theta \Psi$ and $\Delta = \Delta^*$,
the action of TRS in the system is trivial - as we implemented already for real 
Irreps. Of course for the case of complex and pseudoreal Irreps, it is not so 
simple but we can make some adaptations to make it simpler.

If we manage to realify every point group Irrep and compute using them a Rep for 
the space group in a given (now real) basis, then the application of TRS will be
as simple as the one showed above for a non-unitary operation $\mathcal{B} = 
\Theta T$, but taking $T = \mathbb{I}$.