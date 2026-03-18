## Example: 1D bipartite lattice with inversion

Assume we have two sites in a one dimensional lattice of parameter $a=1$ where we place an
inversion-even orbital at the origin denoted by (1a|A); and an inversion-odd orbital at $x=1/2$
denoted by (1b|B).

!!! note "Notation"
    For brevity we will denote orbitals (1a|A) as $a$ and (1b|B) as $b$.
    Additionally, we will indicate with a subscript the unit cell it belongs to.
    For example, $a_0$ will be placed at $x=0$, while $b_1$ will be placed at $x=3/2$ or $a_{-1}$ at $x=-1$.

### Derivation by inspection

These orbitals will transform under inversion symmetry in the following way:

```math
\mathcal{I} a_n = a_{-n}; \quad \mathcal{I} b_n = -b_{-n-1}.
```

Then the most general inversion-symmetric Hamiltonian for first nearest neighbors is:

```math
\mathcal{H} = \sum_n t (a_n^\dagger b_n - a_n^\dagger b_{n-1}) + \text{c.c.}.
```

It is easy to check that this Hamiltonian is inversion symmetric.

If we consider the following Fourier transform: $a_n = \frac{1}{\sqrt{N}} \sum_n e^{-ikn} a_k$,
and $b_n = \frac{1}{\sqrt{N}} \sum_n e^{-ik(n+1/2)} b_k$
then the Hamiltonian in $k$-space will look like:

```math
\mathcal{H} = \sum_k 2it \sin(k/2) a_k^\dagger b_k + \text{c.c}.
```

### Derivation from our method

First, consider the translation $t=0$. Then, remember that $\Delta_{\alpha\to\beta+R} =
\mathbf{q}_\beta + \mathbf{R} - \mathbf{q}_\alpha$.

- Possible $\Delta$'s: $\Delta_{a\to a} = \Delta_{b\to b} = 0; \quad \Delta_{a\to b} =
    -\Delta_{b\to a} = -1/2$.
- Orbits of $\Delta$'s: $\{\Delta_{a\to a}\} = \{\Delta_{b\to b}\} = \{0\}; \quad
    \{\Delta_{a\to b}\} = \{\Delta_{b\to a}\} = \{1/2,-1/2\}$.

Then the non-symmetrized Hamiltonian for this translation will be:

```math
H_{k,0} = \begin{pmatrix} t_{a\to a} & t_{b\to a}^1e^{ik/2}+t_{b\to a}^2e^{-ik/2} \\
t_{a\to b}^1e^{ik/2}+t_{a\to b}^2e^{-ik/2} & t_{b\to b} \end{pmatrix}.
```

Let's proceed now to symmetrize this Hamiltonian. Since we only have inversion we only need
to check that: $H_{k,0} = \rho(\mathcal{I})H_{k,0}\rho^{-1}(\mathcal{I})$, where:

```math
\rho(\mathcal{I}) = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}.
```

Then, that constraint imposes that: $\begin{cases} t_{b\to a}^1 = -t_{b\to a}^2 \\
t_{a\to b}^1 = -t_{a\to b}^2 \end{cases}$, so the Hamiltonian will look like:

```math
H_{k,0} = \begin{pmatrix} t_{a\to a} & 2it_{b\to a} \sin(k/2) \\
2it_{a\to b} \sin(k/2) & t_{b\to b} \end{pmatrix},
```

which is exactly the Hamiltonian deduced with the previous method but with onsite terms and without hermiticity imposed.
