# Phase conventions in `symmetry_eigenvalues`: reconciling with Crystalline.jl

This document explains the phase convention mismatch between the Convention 1 derivation in
[`theory.md`](../theory.md) and Crystalline.jl's `calc_bandreps` and `lgirreps` and how
`symmetry_eigenvalues` in SymmetricTightBinding.jl corrects for it, to align with Crystalline.jl's convention.

## Background: what `theory.md` derives

The derivation in `theory.md` (§ "Transformation properties under symmetry operations")
gives the Convention 1 symmetry eigenvalue for band $n$ at $\mathbf{k}$ under a little-group
operation $g = \{R|\mathbf{v}\}$:

```math
\chi_{n\mathbf{k}}(g)
= (\Theta_\mathbf{G} \, \mathbf{w}_{n\mathbf{k}})^\dagger \, \mathbf{D}_\mathbf{k}(g) \, \mathbf{w}_{n\mathbf{k}}
```

where:

- $\mathbf{w}_{n\mathbf{k}}$ is the eigenvector of $H(\mathbf{k})$ in the Convention 1 coefficient basis
- $\mathbf{G} = g\mathbf{k} - \mathbf{k}$ is a reciprocal lattice vector (since $g$ is in the little group $G_\mathbf{k}$)
- $[\Theta_\mathbf{G}]_{II} = e^{-i\mathbf{G}\cdot\mathbf{q}_I}$ (so $\Theta_\mathbf{G}^\dagger = \Theta_{-\mathbf{G}}$ contributes $e^{+i\mathbf{G}\cdot\mathbf{q}_I}$ upon conjugation)
- $\mathbf{D}_\mathbf{k}(g) = e^{-i(g\mathbf{k})\cdot\mathbf{v}} \, \rho(h)$ is the site-induced space group representation (Convention 1)

The `SiteInducedSGRepElement` functor in `site_representations.jl` faithfully
implements $\mathbf{D}_\mathbf{k}(g)$:

```julia
Dₖ = cispi(-2dot(gk, v)) * sgrep.ρ   # = e^{-2πi(gk)·v} ρ(h)
```

## What Crystalline.jl's `calc_bandreps` computes

In `Crystalline/src/calc_bandreps.jl` (line 179):

```julia
χᴳₖ += cis(2π*dot(kv′, tα′α′)) * χs[site_symmetry_index]
```

The phase uses `cis(+2π...)` where the Cano/Elcoro paper[^1] uses `cis(-2π...)`.
A comment in the source (lines 180–185) explicitly acknowledges this:

> "The sign in this `cis(...)` above is different from in Elcoro's. I think this is
> consistent with our overall sign convention (see #12), however, and flipping the sign
> causes problems for the calculation of some subductions to `LGIrrep`s."

[^1]: B. Bradlyn et al., "Topological quantum chemistry," Nature **547**, 298 (2017);
L. Elcoro et al., "Double crystallographic groups [...]," J. Appl. Cryst. **50**, 1457 (2017).

This means the band representation characters computed by `calc_bandreps` use conjugated
phase signs for both the $\Theta_\mathbf{G}$ factor and the $\mathbf{D}_\mathbf{k}$ global phase,
relative to the Convention 1 derivation. Crystalline.jl is **internally consistent**: its
`LGIrrep` matrices also use these conjugated signs, so subduction (decomposing characters
into irreps) works correctly within Crystalline — the conjugations cancel when matching
characters against irreps, because both sides are conjugated identically.

## The net effect: complex conjugation

The two sign flips (one in $\Theta_\mathbf{G}$, one in $\mathbf{D}_\mathbf{k}$) both
individually conjugate phases in the character formula, and together they give:

```math
\chi_\text{Crystalline}(g) = \overline{\chi_\text{Convention 1}(g)}
```

That is, the Crystalline.jl convention returns characters that are the **complex conjugate**
of the Convention 1 result.

To see this concretely: the Convention 1 character is
$\chi_\text{Convention 1} = (\Theta_\mathbf{G} \mathbf{w})^\dagger \mathbf{D}_\mathbf{k} \mathbf{w}$.
Under the Crystalline sign convention, $\Theta_\mathbf{G} \to \Theta_{-\mathbf{G}}$ and
$e^{-i(g\mathbf{k})\cdot\mathbf{v}} \to e^{+i(g\mathbf{k})\cdot\mathbf{v}}$, which
individually conjugates both scalar factors in the formula, so the overall character is
conjugated.

## Where the mismatch occurs

The mismatch arises at the **interface** between the two packages:

1. **SymmetricTightBinding.jl** computes symmetry eigenvalues from Hamiltonian
   eigenvectors in Convention 1
2. **Crystalline.jl** expects characters in its conjugated convention (to match its
   `LGIrrep`s)

When `collect_compatible` calls `symmetry_eigenvalues` and passes the result to Crystalline's
`collect_compatible(symeigsv, cbr.brs)`, the conventions must agree. This mismatch was the
root cause of the incorrect irrep labels reported in PR #89.

## The fix

Since $\chi_\text{Crystalline} = \overline{\chi_\text{Convention 1}}$, the correction is simply to
take the complex conjugate of the Convention 1 formula. In `symmetry_eigenvalues`:

```julia
Θᴳ = reciprocal_translation_phase(orbital_positions(ptbm), G)   # = Θ_G (positive G)
ρ = sgrep(k)   # = D_k(g) = e^{-2πi(gk)·v} ρ(h) (Convention 1)
for (n, v) in enumerate(eachcol(vs))
    v_kpG = Θᴳ * v
    χ = dot(v_kpG, ρ, v)          # Convention 1: (Θ_G w)† D_k w
    χ_Crystalline = conj(χ)       # [⚠️ phase]: convert to Crystalline.jl's convention
    symeigs[j, n] = χ_Crystalline
end
```

Note: the eigenvectors `vs` are from `solve(ptbm, k; bloch_phase=Val(false))`, i.e.,
they are eigenvectors of $H(\mathbf{k})$ in the Convention 1 coefficient basis. The
Hamiltonian phase convention affects the eigenvectors (particularly at non-TRIM
$\mathbf{k}$-points where $\mathbf{k} \neq -\mathbf{k}$), so the symmetry eigenvalue
result is coupled to the Hamiltonian phase convention. The code uses the Convention 1
Hamiltonian phase ($e^{+i\mathbf{k}\cdot\boldsymbol{\delta}}$, i.e., `cispi(-2k·δ)` in
`types.jl`); changing this sign would require re-deriving the correction here.

## Should this be fixed in Crystalline.jl instead?

There are three options, each with different trade-offs:

### Option A: Complex-conjugate in `symmetry_eigenvalues` (current approach)

`symmetry_eigenvalues` computes the Convention 1 character and conjugates it before
returning.

- **Pro:** Minimal blast radius; no Crystalline.jl changes needed; the correction
  ($\overline{\chi_\text{Convention 1}} = \chi_\text{Crystalline}$) is concise and easy to
  understand
- **Con:** The returned characters differ by a global complex conjugation from what
  `theory.md` derives; downstream callers must be aware of this convention

### Option B: Fix the root cause in Crystalline.jl (best, but harder)

Change Crystalline.jl to use the conventional minus sign in its exponential phase factors
(fixing https://github.com/thchr/Crystalline.jl/issues/12), then remove the conjugation
correction here entirely.

- **Pro:** Both packages use Convention 1 throughout; the code matches the literature
- **Con:** It is not clear what fixing issue #12 means for irrep label assignments across
  Crystalline.jl.
