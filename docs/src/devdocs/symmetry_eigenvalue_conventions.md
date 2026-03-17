# Phase conventions in `symmetry_eigenvalues`: reconciling with Crystalline.jl

This document explains the phase convention mismatch between the physical derivation in
[`theory.md`](../theory.md) and Crystalline.jl's `calc_bandreps`, and how
`symmetry_eigenvalues` in SymmetricTightBinding.jl corrects for it.

## Background: what `theory.md` derives

The derivation in `theory.md` (§ "Transformation properties under symmetry operations")
gives the **physical** symmetry eigenvalue for band $n$ at $\mathbf{k}$ under a little-group
operation $g = \{R|\mathbf{v}\}$:

```math
\langle\psi_{n\mathbf{k}}|\hat{g}|\psi_{n\mathbf{k}}\rangle
= \sum_{IJ} (w_{I,n\mathbf{k}})^* \, e^{+i\mathbf{G}\cdot\mathbf{q}_\alpha} \, [\mathbf{D}_\mathbf{k}(g)]_{IJ} \, w_{J,n\mathbf{k}}
```

where:

- $\mathbf{w}_{n\mathbf{k}}$ is the eigenvector of $H(\mathbf{k})$
- $\mathbf{G} = g\mathbf{k} - \mathbf{k}$ is a reciprocal lattice vector (since $g$ is in the little group $G_\mathbf{k}$)
- The factor $e^{+i\mathbf{G}\cdot\mathbf{q}_\alpha}$ arises from the conjugation of $\Theta_\mathbf{G}$, defined as $[\Theta_\mathbf{G}]_{II} = e^{-i\mathbf{G}\cdot\mathbf{q}_\alpha}$
- $\mathbf{D}_\mathbf{k}(g) = e^{-i(g\mathbf{k})\cdot\mathbf{v}} \, \rho(h)$ is the site-induced space group representation, with $\rho(h)$ the momentum-independent matrix part

In vectorized form (the boxed formula in `theory.md`):

```math
\langle\psi_{n\mathbf{k}}|\hat{g}|\psi_{n\mathbf{k}}\rangle
= (\Theta_\mathbf{G} \, \mathbf{w}_{n\mathbf{k}})^\dagger \, (\mathbf{D}_\mathbf{k}(g) \, \mathbf{w}_{n\mathbf{k}})
```

All signs follow consistently from the Convention 1 Fourier transform used throughout the
package. The `SiteInducedSGRepElement` functor in `site_representations.jl` faithfully
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

This means the band representation characters computed by `calc_bandreps` are **complex
conjugated** in their global phases relative to the physical convention. However,
Crystalline.jl is **internally consistent**: its `LGIrrep` matrices also use this conjugated
convention, so subduction (decomposing characters into irreps) works correctly — the
conjugation cancels when matching characters against irreps, because both sides are
conjugated.

## Where the mismatch occurs

The mismatch arises at the **interface** between the two packages:

1. **SymmetricTightBinding.jl** computes symmetry eigenvalues from actual Hamiltonian
   eigenvectors using the **physical** convention (matching `theory.md`)
2. **Crystalline.jl** expects characters in its **conjugated** convention (matching its
   `LGIrrep`s)

When `collect_compatible` calls `symmetry_eigenvalues` and passes the result to Crystalline's
`collect_compatible(symeigsv, cbr.brs)`, the conventions disagree. This produced incorrect
irrep labels at many high-symmetry $\mathbf{k}$-points (the bug tracked in PR #89).

## The two components of the mismatch

The conjugation affects two independent phase factors in the symmetry eigenvalue formula:

### Component 1: the $\Theta_\mathbf{G}$ factor

| | Physical (`theory.md`) | Crystalline convention |
|---|---|---|
| Phase on position $\mathbf{q}$ | $e^{+i\mathbf{G}\cdot\mathbf{q}}$ (from $\Theta_\mathbf{G}^\dagger$) | $e^{-i\mathbf{G}\cdot\mathbf{q}}$ (i.e., $\Theta_\mathbf{G}$ unconjugated) |

For operations where $\mathbf{G} = 0$ (when $g\mathbf{k} = \mathbf{k}$ exactly), this
component is trivially 1 and doesn't matter. It matters at $\mathbf{k}$-points like K in p3,
where 3-fold rotations give $g\mathbf{k} = \mathbf{k} + \mathbf{G}$ with nonzero
$\mathbf{G}$.

### Component 2: the global phase from the translation part of $g$

| | Physical (`theory.md`) | Crystalline convention |
|---|---|---|
| Phase on translation $\mathbf{v}$ | $e^{-i(g\mathbf{k})\cdot\mathbf{v}}$ | $e^{+i(g\mathbf{k})\cdot\mathbf{v}}$ |

For operations with $\mathbf{v} = 0$ (pure rotations), this is trivially 1. It matters for
screw axes and glide planes (e.g., the $4_1$ screw in SG 93, 141, etc.).

### Failure pattern explained

This two-component structure explains the pattern of test failures across space groups:

- **Original code (both components wrong):** All 2D p3/p6 failures at the K-point
  (Component 1), plus all 3D screw/glide failures (Component 2)
- **$\Theta_\mathbf{G}$ fix only (Component 1 fixed):** 2D fixed, but 3D screw/glide
  operations still fail
- **Both components fixed:** Everything passes except pre-existing centered-lattice bugs
  (a separate issue, likely in Crystalline.jl)

## The fix

With both corrections, `symmetry_eigenvalues` computes:

```math
\mathbf{w}_{n\mathbf{k}}^\dagger \,
\Theta_\mathbf{G} \,
\tilde{\mathbf{D}}_\mathbf{k}(g) \,
\mathbf{w}_{n\mathbf{k}}
```

where $\tilde{\mathbf{D}}_\mathbf{k}(g) = e^{+2\pi i(g\mathbf{k})\cdot\mathbf{v}} \, \rho(h)$ uses the conjugated global phase but the same (unconjugated) matrix part $\rho$.

Note that only the scalar phases are conjugated, **not** the representation matrix $\rho$.
This matches what `calc_bandreps` does: it conjugates the exponential phase factor (line 179
uses `cis(+2π...)`) but does not conjugate the site-irrep character $\chi_s$.

In the code, this is implemented as:

```julia
# Component 1: use conj(Θ_G) = Θ_{-G} instead of Θ_G
Θᴳ_conj = reciprocal_translation_phase(orbital_positions(ptbm), -G)

# Component 2: flip the sign of the global phase e^{-2πi(gk)·v} → e^{+2πi(gk)·v}
v_g = translation(g)
phase_correction = cispi(4dot(gk, v_g))
ρ = phase_correction * sgrep(k)

# Compute w† Θ_G D̃_k w (matching Crystalline.jl's convention)
for (n, v) in enumerate(eachcol(vs))
    v_kpG = Θᴳ_conj * v
    symeigs[j, n] = dot(v_kpG, ρ, v)
end
```

## Should this be fixed in Crystalline.jl instead?

There are three options, each with different trade-offs:

### Option A: Monkey-patch in `symmetry_eigenvalues` (current approach)

The `theory.md` convention is "correct physics" and the code corrects at the interface with
Crystalline.jl.

- **Pro:** Minimal blast radius; no Crystalline.jl changes
- **Con:** The code's comments reference the `theory.md` formula but compute something
  different; the correction factor (`cispi(4dot(gk, v_g))`) is a wart that's easy to get
  confused about

### Option B: Fix the root cause in Crystalline.jl

Change `calc_bandreps` to use `cis(-2π*dot(kv′, tα′α′))` (matching the Elcoro paper).

- **Pro:** Everything matches the physics and the literature
- **Con:** The comment in `calc_bandreps` warns that "flipping the sign causes problems for
  the calculation of some subductions to `LGIrrep`s." This means the `LGIrrep` matrices in
  Crystalline.jl are also in the conjugated convention. Fixing `calc_bandreps` alone would
  break the subduction machinery; you'd need to also fix the `LGIrrep` phase conventions.
  This is a **large, risky change** to a foundational package with high regression risk.

### Option C: Change the `SiteInducedSGRepElement` convention (recommended)

The key observation is that **the functor's global phase is only used in
`symmetry_analysis.jl`** — the constraint-building code in `tightbinding.jl` uses
`sgrep_induced_by_siteir_excl_phase`, which returns just $\rho$ without any global phase.

So changing the functor's phase convention from $e^{-i(g\mathbf{k})\cdot\mathbf{v}}$ to
$e^{+i(g\mathbf{k})\cdot\mathbf{v}}$ would:

- Remove Component 2 of the monkey-patch from `symmetry_eigenvalues`
- Not affect Hamiltonian construction at all
- Make the `SiteInducedSGRepElement` explicitly match Crystalline.jl's convention

You'd still need the Component 1 fix (the $\Theta_\mathbf{G}$ sign), but that becomes a more
natural choice: implement $\mathbf{w}^\dagger \Theta_\mathbf{G} \tilde{\mathbf{D}}_\mathbf{k} \mathbf{w}$
instead of $(\Theta_\mathbf{G} \mathbf{w})^\dagger \mathbf{D}_\mathbf{k} \mathbf{w}$, documented clearly
as "we use $\Theta_\mathbf{G}$ (not $\Theta_\mathbf{G}^\dagger$) to match Crystalline.jl's
character convention."

In this option, `theory.md` should also gain a note (or a new devdoc section) explaining:
"The physical formula uses $e^{-i(g\mathbf{k})\cdot\mathbf{v}}$, but Crystalline.jl's
`calc_bandreps` and `LGIrrep`s use the conjugated phase $e^{+i(g\mathbf{k})\cdot\mathbf{v}}$.
Since `symmetry_eigenvalues` must match Crystalline.jl's irreps for subduction to work, we
adopt the conjugated convention for $\mathbf{D}_\mathbf{k}$."

## The centering fix: `primitivize` with `modw=false`

Beyond the phase convention mismatch, a separate bug affected all centered lattices (C- and
I-centered space groups: SG 68, 88, 141, 142, 214, 220, 230). The original
`collect_compatible` primitivized little groups via `primitivize(group(first(lgirs)))`, which
calls `primitivize(g::LittleGroup)` with the default `modw=true`, reducing translations
modulo the primitive lattice. For centered lattices, the primitivization transforms
conventional translations $\mathbf{v}_\text{conv}$ to primitive translations
$P^{-1}\mathbf{v}_\text{conv}$, then reduces mod 1. The discarded lattice vector
$\mathbf{R}$ changes the phase $e^{2\pi i(g\mathbf{k})\cdot\mathbf{v}}$ by a factor
$e^{2\pi i(g\mathbf{k})\cdot\mathbf{R}} \neq 1$ at non-$\Gamma$ high-symmetry points.

**Concrete example:** SG 88 ($I4_1/a$), operation $\{-4^+_{001}|\frac{1}{4},\frac{3}{4},\frac{3}{4}\}$ at the P-point $[\frac{1}{2},\frac{1}{2},\frac{1}{2}]$:
- Conventional translation: $\mathbf{v} = [\frac{1}{4}, \frac{3}{4}, \frac{3}{4}]$
- Primitivized (full): $P^{-1}\mathbf{v} = [\frac{3}{2}, 1, 1]$
- Primitivized (reduced mod 1): $[\frac{1}{2}, 0, 0]$ — discards lattice vector $\mathbf{R} = [1, 1, 0]$
- $g\mathbf{k} \cdot \mathbf{R} = -\frac{1}{4}$, giving phase error $e^{2\pi i \times 0.25} = i$

The fix uses `primitivize(::Collection{LGIrrep})`, which internally passes `modw=false`,
preserving full (unreduced) translations:

```julia
clgirsv = irreps(cbr)
lgirsv = primitivize.(clgirsv)
lgs = group.(lgirsv)
```

## Note: the Hamiltonian Fourier phase is independent

The Hamiltonian evaluation phase ($e^{-2\pi i \mathbf{k}\cdot\boldsymbol{\delta}}$ in
`evaluate_tight_binding_term!`) does **not** affect symmetry analysis results. This was
verified during the investigation: the phase sign only transposes $H(\mathbf{k}) \to H(-\mathbf{k})$,
and eigenvectors at the high-symmetry $\mathbf{k}$-points used by `symmetry_eigenvalues` are
obtained from `solve(ptbm, k)` which always uses the actual $\mathbf{k}$. Changing the
Hamiltonian phase sign is therefore not a route to fixing symmetry analysis.

## Test results with all fixes

Across all 230 3D space groups (and all 1D and 2D space groups):

- **All 1D:** pass
- **All 2D:** pass (including p3, p6, which previously failed at the K-point)
- **All 3D:** pass (including all centered lattices)
