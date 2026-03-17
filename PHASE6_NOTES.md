# Phase 6 Investigation Notes: Symmetry Analysis Fix

## Problem Summary

`symmetry_eigenvalues` returns incorrect irrep labels at many high-symmetry k-points.
The failures are especially at K/KA/H/HA (3-fold symmetry), but also at other k-points
in 3D space groups (P, PA, Y, T, etc.).

## Key formula (from theory.md)

The symmetry eigenvalue for band n at k-point k under operation g is:
```
⟨ψ_{nk}|ĝ|ψ_{nk}⟩ = (Θ_G w_{nk})† (D_k(g) w_{nk})
                     = w† Θ_G† D_k(g) w
```
where:
- `w_{nk}` = eigenvector of H_k
- `D_k(g) = e^{-2πi(gk)·v} ρ(h)` is the site-induced SG rep
- `Θ_G = diag(e^{-2πi G·q_α})` with G = gk - k
- The code uses `reciprocal_translation_phase(positions, G)` for Θ_G

## Current state of the branch

Changes already made (in original diff, before this exploration):
1. Flipped sign in Hamiltonian evaluation phase: `cispi(-2k·δ)` → `cispi(+2k·δ)`
   (types.jl:419, tightbinding.jl, gradients.jl)
2. Updated print code for conjugate blocks
3. Ported better test infrastructure to test/symmetry_analysis.jl
4. Commented out berry.jl test that compares against manual Haldane model

## Key observations from PR #89

1. Complex conjugating the computed symmetry eigenvalues fixes K/KA/H/HA failures
   (all 2D failures) but not all 3D failures (e.g., SG 68 remains)
2. Changing sign in Crystalline.jl's `calc_bandreps` also fixes K/KA but breaks
   integer-valued subduction elsewhere
3. The issue exists even for k-independent Hamiltonians (p3, (1c|A), on-site only)

## Root cause analysis

The mismatch is between two conventions:
- **Physical (Convention 1)**: `⟨ψ|ĝ|ψ⟩ = w† Θ_G† D_k(g) w` where
  `D_k(g) = e^{-2πi(gk)·v} ρ(h)` and `Θ_G = diag(e^{-2πiG·q_α})`
- **calc_bandreps (Crystalline.jl)**: Uses `cis(+2π·dot(gk, t))` instead of `cis(-2π...)`
  (acknowledged in Crystalline.jl issue #12, calc_bandreps.jl line 180-185 comment).
  Effectively computes the complex conjugate of the physical phases.

The mismatch has TWO components:
1. **Θ_G factor**: `e^{-2πiG·q}` (physical) vs `e^{+2πiG·q}` (calc_bandreps)
2. **Global phase**: `e^{-2πi(gk)·v}` (physical, from `SiteInducedSGRepElement` functor)
   vs `e^{+2πi(gk)·v}` (calc_bandreps)

For operations with v=0 (pure rotations), only component 1 matters. This is why the 2D
plane groups p3/p6 (which have only rotations at K) are fixed by the Θ_G fix alone.

## Exploration log

### Finding 0: Hamiltonian sign flip is irrelevant to symmetry analysis

Verified that the Hamiltonian phase sign flip (`cispi(-2k·δ)` → `cispi(+2k·δ)` in
types.jl/tightbinding.jl/gradients.jl) has NO effect on symmetry analysis results —
identical test output with or without it. This makes sense: the sign flip just transposes
H(k) → H(-k), and eigenvectors at the high-symmetry k-points used by `symmetry_eigenvalues`
are obtained from `solve(ptbm, k)` which always uses the actual k.

Reverted the Hamiltonian sign flip (and related print/gradients changes) back to original.
Re-enabled the Haldane model comparison test in berry.jl.

### Attempt 1: Fix Θ_G sign only (use -G instead of G)

Changed `symmetry_analysis.jl` line 117:
```julia
Θᴳ_conj = reciprocal_translation_phase(orbital_positions(ptbm), -G) # = conj(Θᴳ)
```
(was using `G` before, now uses `-G` to get `conj(Θ_G)` = `e^{+2πiG·q}`)

**Result (max_sgnum_3d=130):**
- Fixed: All 1D, all 2D (SG 13, 14 — p3, p6)
- Still failing: SG 68, 88, 93, 98, 108, 112, 116, 118, 122
- Pattern: remaining failures all involve operations with non-zero translation v
  (screws, glides), where the second convention mismatch (global phase) matters

### Attempt 2: Also fix global phase sign (phase correction in symmetry_eigenvalues)

Added phase correction within `symmetry_eigenvalues` to flip the sign of `(gk)·v`:
```julia
v_g = translation(g)
phase_correction = cispi(4dot(gk, v_g))  # flips e^{-2πi(gk)·v} → e^{+2πi(gk)·v}
ρ = phase_correction * sgrep(k)
```
Applied within `symmetry_eigenvalues` rather than modifying the `SiteInducedSGRepElement`
functor (which is also used for constraint building in `tightbinding.jl`).

**Result (max_sgnum_3d=130):**
- Fixed additionally: SG 93, 98, 108, 112, 116, 118, 122
- Still failing: SG 68, 88 only

**Result (max_sgnum_3d=230, full scan):**
- Failing SGs: 68, 88, 141, 142, 214, 220, 230 (7 total out of 230)
- Two failure modes:
  1. "symmetry vector discrepancy" (SG 68 positions 8h/8d/8c, SG 220 positions 12a/12b):
     irreps are swapped (e.g., Y⁺↔Y⁻ parity, T₁↔T₂, P₁↔P₂)
  2. "failed to collect" (SG 68 4a/4b, SG 88 4a/4b, SG 141 4a/4b, SG 142 8b/16e,
     SG 214 8a/8b/12c/12d, SG 230 16b/24c): no valid irrep decomposition found at all

**Pattern in remaining failures:**
- All failing SGs have non-primitive centering: SG 68 is C-centered, rest are I-centered
- Many involve 4₁ screw axes (SG 88=I4₁/a, 141=I4₁/amd, 142=I4₁/acd, 214=I4₁32)
- SG 220=I-43d and SG 230=Ia-3d are body-centered cubic
- Likely a separate class of bug related to centering and/or orbit computation

### Baseline comparison: remaining failures are PRE-EXISTING

Ran the same 7 failing SGs (68, 88, 141, 142, 214, 220, 230) with the **original code**
(no Θ_G fix, no phase correction). Results:

| SG  | Original failures | With fix | Notes |
|-----|------------------|----------|-------|
| 68  | 14 (6 disc + 8 fc) | 14 (6 disc + 8 fc) | Same failures |
| 88  | 4 fc             | 4 fc     | Same failures |
| 141 | 8 fc             | 8 fc     | Same failures |
| 142 | 4 fc             | 4 fc     | Same failures |
| 214 | 10 fc            | 8 fc     | Fix resolves 8a/E, 8b/E |
| 220 | 6 (all fc)       | 4 (all disc) | Fix resolves 16c; 12a/12b change from fc→disc |
| 230 | 9 fc             | 4 fc     | Fix resolves 16a (4 EBRs), 16b/E |

**Conclusion:** All remaining failures are pre-existing bugs, likely in Crystalline.jl's
handling of centered lattices (all failing SGs are C- or I-centered). My fix actually
*improves* the situation for SG 214, 220, and 230. These pre-existing failures are a
separate class of bug from the phase convention issue and should be tracked independently.

### Summary of fix

The fix to `symmetry_eigenvalues` has two components, both in `src/symmetry_analysis.jl`:

1. **Θ_G sign**: Use `-G` instead of `G` in `reciprocal_translation_phase`, giving
   `conj(Θ_G) = e^{+2πiG·q}` instead of `Θ_G = e^{-2πiG·q}`
2. **Global phase**: Multiply by `cispi(4dot(gk, v))` to flip the sign of the
   `e^{-2πi(gk)·v}` term from the `SiteInducedSGRepElement` functor to
   `e^{+2πi(gk)·v}`, matching `calc_bandreps`' convention

Both corrections are needed to match Crystalline.jl's `calc_bandreps` convention, which
uses `cis(+2π...)` where the physics gives `cis(-2π...)` (Crystalline.jl issue #12).

**Test results with fix (all 230 SGs):**
- All 1D: PASS
- All 2D: PASS
- 3D: 223/230 pass, 7 fail (all pre-existing, centered lattice bugs)

### Cleanup: reverted Hamiltonian sign flip

The WIP commit had changed `cispi(-2k·δ)` → `cispi(+2k·δ)` in types.jl, tightbinding.jl,
and gradients.jl. This was the original hypothesis for fixing the symmetry analysis bug,
but it turned out to be wrong — it broke berry.jl tests (Berry phase = -π instead of π,
wrong 3D Berry curvature values, wrong Chern numbers). Reverted tightbinding.jl and
gradients.jl back to original `cispi(-2k·δ)`. types.jl was already reverted.

Also added `using LinearAlgebra: dot` to berry.jl for the Haldane model comparison test.

### Verification: fix doesn't break other tests

- Berry curvature tests: 19/19 PASS
- Chern number tests: 27/27 PASS (including Haldane model comparison)
- Symmetry analysis: same results as before (all pass except pre-existing centered failures)

### Attempt 3: Fix centered-lattice failures via `primitivize(::Collection{LGIrrep})`

The remaining 7 failing SGs (68, 88, 141, 142, 214, 220, 230) were NOT pre-existing
Crystalline.jl bugs — they were caused by a second bug in `collect_compatible`.

**Root cause:** `collect_compatible` primitivized the little groups via
`primitivize(group(first(lgirs)))`, which calls `primitivize(g::LittleGroup)` with the
default `modw=true`. This reduces translations modulo the primitive lattice. For centered
lattices, the primitivization transforms conventional translations `v_conv` to primitive
translations `P⁻¹ v_conv`, then reduces mod 1. The discarded lattice vector `R` changes
the phase `e^{2πi(gk)·v}` by a factor `e^{2πi(gk)·R}`, which is ≠ 1 when `gk` is not an
integer vector (i.e., at most non-Γ high-symmetry points).

**Concrete example:** SG 88, operation `{-4⁺₀₀₁|¼,¾,¾}` at the P-point `[½,½,½]`:
- Conventional translation: `v = [1/4, 3/4, 3/4]`
- Primitivized (full): `P⁻¹ v = [3/2, 1, 1]`
- Primitivized (reduced mod 1): `[1/2, 0, 0]` — discards lattice vector `R = [1, 1, 0]`
- `gk · R = -1/4`, giving phase error `exp(2πi × 0.25) = i`

**Fix:** Use Crystalline's `primitivize(::Collection{LGIrrep})` which internally calls
`primitivize(g::LittleGroup, false)` with `modw=false`, preserving full (unreduced)
translations. Changed `collect_compatible` from:
```julia
lgirsv = irreps(cbr)
lgs = [primitivize(group(first(lgirs))) for lgirs in lgirsv]
```
to:
```julia
clgirsv = irreps(cbr)
lgirsv = primitivize.(clgirsv)
lgs = group.(lgirsv)
```

**Result: all 230 SGs pass, all dimensions, all EBRs. Zero failures.**

## Summary of error classes and fixes

The symmetry analysis bug (PR #89) had three distinct error classes:

1. **Θ_G sign** — Used `Θ_G` instead of `conj(Θ_G)` to match `calc_bandreps`' trace
   convention. Fix: use `-G` in `reciprocal_translation_phase`. Resolved all 2D failures.

2. **Global phase** — The `SiteInducedSGRepElement` functor computes `e^{-2πi(gk)·v}`
   (physical Convention 1), but `calc_bandreps` uses `e^{+2πi(gk)·v}` (Crystalline.jl
   issue #12). Fix: multiply by `cispi(4dot(gk, v))` in `symmetry_eigenvalues`. Resolved
   3D failures at k-points with non-symmorphic operations.

3. **Translation reduction under primitivization** — `primitivize(::LittleGroup)` with
   default `modw=true` discards lattice vectors from translations, corrupting phase factors
   `e^{2πi(gk)·v}` at non-Γ k-points in centered lattices. Fix: use
   `primitivize(::Collection{LGIrrep})` which passes `modw=false`. Resolved all remaining
   centered-lattice failures (SG 68, 88, 141, 142, 214, 220, 230).
