# SymmetricTightBinding.jl

## Overview

A Julia package for constructing symmetry-constrained tight-binding Hamiltonians in any
crystallographic space group (1D, 2D, 3D). Built on the framework of topological quantum
chemistry (TQC) and the theory of band representations.

Given a set of elementary band representations (EBRs) from Crystalline.jl and a list of
hopping ranges, the package:
1. Enumerates symmetry-related hopping orbits
2. Constructs the M-matrix encoding all allowed Hamiltonian terms
3. Enforces spatial symmetry, time-reversal, and hermiticity constraints via nullspace +
   Zassenhaus intersection
4. Returns a `TightBindingModel` whose free parameters can be tuned to produce concrete
   Hamiltonians

**Repository:** https://github.com/CrystallineOrg/SymmetricTightBinding.jl
**Depends on:** [Crystalline.jl](https://github.com/thchr/Crystalline.jl) (local copy at `../Crystalline`)
**Authors:** Antonio Morales Perez, Thomas Christensen

## Key types and API

### Core pipeline
- `tb_hamiltonian(cbr, Rs; antihermitian)` — top-level entry point; returns `TightBindingModel`
- `obtain_symmetry_related_hoppings(Rs, br_a, br_b)` — enumerates hopping orbits
- `sgrep_induced_by_siteir(br, op)` — site-symmetry induced space group representation

### Types (defined in `src/types.jl`)
- `HoppingOrbit{D}` — orbit of symmetry-related hopping vectors {delta_i}
- `TightBindingBlock{D}` — single block of the Hamiltonian (br1 -> br2 hopping)
- `TightBindingTerm{D}` — block embedded in full matrix, with hermiticity info
- `TightBindingModel{D}` — collection of terms; functor `tbm(cs)` -> parameterized model
- `ParameterizedTightBindingModel{D}` — model with coefficients; functor `ptbm(k)` -> H(k)

### Analysis tools
- `spectrum(ptbm, ks)` — band eigenvalues over k-path (`src/spectrum.jl`)
- `symmetry_eigenvalues(ptbm, ops, k, sgreps)` — irrep content at high-symmetry k (`src/symmetry_analysis.jl`)
- `collect_compatible(ptbm)` / `collect_irrep_annotations(ptbm)` — band symmetry labels
- `berrycurvature(ptbm, k, n)` / `chern(ptbm, n, Nk)` / `chern_fukui(ptbm, n, Nk)` (`src/berry.jl`)
- `gradient_wrt_hopping` / `gradient_wrt_momentum` (`src/gradients.jl`)
- `subduced_complement(tbm, sgnum_H)` — new terms from symmetry breaking (`src/symmetry_breaking.jl`)

### Utilities
- `pin_free!(brs, idx2abc)` — fix free Wyckoff parameters (`src/utils.jl`)
- `solve(ptbm, k; bloch_phase)` — eigenvalues + eigenvectors (in `src/types.jl`)
- `fit(tbm, E_ref, ks)` — least-squares fitting via Optim.jl extension

### Extensions
- `SymmetricTightBindingMakieExt` — Makie recipes for hopping orbits and band structures
- `SymmetricTightBindingOptimExt` — model fitting via Optim.jl (multi-start LBFGS)

## Fourier convention

The package uses **Convention 1** (PythTB-style) throughout: the Bloch basis includes
position phases, i.e., `exp(ik*(t + q_alpha))`. This differs from much of the literature
(Convention 2, lattice-phase only). See `docs/src/theory.md` Appendix A and
`docs/src/devdocs/fourier.md`. The `solve(...; bloch_phase=Val(true))` option converts eigenvectors
to Convention 2.

## Project structure

```
src/
  SymmetricTightBinding.jl  # module definition, exports, constants
  types.jl                  # core data structures + evaluation functors
  tightbinding.jl           # constraint pipeline: M-matrix, nullspace, sparsification
  site_representations.jl   # coset-decomposition-based induced representations
  symmetry_analysis.jl      # irrep content at high-symmetry k-points
  spectrum.jl               # band eigenvalue computation
  berry.jl                  # Berry curvature, Chern numbers (Kubo + Fukui)
  gradients.jl              # dH/dc_i and dH/dk_i
  hermiticity.jl            # hermiticity/anti-hermiticity constraint intersection
  timereversal.jl           # TRS constraint: H(k) = H*(-k)
  zassenhaus.jl             # Zassenhaus algorithm for subspace intersection
  symmetry_breaking.jl      # subduced complement for translationengleiche subgroups
  utils.jl                  # orbital positions, pin_free!, split_complex, etc.
  show.jl                   # pretty-printing for all types
ext/
  SymmetricTightBindingMakieExt.jl   # Makie visualization recipes
  SymmetricTightBindingOptimExt.jl   # Optim.jl fitting
test/
  runtests.jl               # test runner
  pg_tb_hamiltonian.jl      # plane group tests (graphene, p4mm, p3, p2)
  sg_tb_hamiltonian.jl      # space group tests (SG 2, 16, 47, 225; 1D SG 2)
  site_representations.jl   # representation matrix tests
  symmetry-breaking.jl      # subduced complement tests
  symmetry_analysis.jl      # comprehensive symmetry tests (commented out; see PR #89)
  symmetry_analysis_stopgap.jl # interim symmetry vector tests (see PR #89)
  berry.jl                  # Berry curvature + Chern number tests
  spectrum.jl               # band eigenvalue tests
  show.jl                   # regression tests for show/summary methods
  gradients.jl              # hopping/momentum gradient tests
docs/src/
  tutorial.md               # graphene walkthrough
  theory.md                 # mathematical framework (polished)
  band-symmetry.md          # symmetry analysis example
  symmetry-breaking.md      # symmetry reduction example
  berry.md                  # Haldane model, Chern numbers
  devdocs/
    README.md               # index of developer docs
    trs_notes.md            # co-representation theory, TRS quantization, realification
    fourier.md              # Convention 1 vs 2 quick reference
    1d_example.md           # worked 1D bipartite lattice example
```

## Running tests

```bash
cd /path/to/SymmetricTightBinding
julia --project -e 'using Pkg; Pkg.test()'
```

Or interactively:
```julia
using Pkg; Pkg.activate(".")
include("test/runtests.jl")
```

Tests take ~2 minutes.

## Known issues

1. **ANTIHERMITIAN `solve` not implemented** (`src/types.jl:455`).

2. **Performance** (issue #44): `tb_hamiltonian` is slow for high space group numbers
   (SG ~200+ can take minutes per EBR).

3. **Test coverage gaps:** No tests for extensions (Optim, Makie).

4. **Symmetry eigenvalue convention mismatch (workaround in place):**
   `symmetry_eigenvalues` monkey-patches two phase corrections to match Crystalline.jl's
   `calc_bandreps` convention (Crystalline.jl issue #12). Code locations are marked with
   `[⚠️ phase]`. See `docs/src/devdocs/symmetry_eigenvalue_conventions.md` for details and
   recommended future cleanup (Option C: change `SiteInducedSGRepElement` convention).

## Improvement plan

See `PLAN.md` for the phased work plan. Completed: Phase 1 (CLAUDE.md), Phase 2 (TODO
audit), Phase 3 (devdocs restructuring), Phase 4 (tests/coverage), Phase 6 (symmetry
analysis fix). Remaining: Phase 5 (refactoring).

## Code conventions and preferences

- Follow existing docstring style (see e.g. `HoppingOrbit`, `TightBindingBlock` in `types.jl`)
- Use Unicode math where it improves readability (e.g., `δ`, `ρ`, `Θ`)
- Numeric tolerances are defined as module-level constants (`NULLSPACE_ATOL_DEFAULT`, etc.)
- The codebase uses Convention 1 Fourier transform everywhere; be careful with phases
- Prefer using Crystalline.jl's API over reimplementing group-theoretic operations
- Don't add new functionality without tests
- Always run tests after making changes
- Keep comments brief but include enough detail for non-obvious physics/math
- Work in separate branches for each groupable change (for PR review with collaborators)

## Gotchas

- **`ptbm(k)` returns mutable scratch:** `ParameterizedTightBindingModel` reuses an internal
  buffer. When comparing evaluations at different k-points, `copy()` the result.
- **`@composite` macro hygiene:** The `@composite` macro evaluates in Crystalline's module
  scope, so loop variables are not visible. Use `CompositeBandRep(coefs, brs)` constructor
  in programmatic contexts (loops, comprehensions).
- **`TightBindingElementString` is not exported:** import explicitly if needed in tests via
  `using SymmetricTightBinding: TightBindingElementString`.
