# SymmetricTightBinding.jl — Improvement plan

A plan for refactoring, testing, documentation, and bug-fixing work on
SymmetricTightBinding.jl. Work is a collaboration between Thomas Christensen and Antonio
Morales. Each phase is worked in a separate branch and PR'd for review.

Status key: `[ ]` = not started, `[~]` = in progress, `[x]` = done

---

## Phase 1: CLAUDE.md `[x]`

**Branch:** `main` (committed directly)

Created `CLAUDE.md` summarizing API, project structure, conventions, known issues, and
developer preferences.

---

## Phase 2: TODO audit `[x]`

**Branch:** `todo-audit`

Catalog all TODOs in the codebase, distinguish big-picture items from minor nits, and
document them. Remove or resolve any that are trivially fixable.

### TODOs found

| Location | Text | Category |
|---|---|---|
| `src/types.jl:223` | Could use knowledge that `-delta` is in orbit | Minor optimization |
| `src/types.jl:455` | ANTIHERMITIAN model `solve` not implemented | **Feature gap** |
| `src/symmetry_analysis.jl:112,117` | Preallocate Theta_G; avoid allocations | Performance nit |
| `ext/MakieExt.jl:50` | Actually do Observable updates in recipe | Feature gap |
| `ext/MakieExt.jl:336` | How to pass kwargs with Makie SpecApi | Minor API gap |
| `ext/OptimExt.jl:142` | Improve initial guess in fitting | Enhancement |
| `test/sg_tb_hamiltonian.jl:6` | Implement 3D space group tests | **Major test gap** |
| `test/pg_tb_hamiltonian.jl:51,56` | Finish structure test; add more PG cases | **Test gap** |

### Big-picture items to track

1. ANTIHERMITIAN `solve` (and downstream spectrum/analysis support)
2. 3D space group Hamiltonian tests (file is empty)
3. More plane group test cases beyond graphene
4. Symmetry analysis correctness (Phase 6)
5. Performance for high SG numbers (issue #44)

### Outcome

All TODOs catalogued above. None warrant immediate action — they are either minor
optimizations, feature gaps that need design work (e.g., ANTIHERMITIAN support), or test
gaps addressed in Phase 4. One stale TODO removed from `test/pg_tb_hamiltonian.jl:51`
(referenced structure that has since been implemented). The rest remain in code as-is.

---

## Phase 3: Devdocs restructuring `[x]`

**Branch:** `devdocs-restructure`

### Current state

- `devdocs/devdocs.md` — early notes; photonic framing; strict subset of `docs.md`
- `devdocs/docs.md` — most complete derivation; photonic framing; incomplete hermiticity section
- `devdocs/fourier.md` — short, focused; Convention 1 vs 2
- `devdocs/trs_notes.md` — co-representation theory; heavy overlap with `docs.md`
- `docs/src/theory.md` — polished user-facing synthesis; no photonic framing

### Plan

Move devdocs into the Documenter tree at `docs/src/devdocs/` so they live alongside the
user-facing docs but remain developer-oriented.

**Step 1 — Tidy and consolidate:**
- Delete `devdocs/devdocs.md` (fully superseded by `docs.md`)
- Move remaining files to `docs/src/devdocs/`
- Rename `docs.md` -> `derivations.md`; strip photonic-specific framing; fill in the
  incomplete hermiticity section
- Merge `trs_notes.md` content into `derivations.md` as a co-representations section
  (the unique content `trs_notes.md` adds is the Wigner/Frobenius classification and
  the W = V * Lambda^{1/2} * V^T trick for physically real forms)
- Keep `fourier.md` as-is (it's focused and self-contained)
- Add a short `docs/src/devdocs/README.md` explaining what each file covers

**Step 2 — Cross-link:**
- Add references from `docs/src/theory.md` to the devdocs for readers wanting more detail
- Ensure `docs/make.jl` includes the devdocs pages (if desired in built docs) or excludes
  them intentionally

**Principle:** `docs/src/theory.md` = interested-user-facing exposition;
`docs/src/devdocs/` = developer reference for conventions, derivations, and physics details.

### Outcome

Devdocs restructured in 6 commits (3a–3f):
- Deleted `devdocs/devdocs.md` (superseded); removed photonic framing from remaining files
- Extracted 1D bipartite lattice example into standalone `1d_example.md`
- Deleted old `devdocs/docs.md`, transferring unique content to `trs_notes.md` and `theory.md`
- Moved `devdocs/` to `docs/src/devdocs/`; added `README.md` index
- Fixed typos, grammar, and math errors (sign error in theory.md Fourier derivation,
  missing 1/√N normalizations)
- Polished prose for public readability: neutral voice, removed intimate discussion notes,
  fixed non-standard constructions

---

## Phase 4: Tests & coverage `[x]`

**Branch:** `tests/coverage-improvements`

Starting coverage: 55%, 116 tests. One branch for all test improvements.

### 4A: `spectrum` tests (easy win)
- `spectrum(ptbm, ks)` returns correct shape (Nk x N)
- Eigenvalues at Gamma match direct `eigvals(Hermitian(ptbm([0,0])))`
- Use existing graphene model from `pg_tb_hamiltonian.jl`

### 4B: `show` method regression tests
- Full output regression tests for `HoppingOrbit`, `TightBindingTerm`,
  `TightBindingModel`, `ParameterizedTightBindingModel`, `TightBindingElementString`
- Uses DeepDiffs.jl for readable failure diffs (pattern from Crystalline.jl's test/show.jl)
- Tests `TightBindingElementString` ANSI color output for active/inactive/zero states

### 4C: Symmetry analysis tests
- Add the two graphene examples from PR #89 (plane group 17) as proper `@test` assertions
  — these are known to pass
- Run the full EBR scan from PR #89 offline to identify which cases currently pass; add
  those as `@test`
- Pick a curated subset of failing cases (preferring simpler/lower-dimensional groups;
  representatives from each failure class) and add as `@test_broken`
- Uncomment `include("symmetry_analysis.jl")` in `runtests.jl`

### 4D: Docs examples as tests
- `docs/src/tutorial.md` graphene example: add spectrum evaluation check (partially covered
  already in `pg_tb_hamiltonian.jl`)
- `docs/src/band-symmetry.md`: covered by 4C above
- `docs/src/symmetry-breaking.md`: already covered in `test/symmetry-breaking.jl`
- `docs/src/berry.md`: already well-covered in `test/berry.jl`

### 4E: Extension tests (requires extra deps in test env)
- `test/ext_optim.jl`: test `fit` on Haldane model with known target coefficients; verify
  recovery within tolerance
- Makie extension: smoke test only (call `plot`, check no error); may skip if CI display
  setup is too annoying

### 4F: More plane/space group Hamiltonian tests
- Add 2-3 more plane group cases to `pg_tb_hamiltonian.jl` (e.g., p4mm, p3)
- Add at least one 3D space group case to `sg_tb_hamiltonian.jl` (e.g., SG 221 / Pm-3m,
  which already has handwritten reps in `test/site_representations.jl`)

### Outcome

Tests expanded from 116 to 293 (281 passing, plus 12 expected `@test_broken`).
Coverage improved from 55% to 95% (measured on `src/` files). New files and changes:

- **`test/spectrum.jl`** (4A): 18 tests covering single/multi-k evaluation, shape,
  correctness vs direct `eigvals`, transform keyword, Dirac degeneracy at K
- **`test/show.jl`** (4B): 13 regression tests for `show`/`summary` of all core types
  (`HoppingOrbit`, `TightBindingTerm`, `TightBindingModel`,
  `ParameterizedTightBindingModel`, `TightBindingElementString`), using DeepDiffs.jl for
  readable failure output (pattern from Crystalline.jl)
- **`test/symmetry_analysis_stopgap.jl`** (4C): 114 tests (104 pass, 10 broken) covering
  all 2D plane groups with special Wyckoff positions. Interim/stopgap tests until PR #89 is
  resolved. Systematic scan identified failures at SG 13 (1b/1c), SG 14 (1b/1c), and
  SG 16 (3c) — all involving K-point 3-fold symmetry phase issue. SG 16 3c EBRs are flaky
  (coefficient-dependent) so are skipped rather than marked broken.
- **`test/gradients.jl`** (new, not in original plan): 18 tests for
  `gradient_wrt_hopping`, `gradient_wrt_momentum`, `energy_gradient_wrt_hopping`, including
  finite-difference validation and degenerate-band handling
- **`test/pg_tb_hamiltonian.jl`** (4F): added p4mm (C₄ spectrum symmetry), p3 (Hermiticity),
  p2 (TRS H(k)=H*(-k)); removed stale TODO, added Hermiticity checks to graphene
- **`test/sg_tb_hamiltonian.jl`** (4F): populated with SG 2 (P-1, Hermiticity), SG 16
  (P222, TRS), SG 225 (Fm-3m), 1D SG 2 (periodicity), SG 47 (Pmmm, multi-EBR composite)
- **`test/runtests.jl`**: added includes for spectrum, show, gradients, symmetry_analysis_stopgap

Deferred to future work: 4E (extension tests — requires adding Optim/Makie to test deps)

---

## Phase 5: Maintenance & refactoring `[ ]`

**Branch:** `refactor/cleanup`

Done *after* Phase 4 so tests catch regressions. Incremental, file-by-file review.

### Areas to examine

- **Code duplication in `tightbinding.jl`:** the constraint matrix assembly functions
  (`representation_constraint_matrices`, `reciprocal_constraints_matrices`,
  `constraint_hermiticity`) share structural patterns that may be factorable
- **Dead/outdated comments:** several reference early development or photonic context
- **`types.jl` `_getindex` logic:** the pretty-printing path is complex; see if it can be
  simplified without losing functionality
- **`show.jl`:** review for unnecessary complexity
- **Crystalline.jl API check:** verify we're not reimplementing functions that exist upstream
  (orbit handling, translation utilities, etc.)
- **Unused code:** any functions that are defined but never called

### Approach

Read each source file carefully, propose changes, verify tests pass after each change.
Keep refactoring PRs small and reviewable.

---

## Phase 6: Symmetry analysis correctness (PR #89) `[x]`

**Branch:** `explore-fix-hamiltonian-phase-v2` (PR #104, targeting `tests/coverage-improvements`)

### Root cause

Two distinct error classes, all in `src/symmetry_analysis.jl`:

1. **Convention mismatch (net effect: complex conjugation)**: `symmetry_eigenvalues`
   computed the Convention 1 character `χ_C1 = (Θ_G w)† D_k w`, but Crystalline.jl's
   `calc_bandreps` and `lgirreps` use a conjugated sign convention for both the `Θ_G`
   factor and the global phase of `D_k` (cf. issue #12). These two conjugations together
   give `χ_Crystalline = conj(χ_C1)`. Fix: compute `χ_C1` as before, then return
   `conj(χ_C1)`. This resolved all 2D and 3D failures from the phase convention mismatch.

2. **Translation reduction under primitivization**: `primitivize(::LittleGroup)` with default
   `modw=true` discards lattice vectors from translations, corrupting phases at non-Γ
   k-points in centered lattices. Fix: use `primitivize(::Collection{LGIrrep})` which passes
   `modw=false`. Resolved all centered-lattice failures (SG 68, 88, 141, 142, 214, 220, 230).

### Additional changes

- New devdoc `docs/src/devdocs/symmetry_eigenvalue_conventions.md` explaining the mismatch
  between Convention 1 and Crystalline.jl's convention, the implemented fix, and options
  for future cleanup
- `[⚠️ phase]` code annotations at convention-sensitive locations in `symmetry_analysis.jl`
- Symmetry analysis tests rewritten with deterministic RNG and full 230-SG coverage
- Stopgap `@test_broken` for p3/p6 flipped to `@test`

### Note on Hamiltonian phase convention

The symmetry eigenvalue correction is coupled to the Hamiltonian Fourier phase in
`types.jl`. The correct Convention 1 phase is `cispi(-2k·δ)` (effective phase
`e^{+ik·δ}`, with δ = b+R−a). The fix above was derived and tested with this convention.
Changing the Hamiltonian phase sign would require re-deriving the correction in
`symmetry_analysis.jl`.

### Relevant upstream issues
- Crystalline.jl issue #12 (sign convention in `calc_bandreps`)
