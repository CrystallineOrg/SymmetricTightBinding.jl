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

## Phase 3: Devdocs restructuring `[ ]`

**Branch:** `docs/devdocs-restructure`

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

---

## Phase 4: Tests & coverage `[ ]`

**Branch:** `tests/coverage-improvements`

Current coverage: 55%. One branch for all test improvements.

### 4A: `spectrum` tests (easy win)
- `spectrum(ptbm, ks)` returns correct shape (Nk x N)
- Eigenvalues at Gamma match direct `eigvals(Hermitian(ptbm([0,0])))`
- Use existing graphene model from `pg_tb_hamiltonian.jl`

### 4B: `show` method tests
- `sprint(show, obj)` for `HoppingOrbit`, `TightBindingBlock`, `TightBindingTerm`,
  `TightBindingModel`, `ParameterizedTightBindingModel`
- Verify no errors and basic content presence (e.g., check for expected substrings)

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

## Phase 6: Symmetry analysis correctness (PR #89) `[ ]`

**Branch:** new branch (or continue from `symeigs-test`)

The hardest problem. `symmetry_eigenvalues` returns incorrect irrep labels at many
high-symmetry k-points, especially K/KA/H/HA (3-fold symmetry).

### Root cause hypothesis

The Theta_G phase factor `exp(-2*pi*i * G * q_i)` in `symmetry_eigenvalues` has a
sign/convention mismatch. Conjugating the phase fixes K-point failures but introduces new
ones elsewhere — so it's not a simple sign flip.

### Approach

**Step 1 — Anchor on a minimal test case:**
- Plane group 13, EBR `(1c|A)`, on-site only (zero hopping range)
- All space group representations are the identity, so the only source of non-trivial
  irreps is Theta_G
- Hand-compute the expected symmetry eigenvalues and compare with `symmetry_eigenvalues`

**Step 2 — Trace the phase pipeline:**
- Convention 1 Bloch basis definition
- `sgrep_induced_by_siteir` output
- `symmetry_eigenvalues` Theta_G application
- Crystalline.jl's `calc_bandreps` convention (and its issue #12)

**Step 3 — Classify failures:**
- Which k-points are affected?
- Is the pattern consistent (e.g., always complex-conjugated)?
- Are there cases where the *structure* of the symmetry vector is wrong (not just permuted
  labels)?

**Step 4 — Fix and validate:**
- Propose a fix (may involve changes in both SymmetricTightBinding and Crystalline)
- Run the full EBR scan to verify; convert `@test_broken` from Phase 4C to `@test`

### Relevant upstream issues
- Crystalline.jl issue #12 (sign convention in `calc_bandreps`)
- Possibly other Crystalline.jl issues identified during Step 2
