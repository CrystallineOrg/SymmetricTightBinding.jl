# Changelog

## v0.2.0 (unreleased)

### Breaking changes

- `sgrep_induced_by_siteir` is renamed to `site_induced_sgrep` (and the internal
  `sgrep_induced_by_siteir_excl_phase` to `site_induced_sgrep_excl_phase`) (#140).
- `obtain_symmetry_related_hoppings` is no longer exported.
- The orbitals of each band representation are now ordered partner-function-major (i.e.,
  sites iterate fastest; then site partner functions; then BRs) rather than site-major, so
  that spinful models take a pseudospin-block form (#141).
  This changes the row/column order of `H(k)`, eigenvectors, `orbital_positions`,
  `site_induced_sgrep` matrices, and gradient matrices for band representations with both
  multiple sites and multidimensional site irreps; the relation between previously
  established tight-binding model terms and their coefficient vectors `cs` is expected
  to be unaffected.
- The `Mm` field of `TightBindingBlock` (the M-matrix encoding of a block) is now a
  `Matrix{SparseMatrixCSC{Int, Int}}`, indexed as `Mm[s,t][i,j]`, rather than a dense
  `Array{Int, 4}` indexed as `Mm[i,j,s,t]`. This substantially speeds up constraint
  assembly and hence also `tb_hamiltonian` while also reducing assembly allocations for
  large models.
