# Developer notes overview

These notes are for developers working on the SymmetricTightBinding.jl internals. For the
user-facing theory exposition, see [`theory.md`](../theory.md).

## Files

- **[`trs_notes.md`](trs_notes.md)** — Co-representation theory (following Wooten Ch. 12),
  grey group specialization, realification, and quantization of TRS. See also `theory.md`
  for the polished user-facing treatment of symmetry constraints and time-reversal.

- **[`fourier.md`](fourier.md)** — Quick reference for Convention 1 (PythTB-style, position
  phase included) vs Convention 2 (lattice phase only), including second-quantized notation.
  See also `theory.md` Appendix A for a more thorough treatment.

- **[`1d_example.md`](1d_example.md)** — Worked example of a 1D bipartite lattice with
  inversion, illustrating the M-matrix approach by comparing a hand-derived Hamiltonian with
  the symmetry-constrained result.

- **[`symmetry_eigenvalue_conventions.md`](symmetry_eigenvalue_conventions.md)** — Phase
  convention mismatch between the physical derivation (`theory.md`) and Crystalline.jl's
  `calc_bandreps`, how `symmetry_eigenvalues` corrects for it, the centered-lattice
  primitivization fix, and options for future cleanup.
