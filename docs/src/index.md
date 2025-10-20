# SymmetricTightBinding.jl

SymmetricTightBinding.jl provides tools for the construction and manipulation of tight-binding models. The main novelty -- and principal strength -- of the package is that every such model is associated with, and specified by, a set of [band representations](https://doi.org/10.1146/annurev-conmatphys-041720-124134).

Put more simply, SymmetricTightBinding.jl will automatically generate all possible tight-binding Hamiltonians that are compatible with a global (space group) symmetry, as well as a selection of orbitals with specified local symmetries (i.e., transforming as specific site symmetry irreps), each situated at specified positions in the unit cell (i.e., at specific Wyckoff positions).

The underlying physics is that the Bloch Hamiltonian of a Wannierizable set of bands must transform under under a site-symmetry induced representation (also called band representation) $D(g)$ for operations $g$ in the associated space group. That is, the Bloch Hamiltonian $\mathbf{h}(\mathbf{k})$ must be symmetric in the sense:

```math
\mathbf{D}_\mathbf{k}(g) \mathbf{h}_{\mathbf{k}} \mathbf{D}_\mathbf{k}^\dagger(g) = \mathbf{h}_{g\mathbf{k}}
```

where $\mathbf{D}_{\mathbf{k}}(g)$ is the (momentum-)block-diagonal part of the Fourier transformed band representation $D(g)$. This is simply the familiar operator relation $g \hat{h}_{\mathbf{k}} g^{-1} = \hat{h}_{g\mathbf{k}}$ cast into the basis of local orbitals.

## Installation

The package is registered in the Julia General registry and can be installed from the `pkg>` command line (entered by pressing `]` in the Julia REPL):

```jl
pkg> add SymmetricTightBinding, Crystalline
```

SymmetricTightBinding.jl is designed to work as a companion package to [Crystalline.jl](https://github.com/thchr/Crystalline.jl); so we add Crystalline.jl in the above as well.
The packages can subsequently be loaded at the Julia REPL:
```@repl
julia> using SymmetricTightBinding, Crystalline
```

## Table of contents

```@contents
Pages = ["index.md", "tutorial.md", "band-symmetry.md", "symmetry-breaking.md", "api.md", "internal-api.md"]
Depth = 2
```