"""
    Wannier90TightBindingModel{D}

The real-space hopping matrices ``H(\\mathbf{R})`` of a parameterized tight-binding model,
in the layout of Wannier90's
[`_hr.dat`](https://wannier90.readthedocs.io/en/latest/user_guide/wannier90/files/#seedname_hrdat)
files: containing one matrix ``H(\\mathbf{R})`` per lattice translation ``\\mathbf{R}``,
with elements ``H_{ij}(\\mathbf{R}) = ⟨i,\\mathbf{0}|\\hat{H}|j,\\mathbf{R}⟩``, i.e., the
amplitude for hopping from orbital `j` in the unit cell at ``\\mathbf{R}`` to orbital `i` in
the home unit cell.

Construct from a parameterized model `ptbm` via `Wannier90TightBindingModel(ptbm)`; write
to disk via [`write_wannier90`](@ref).

## Fields
- `translations :: Vector{SVector{D, Int}}`: featured lattice translations ``\\mathbf{R}``,
  in the (primitive) lattice basis, sorted lexicographically.
- `hopping_matrices :: Vector{Matrix{ComplexF64}}`: `hopping_matrices[n]` is
  ``H(\\mathbf{R})`` for ``\\mathbf{R}`` = `translations[n]`.
- `positions :: Vector{DirectPoint{D}}`: orbital positions (Wannier centres), in the
  (primitive) lattice basis.

## Evaluation
A `Wannier90TightBindingModel` `w90tbm` can be evaluated at a momentum `k` as `w90tbm(k)`,
returning the Bloch Hamiltonian in SymmetricTightBinding's Fourier convention (Convention
1), i.e., ``H_{ij}(\\mathbf{k}) = \\sum_{\\mathbf{R}} H_{ij}(\\mathbf{R})
\\mathrm{e}^{2π\\mathrm{i}\\mathbf{k}\\cdot(\\mathbf{R} + \\mathbf{r}_j - \\mathbf{r}_i)}``,
same as `ptbm(k)`.

## Fourier conventions
The hopping matrices ``H(\\mathbf{R})`` are real-space quantities and do not depend on a
Fourier convention. The associated full Hamiltonian is:

```math
\\hat{H} = \\sum_{ij} \\sum_{\\mathbf{R}} H_{ij}(\\mathbf{R}) \\sum_{\\mathbf{R}'}
    \\hat{c}^†_{i,\\mathbf{R}'} \\hat{c}_{j,\\mathbf{R}'+\\mathbf{R}},
```

with ``\\hat{c}^†_{i,\\mathbf{R}}`` creating a particle in orbital `i` of the unit cell at
``\\mathbf{R}``. A Fourier convention does enter when transforming this to momentum space,
i.e., for evaluating the associated Bloch Hamiltonian ``\\mathbf{H}(\\mathbf{k})`` at
momentum ``\\mathbf{k}``.
SymmetricTightBinding uses Convention 1, placing orbitals at the Wannier centres
``\\mathbf{r}_i`` of `_centres.xyz`, as described under "Evaluation" above.
Wannier90 itself transforms with lattice phases only (Convention 2):
``H_{ij}^{(2)}(\\mathbf{k}) = \\sum_{\\mathbf{R}} H_{ij}(\\mathbf{R})
\\mathrm{e}^{2π\\mathrm{i}\\mathbf{k}\\cdot\\mathbf{R}}``.
The two conventions are related by a diagonal unitary ``Θ_{\\mathbf{k}} =
\\mathrm{diag}(\\mathrm{e}^{-2π\\mathrm{i}\\mathbf{k}\\cdot\\mathbf{r}_i})`` such that
``H^{(1)}(\\mathbf{k}) = Θ_{\\mathbf{k}} H^{(2)}(\\mathbf{k}) Θ_{\\mathbf{k}}^†``.
Which convention a reader of `_hr.dat` adopts is a choice: e.g., PythTB uses Convention 1
and TBmodels offers both.
"""
struct Wannier90TightBindingModel{D}
    translations::Vector{SVector{D, Int}}
    hopping_matrices::Vector{Matrix{ComplexF64}}
    positions::Vector{DirectPoint{D}}
end
orbital_positions(w90tbm::Wannier90TightBindingModel) = w90tbm.positions
orbital_count(w90tbm::Wannier90TightBindingModel) = length(w90tbm.positions)
Crystalline.dim(::Wannier90TightBindingModel{D}) where D = D

function Wannier90TightBindingModel(
    aptbm::AbstractParameterizedTightBindingModel{D}
) where D
    positions = orbital_positions(aptbm)
    acc = Dict{SVector{D, Int}, Matrix{ComplexF64}}()
    _add_hopping_matrices!(acc, aptbm, positions)
    translations = sort!(collect(keys(acc)))
    return Wannier90TightBindingModel{D}(translations, [acc[R] for R in translations],
                                         positions)
end

function _add_hopping_matrices!(
    acc::Dict{SVector{D, Int}, Matrix{ComplexF64}},
    ptbm::ParameterizedTightBindingModel{D},
    positions::AbstractVector{DirectPoint{D}}
) where D
    N = length(positions)
    for (tbt, c) in zip(ptbm.tbm.terms, ptbm.cs)
        _add_hopping_matrices!(acc, tbt, c, positions, N)
    end
    return acc
end
function _add_hopping_matrices!(
    acc::Dict{SVector{D, Int}, Matrix{ComplexF64}},
    pctbm::ParameterizedCompositeTightBindingModel{D},
    positions::AbstractVector{DirectPoint{D}}
) where D
    # `cs` holds the coefficients of the Hermitian terms, then of the anti-Hermitian terms;
    # as in `pctbm(k)`, the two loops are split to keep each type-stable
    N = length(positions)
    tbm_h = pctbm.tbm.tbm_h # iterates `tbt :: TightBindingTerm{…, HERMITIAN, …}`
    for (tbt, c) in zip(tbm_h.terms, @view pctbm.cs[1:length(tbm_h)])
        _add_hopping_matrices!(acc, tbt, c, positions, N)
    end
    tbm_a = pctbm.tbm.tbm_a # iterates `tbt :: TightBindingTerm{…, ANTIHERMITIAN, …}`
    for (tbt, c) in zip(tbm_a.terms, @view pctbm.cs[length(tbm_h)+1:end])
        _add_hopping_matrices!(acc, tbt, c, positions, N)
    end
    return acc
end

# add the hopping matrices of `c * tbt` into `acc`, mirroring
# `evaluate_tight_binding_term!`: there, each `(n, i, j)` element contributes
# `MmtC[n,i,j] e^{-2πik·δₙ}` to `H[i,j](k)`, with `δₙ = rᵢ - rⱼ - R`; i.e., it is the
# hopping from orbital `j` in cell `R` to orbital `i` in the home cell
function _add_hopping_matrices!(
    acc::Dict{SVector{D, Int}, Matrix{ComplexF64}},
    tbt::TightBindingTerm{D, S},
    c::Real,
    positions::AbstractVector{DirectPoint{D}},
    N::Int
) where {D, S}
    iszero(c) && return acc
    block = tbt.block
    block_i, block_j = tbt.block_ij
    is = tbt.axis[Block(block_i)] # global row indices
    js = tbt.axis[Block(block_j)] # global col indices
    δs = constant.(orbit(block.h_orbit)) # NB: assumes no free parameters in `δ`
    _init_matrix_f = let N = N; () -> zeros(ComplexF64, N, N); end
    for (local_j, j) in enumerate(js)
        for (local_i, i) in enumerate(is)
            for n in eachindex(δs)
                v = c * block.MmtC[n, local_i, local_j]
                abs(v) < PRUNE_ATOL_DEFAULT && continue
                R = _integer_translation(positions[i] - positions[j] - δs[n])
                get!(_init_matrix_f, acc, R)[i, j] += v
                if S !== NONHERMITIAN && block_i ≠ block_j # hermiticity-related block
                    # NB: `H[j,i](k)` gets `±conj(v) e^{+2πik·δₙ}`, i.e., a hopping by `-R`
                    v′ = S === ANTIHERMITIAN ? -conj(v) : conj(v)
                    get!(_init_matrix_f, acc, -R)[j, i] += v′
                end
            end
        end
    end
    return acc
end

function _integer_translation(R::StaticVector{D, <:Real}) where D
    R′ = round.(Int, R)
    if !isapprox(R, R′; atol = VEC_CMP_ATOL)
        error(lazy"hopping implies a non-integer lattice translation $R; orbital positions and hopping vectors are inconsistent")
    end
    return SVector{D, Int}(R′)
end

function (w90tbm::Wannier90TightBindingModel{D})(k::ReciprocalPointLike{D}) where D
    length(k) == D || error(lazy"momentum `k` must be a $D-dimensional vector")
    N = orbital_count(w90tbm)
    H = zeros(ComplexF64, N, N)
    for (R, Hᴿ) in zip(w90tbm.translations, w90tbm.hopping_matrices)
        phase = cispi(2 * dot(k, R))
        H .+= phase .* Hᴿ
    end
    # `H` is now Wannier90's lattice-phase-only Convention 2; add the orbital-position
    # phases to return the package's Convention 1, i.e., H[i,j] *= e^{2πik·(rⱼ - rᵢ)}
    Θ = reciprocal_translation_phase(orbital_positions(w90tbm), k)
    return Θ * H * Θ'
end

# ---------------------------------------------------------------------------------------- #

"""
    write_wannier90(prefix::AbstractString, ptbm; lattice = nothing, digits = 12)
                                                                --> NamedTuple

Write the parameterized tight-binding model `ptbm` to the
[`_hr.dat`](https://wannier90.readthedocs.io/en/latest/user_guide/wannier90/files/#seedname_hrdat)
file format of [Wannier90](https://wannier.org), to enable loading by other packages (e.g.,
[PythTB](https://www.physics.rutgers.edu/pythtb/), [TBmodels](https://tbmodels.greschd.ch),
[WannierTools](https://github.com/quanshengwu/wannier_tools),
[WannierBerri](https://wannier-berri.org),
[WannierNLQG.jl](https://github.com/ZhuochengLu/WannierNLQG)).

Writes the following files:
- `<prefix>_hr.dat`: the real-space hopping matrices ``H(\\mathbf{R})``; always written.
- `<prefix>_centres.xyz` and `<prefix>.win`: the orbital positions (as Wannier centres)
  and the lattice vectors, both in Cartesian coordinates, in the
  [`_centres.xyz`](https://wannier90.readthedocs.io/en/latest/user_guide/wannier90/files/#seedname_centresxyz)
  and [`.win`](https://wannier90.readthedocs.io/en/latest/user_guide/wannier90/files/#seednamewin)
  formats. Written only if `lattice !== nothing`. Some readers such as PythTB need them
  alongside `_hr.dat`.

Returns the paths of the written files as a `NamedTuple` of `String`s with keys `path_hr`,
`path_centres`, and `path_win` (files not written have empty paths).

`ptbm` is an [`AbstractParameterizedTightBindingModel`](@ref) (or an internal
[`Wannier90TightBindingModel`](@ref)). `prefix` may include a directory.

## Keyword arguments
- `lattice`: the *primitive* direct basis (`DirectBasis{D}`) that the model's positions
  and translations refer to; for a conventional basis `Rs` of space group `sgnum`, this is
  `primitivize(Rs, centering(sgnum, D))`. Units are implicitly assumed to be Å.
- `digits`: decimal digits of the matrix elements in `_hr.dat`.

## Conventions and limitations
- `_hr.dat` lists ``H_{ij}(\\mathbf{R}) = ⟨i,\\mathbf{0}|\\hat{H}|j,\\mathbf{R}⟩``.
  Wannier90 itself constructs Bloch Hamiltonians from this in Convention 2 (lattice phases
  only), i.e., ``H_{ij}(\\mathbf{k}) = \\sum_{\\mathbf{R}} H_{ij}(\\mathbf{R})
  \\mathrm{e}^{2π\\mathrm{i}\\mathbf{k}\\cdot\\mathbf{R}}``. Not all readers of `_hr.dat`
  files adopt Convention 2, however: e.g., PythTB adopts Convention 1 (orbital phases
  included), consistent with SymmetricTightBinding's choice.
- All Wigner–Seitz degeneracy weights are 1: each ``\\mathbf{R}`` appears only once.
- Wannier90's formats are three-dimensional. For `D < 3`, translations are padded with
  zeros, and `lattice` is extended by unit vectors along the missing Cartesian axes, i.e.,
  the model is placed in the ``xy``-plane (`D = 2`) or along the ``x``-axis (`D = 1`).
- Matrix elements are written in fixed-point notation, as Wannier90 does, with `digits`
  decimals; elements smaller than `10^(-digits)` are lost.
"""
function write_wannier90(
    prefix::AbstractString,
    w90tbm::Wannier90TightBindingModel{D};
    lattice::Union{Nothing, DirectBasis{D}} = nothing,
    digits::Integer = 12,
) where D
    # header with a timestamp in Wannier90's own format, e.g., "… on 04Nov2026 at 10:12:33"
    date = Dates.format(Dates.now(), "dduyyyy at HH:MM:SS")
    header = "Created by SymmetricTightBinding.jl on " * date # shared by all written files

    # write the `_hr.dat` file
    path_hr = prefix * "_hr.dat"
    open(io -> _write_wannier90_hr(io, w90tbm, header, digits), path_hr, "w")

    # optionally write `_centres.xyz` and `.win` files, if a lattice is provided
    paths = if !isnothing(lattice)
        Rm = _embed_lattice_3d(lattice)
        open(prefix * "_centres.xyz", "w") do io
            _write_wannier90_centres(io, w90tbm, header, Rm)
        end
        open(io -> _write_wannier90_win(io, w90tbm, header, Rm), prefix * ".win", "w")
        (; path_hr, path_centres = prefix * "_centres.xyz", path_win = prefix * ".win")
    else
        (; path_hr, path_centres = "", path_win = "")
    end

    return paths
end
function write_wannier90(
    prefix::AbstractString,
    ptbm::AbstractParameterizedTightBindingModel;
    kws...
)
    w90tbm = Wannier90TightBindingModel(ptbm)
    return write_wannier90(prefix, w90tbm; kws...)
end

function _pad_to_3d_vec(v::AbstractVector{T}) where T
    return SVector{3, T}(ntuple(i -> i ≤ length(v) ? v[i] : zero(T), Val(3)))
end

# the lattice vectors as columns of a 3×3 matrix (cf. `stack(lattice)`), extending a `D < 3`
# lattice by unit vectors along the missing Cartesian axes
function _embed_lattice_3d(lattice::DirectBasis{D}) where D
    Rm = MMatrix{3, 3, Float64, 9}(I)
    Rm[1:D, 1:D] .= stack(lattice)
    return Rm
end

function _write_wannier90_hr(
    io::IO, w90tbm::Wannier90TightBindingModel, header::String, digits::Integer
)
    N = orbital_count(w90tbm)
    Nᴿ = length(w90tbm.translations)
    # 1st line: header
    println(io, header)
    # 2nd line: Hamiltonian dimension (orbital count)
    println(io, N)
    # 3rd line: number of translations (hopping matrices H(R))
    println(io, Nᴿ)
    # 4th line and below: degeneracy weights (1 per translation, 15 per line): all unity
    for chunk in Iterators.partition(1:Nᴿ, 15)
        println(io, join(fill("1", length(chunk)), " "))
    end
    # subsequent lines: hopping matrices H(R) given as R₁ R₂ R₃ i j Re(Hᵢⱼ(R)) Im(Hᵢⱼ(R))
    # one line per matrix element Hᵢⱼ(R), with i iterating fastest
    fmt = Printf.Format("%3d %3d %3d %3d %3d % .$(digits)f % .$(digits)f\n")
    for (R, H) in zip(w90tbm.translations, w90tbm.hopping_matrices)
        R³ = _pad_to_3d_vec(R)
        for j in 1:N, i in 1:N # row index fastest, as Wannier90
            Printf.format(io, fmt, R³[1], R³[2], R³[3], i, j, reim(H[i, j])...)
        end
    end
end

function _write_wannier90_centres(
    io::IO, w90tbm::Wannier90TightBindingModel, header::String, Rm::AbstractMatrix
)
    println(io, orbital_count(w90tbm))
    println(io, "Wannier centres. ", header)
    fmt = Printf.Format("X % .10f % .10f % .10f\n") # `X`: standard Wannier centre/dummy atom symbol
    for r in orbital_positions(w90tbm)
        R³_cart = Rm * _pad_to_3d_vec(r) # Cartesian coordinates
        Printf.format(io, fmt, R³_cart[1], R³_cart[2], R³_cart[3])
    end
end

function _write_wannier90_win(
    io::IO, w90tbm::Wannier90TightBindingModel, header::String, Rm::AbstractMatrix
)
    N = orbital_count(w90tbm)
    println(io, "! ", header)
    println(io, "num_wann  = ", N)
    println(io, "num_bands = ", N)
    println(io)
    println(io, "begin unit_cell_cart")
    println(io, "ang")
    fmt = Printf.Format("% .10f % .10f % .10f\n")
    for n in 1:3
        Printf.format(io, fmt, Rm[1, n], Rm[2, n], Rm[3, n])
    end
    println(io, "end unit_cell_cart")
end
