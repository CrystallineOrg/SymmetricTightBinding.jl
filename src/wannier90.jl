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

# ---------------------------------------------------------------------------------------- #
# Reading `_hr.dat` files

"""
    read_wannier90(prefix::AbstractString, positions::AbstractVector{DirectPoint{D}})
                                                        --> Wannier90TightBindingModel{D}

Read the real-space hopping matrices ``H(\\mathbf{R})`` of a Wannier90
[`<prefix>_hr.dat`](https://wannier90.readthedocs.io/en/latest/user_guide/wannier90/files/#seedname_hrdat)
file, e.g., as obtained from a Wannierization of a DFT calculation.

The returned hopping matrices reproduce Wannier90's own interpolated Bloch Hamiltonian,
i.e., ``\\sum_{\\mathbf{R}} H(\\mathbf{R}) \\mathrm{e}^{2π\\mathrm{i}\\mathbf{k} \\cdot
\\mathbf{R}}`` equals the Hamiltonian that Wannier90 interpolates.
To this end, the hopping matrices are divided by the Wigner–Seitz degeneracy weights listed
in `_hr.dat` and, if present, `<prefix>_wsvec.dat` is read and each matrix element
``H_{ij}(\\mathbf{R})`` is moved to the minimal-image translations
``\\mathbf{R} + \\mathbf{T}`` it lists (relevant if Wannier90 ran with
`use_ws_distance = true`, its default).
Minimal-image translations are generally needed to obtain hopping matrices that respect the
crystal symmetry.

`positions` are the orbital positions (Wannier centres), in the (primitive) lattice basis
that the translations ``\\mathbf{R}`` refer to, and in the order of the Wannier functions.
They are not featured in `_hr.dat` but are needed to evaluate the model in
SymmetricTightBinding's Fourier convention (see [`Wannier90TightBindingModel`](@ref)).
Ideally, they are the Wyckoff positions of the band representation whose orbitals the
Wannier functions represent (e.g., `orbital_positions(cbr)`).

For `D < 3`, translations with nonzero components along the missing dimensions are dropped.
"""
function read_wannier90(
    prefix::AbstractString, positions::AbstractVector{DirectPoint{D}}
) where D
    D ≤ 3 || error("only D ≤ 3 supported")
    N = length(positions) # number of orbitals, i.e., of Wannier functions

    # Wannier90 interpolates the Bloch Hamiltonian from the contents of `_hr.dat` and
    # `_wsvec.dat` as (in its lattice-phase-only Fourier convention):
    #     Hᵢⱼ(k) = ∑_R ∑_{T ∈ 𝒯ᵢⱼ(R)} Hᵢⱼ(R) e^{2πik·(R+T)} / (deg(R) |𝒯ᵢⱼ(R)|),
    # with ∑_R runs over the translations `R` listed in `_hr.dat`, and where
    #   - `deg(R)` is the "degeneracy weight" of `R`, listed in `_hr.dat` (see
    #     `_parse_wannier90_hr` below for its meaning),
    #   - `𝒯ᵢⱼ(R)` are the minimal-image shifts of `R` for the orbital pair `(i, j)`, listed
    #     in `_wsvec.dat` (see `_parse_wannier90_wsvec` below; if this file is absent, or if
    #     Wannier90 ran with `use_ws_distance = false`, 𝒯ᵢⱼ(R) = {0}).
    # So, as a real-space model, each listed element Hᵢⱼ(R) amounts to a hopping amplitude
    # Hᵢⱼ(R) / (deg(R) |𝒯ᵢⱼ(R)|) at each of the translations `R + T`, T ∈ 𝒯ᵢⱼ(R), and we
    # reconstruct the model's hopping matrices by accumulating these amplitudes

    # parse the files: `hr` holds the contents of `_hr.dat` as is, i.e., `N`, and, per
    # listed translation `R`, its degeneracy weight `deg(R)` and hopping matrix `H(R)`;
    # `shifts[(R, i, j)]` holds 𝒯ᵢⱼ(R)
    hr = _parse_wannier90_hr(prefix * "_hr.dat")
    hr.N == N || throw(DimensionMismatch(lazy"`_hr.dat` has $(hr.N) Wannier functions, but $N positions were given"))
    path_wsvec = prefix * "_wsvec.dat"
    shifts = isfile(path_wsvec) ? _parse_wannier90_wsvec(path_wsvec) : nothing

    # reconstruct the model's hopping matrices
    hoppings = Dict{SVector{D, Int}, Matrix{ComplexF64}}() # translation => hopping matrix
    for (R, deg, H) in zip(hr.translations, hr.degeneracies, hr.hopping_matrices)
        for j in 1:N, i in 1:N
            Ts = if isnothing(shifts)
                [zero(SVector{3, Int})] # no `_wsvec.dat`: the hopping stays at `R`
            else
                get(shifts, (R, i, j)) do
                    error(lazy"`_wsvec.dat` lacks an entry for (R, i, j) = ($R, $i, $j)")
                end
            end
            v = H[i, j] / (deg * length(Ts))
            for T in Ts
                R′ = R + T
                # for `D < 3`, drop translations along the missing dimensions
                all(iszero, R′[D+1:3]) || continue
                R′ᴰ = SVector{D, Int}(ntuple(d -> R′[d], Val(D)))
                get!(() -> zeros(ComplexF64, N, N), hoppings, R′ᴰ)[i, j] += v
            end
        end
    end

    translations = sort!(collect(keys(hoppings)))
    return Wannier90TightBindingModel{D}(translations, [hoppings[R] for R in translations],
                                         collect(positions))
end

# Parse a Wannier90 `_hr.dat` file, returning its contents as is: the number of Wannier
# functions `N` and, for each listed lattice translation `R` (in the lattice basis, as a 3D
# integer vector), its degeneracy weight `deg(R)` and its `N×N` hopping matrix `H(R)`.
#
# File layout (see Wannier90's user guide):
#   line 1:                header (e.g., the creation date; unused)
#   line 2:                `N`, the number of Wannier functions
#   line 3:                `Nᴿ`, the number of listed translations
#   next ⌈Nᴿ/15⌉ lines:    the `Nᴿ` degeneracy weights, 15 per line
#   next `Nᴿ⋅N²` lines:    records `R₁ R₂ R₃ i j Re(Hᵢⱼ(R)) Im(Hᵢⱼ(R))`, grouped by `R`:
#                          `N²` consecutive records per translation, one per `(i, j)`
#
# The degeneracy weights: Wannier90 obtains `H(R)` by Fourier transforming the Hamiltonian
# from a k-mesh of `N₁×N₂×N₃` points, so `H(R)` is only determined up to translations of `R`
# by the mesh's supercell (spanned by `N₁a₁`, `N₂a₂`, `N₃a₃`). Wannier90 picks the shortest
# representatives, i.e., the translations `R` in the Wigner–Seitz cell of this supercell.
# A translation on the cell's boundary is equally short as `deg(R) - 1` translations that
# are equivalent to it (i.e., differ from it by supercell translations); these are listed
# too, and the hopping is split evenly over all `deg(R)` of them: hence, the division by
# `deg(R)` in `read_wannier90`.
function _parse_wannier90_hr(path::AbstractString)
    open(path) do io
        _readline_or_eof_error(io, path) # line 1: header
        N = parse(Int, _readline_or_eof_error(io, path)) # line 2
        Nᴿ = parse(Int, _readline_or_eof_error(io, path)) # line 3

        # degeneracy weights, 15 per line
        degeneracies = Int[]
        while length(degeneracies) < Nᴿ
            append!(degeneracies, parse.(Int, split(_readline_or_eof_error(io, path))))
        end
        length(degeneracies) == Nᴿ || error(lazy"`$path` lists $(length(degeneracies)) degeneracy weights; expected $Nᴿ")

        # hopping matrices: a block of `N²` records per translation; we check that each
        # block is complete (every `(i, j)` once) and refers to a single translation
        translations = Vector{SVector{3, Int}}(undef, Nᴿ)
        hopping_matrices = [zeros(ComplexF64, N, N) for _ in 1:Nᴿ]
        for n in 1:Nᴿ
            seen = falses(N, N) # elements `(i, j)` of the `n`th block read so far
            for record in 1:N^2
                fields = split(_readline_or_eof_error(io, path))
                length(fields) == 7 || error(lazy"malformed record in `$path`: expected 7 fields, got \"$(join(fields, ' '))\"")
                R = SVector{3, Int}(ntuple(d -> parse(Int, fields[d]), Val(3)))
                i, j = parse(Int, fields[4]), parse(Int, fields[5])
                if record == 1
                    translations[n] = R
                elseif R ≠ translations[n]
                    error(lazy"unexpected record in `$path`: translation $R in the block of translation $(translations[n]) (expected $(N^2) records per translation)")
                end
                (1 ≤ i ≤ N && 1 ≤ j ≤ N) || error(lazy"orbital indices ($i, $j) out of range 1:$N in `$path`")
                seen[i, j] && error(lazy"duplicate record for (R, i, j) = ($R, $i, $j) in `$path`")
                seen[i, j] = true
                Hᵢⱼʳᵉ, Hᵢⱼⁱᵐ = parse(Float64, fields[6]), parse(Float64, fields[7])
                hopping_matrices[n][i, j] = complex(Hᵢⱼʳᵉ, Hᵢⱼⁱᵐ)
            end
        end
        return (; N, translations, degeneracies, hopping_matrices)
    end
end

# Parse a Wannier90 `_wsvec.dat` file, returning a dictionary that maps each `(R, i, j)`, a
# translation `R` and an orbital pair `(i, j)` of `_hr.dat`, to its minimal-image shifts
# `𝒯ᵢⱼ(R)`: the hopping `Hᵢⱼ(R)` is to be placed at the translations `R + T`, T ∈ 𝒯ᵢⱼ(R),
# rather than at `R`.
#
# Purpose: as for the degeneracy weights (see `_parse_wannier90_hr`), Wannier90 must pick
# among translations `R` that are equivalent modulo the k-mesh's supercell. In `_hr.dat`,
# it picks the shortest `R`; but the actual hopping, from orbital `j` in cell `R` to orbital
# `i` in the home cell, spans the distance |rⱼ + R - rᵢ|, which also depends on the Wannier
# centres `rᵢ` and `rⱼ`. With `use_ws_distance = true`, Wannier90 instead picks, for each
# `(i, j)`, the shifts `T` (supercell translations) that minimize |rⱼ + R + T - rᵢ|,
# splitting the hopping evenly if several do. This matters for symmetry: hoppings related by
# symmetry have equal lengths, so the minimal-image choice treats them alike; the choice
# based on `R` alone need not, unless `rᵢ = rⱼ`. If `use_ws_distance = false`, all shifts
# are zero.
#
# File layout:
#   line 1:                header (e.g., the creation date and `use_ws_distance`; unused)
#   then, per `(R, i, j)`: a line `R₁ R₂ R₃ i j`; a line with the number of shifts `M`; and
#                          `M` lines `T₁ T₂ T₃`
function _parse_wannier90_wsvec(path::AbstractString)
    shifts = Dict{Tuple{SVector{3, Int}, Int, Int}, Vector{SVector{3, Int}}}()
    open(path) do io
        _readline_or_eof_error(io, path) # line 1: header
        while !eof(io)
            line = readline(io)
            isempty(strip(line)) && continue # e.g., a trailing empty line
            fields = parse.(Int, split(line))
            length(fields) == 5 || error(lazy"malformed line in `$path`: expected `R₁ R₂ R₃ i j`, got \"$line\"")
            R = SVector{3, Int}(fields[1], fields[2], fields[3])
            i, j = fields[4], fields[5]
            M = parse(Int, _readline_or_eof_error(io, path)) # number of shifts
            shifts[(R, i, j)] = map(1:M) do _
                SVector{3, Int}(parse.(Int, split(_readline_or_eof_error(io, path))))
            end
        end
    end
    return shifts
end

# the next line of `io`, read from file `path`; erroring if the end of the file is reached
function _readline_or_eof_error(io::IO, path::AbstractString)
    eof(io) && error(lazy"unexpected end of file `$path`")
    return readline(io)
end

# ---------------------------------------------------------------------------------------- #
# Fitting a symmetric tight-binding model to a `Wannier90TightBindingModel`

"""
    fit_wannier90(tbm::TightBindingModel{D}, w90tbm::Wannier90TightBindingModel{D})
                                                    --> ParameterizedTightBindingModel{D}

Determine the coefficients `cs` of the tight-binding model `tbm` that best reproduce the
real-space hopping matrices of `w90tbm` (e.g., read by [`read_wannier90`](@ref) from a
Wannierization of a DFT calculation), returning `tbm(cs)`.

The fit is a linear least-squares problem in real space, which decouples into a small,
independent problem for each hopping orbit of `tbm`: each matrix element
``H_{ij}(\\mathbf{R})`` is a hopping along ``\\boldsymbol{δ} = \\mathbf{r}_i -
\\mathbf{r}_j - \\mathbf{R}``, so it belongs to a single hopping orbit. For an orbit
``\\{\\boldsymbol{δ}_n\\}`` of the block between two band representations, spanned by
terms with amplitudes `MmtC_t[n,i,j]` (for hopping along ``\\boldsymbol{δ}_n`` from
orbital `j` to orbital `i`), the coefficients ``c_t`` of the orbit's terms minimize
``\\sum_{nij} |H^{\\text{w90}}_{ij}(\\mathbf{r}_i - \\mathbf{r}_j -
\\boldsymbol{δ}_n) - \\sum_t c_t \\mathrm{MmtC}_t[n,i,j]|^2``. I.e., the Wannier90
hopping amplitudes are orthogonally projected onto the symmetry-allowed amplitudes of each
orbit. Amplitudes of `w90tbm` outside the hopping orbits of `tbm` are ignored.

The coefficients are real, as for any `tbm`: its terms span the allowed Hamiltonians as a
real vector space (hermiticity and time-reversal involve complex conjugation, and so are
not closed under complex scaling). Accordingly, the least-squares problems are solved over
the reals.

The orbitals of `w90tbm` and `tbm` are assumed to correspond one-to-one: i.e., the Wannier
functions must have the same order, centres, and local basis (i.e., transform identically
under the site-symmetry group) as the orbitals of `tbm`.
"""
function fit_wannier90(
    tbm::TightBindingModel{D}, w90tbm::Wannier90TightBindingModel{D}
) where D
    positions = orbital_positions(tbm)
    N = length(positions)
    N == orbital_count(w90tbm) || throw(DimensionMismatch(lazy"`tbm` has $N orbitals, but `w90tbm` has $(orbital_count(w90tbm))"))
    for (r, r′) in zip(positions, orbital_positions(w90tbm))
        isapprox(r, r′; atol = VEC_CMP_ATOL) || error(lazy"orbital positions of `tbm` and `w90tbm` differ: $r vs. $r′")
    end
    w90_hoppings = Dict(zip(w90tbm.translations, w90tbm.hopping_matrices)) # R => H(R)

    # the least-squares problem decouples into one problem per hopping orbit of `tbm`; each
    # orbit is spanned by one or more terms of `tbm`, which `tb_hamiltonian` lists
    # consecutively: so, each orbit corresponds to a range of term indices `terms` (e.g.,
    # `1:1`, `2:3`, `4:4`, …; see `_orbit_term_ranges` below), and its fit determines the
    # coefficients `cs[terms]`
    cs = zeros(length(tbm))
    for terms in _orbit_term_ranges(tbm)
        cs[terms] = _fit_hopping_orbit(tbm, terms, w90_hoppings).c
    end
    return tbm(cs)
end

# Least-squares fit of a single hopping orbit of `tbm`, spanned by the terms
# `tbm.terms[terms]`, to the Wannier90 hopping matrices `w90_hoppings` (a dictionary
# `R => H(R)`). Returns the fitted coefficients `c` of the terms, as well as the (Frobenius)
# norms of the targeted Wannier90 amplitudes, `w90_norm`, and of the fit's residual,
# `residual`.
#
# All terms of a hopping orbit {δₙ} belong to the same block `(block_i, block_j)` of `tbm`,
# i.e., to the same pair of band representations, whose orbitals have the (global) indices
# `is` and `js`. Term `t` has the amplitude `MmtC_t[n, local_i, local_j]` (from
# `tbm.terms[t].block.MmtC`) for the hopping along δₙ from orbital `j = js[local_j]` in cell
# `R` to orbital `i = is[local_i]` in the home cell, with `R = rᵢ - rⱼ - δₙ`: i.e., term `t`
# contributes `MmtC_t[n, local_i, local_j]` to the matrix element Hᵢⱼ(R) (as in
# `_add_hopping_matrices!`).
#
# The least-squares problem is `A c ≈ b`, where each row is a matrix element Hᵢⱼ(R) that the
# orbit contributes to, and each column is a term of the orbit:
#   - `A[row, t]` is the amplitude of the orbit's `t`th term at the row's matrix element,
#   - `b[row]` is the Wannier90 amplitude at the row's matrix element,
# to be solved for real coefficients `c`.
function _fit_hopping_orbit(
    tbm::TightBindingModel{D, S}, terms::UnitRange{Int}, w90_hoppings::Dict
) where {D, S}
    positions = orbital_positions(tbm)
    tbt = tbm.terms[first(terms)] # the block & hopping orbit are shared by all `terms`
    block_i, block_j = tbt.block_ij
    is = tbt.axis[Block(block_i)] # global indices of the block's rows
    js = tbt.axis[Block(block_j)] # global indices of the block's columns
    δs = constant.(orbit(tbt.block.h_orbit)) # the orbit's hopping vectors δₙ
    MmtCs = [tbm.terms[t].block.MmtC for t in terms]

    # for (anti-)Hermitian models, each term also contributes to the hermiticity-related
    # block `(block_j, block_i)`, which `tbm` does not list separately: an amplitude `v` at
    # Hᵢⱼ(R) implies an amplitude `conj(v)` (Hermitian) or `-conj(v)` (anti-Hermitian) at
    # Hⱼᵢ(-R); we include those matrix elements as rows too
    # NB: for diagonal blocks (`block_i == block_j`), these are already among the orbit's
    #     own matrix elements, since the orbit then contains both δ and -δ
    has_hermiticity_partner = S !== NONHERMITIAN && block_i ≠ block_j
    partner_sign = S === ANTIHERMITIAN ? -1 : 1

    # assemble `A` and `b`, a row at a time: a row for each `(n, local_i, local_j)` at which
    # some term of the orbit has a nonzero amplitude, i.e., for each matrix element that the
    # orbit contributes to (and a row for its hermiticity partner, if any)
    A_rows = Vector{ComplexF64}[]
    b = ComplexF64[]
    for I in CartesianIndices(first(MmtCs))
        amplitudes = [MmtC[I] for MmtC in MmtCs] # the terms' amplitudes at this element
        all(v -> abs(v) < PRUNE_ATOL_DEFAULT, amplitudes) && continue # not in the orbit
        n, local_i, local_j = Tuple(I)
        i, j = is[local_i], js[local_j]
        R = _integer_translation(positions[i] - positions[j] - δs[n])
        push!(A_rows, amplitudes)
        push!(b, _w90_amplitude(w90_hoppings, R, i, j))
        if has_hermiticity_partner
            push!(A_rows, partner_sign * conj(amplitudes))
            push!(b, _w90_amplitude(w90_hoppings, -R, j, i))
        end
    end
    A = stack(A_rows; dims = 1) # rows: matrix elements; columns: terms

    # solve `A c ≈ b` for real `c`: i.e., w/ real & imag parts of `A c ≈ b` as separate rows
    c = [real(A); imag(A)] \ [real(b); imag(b)]

    return (; c, w90_norm = norm(b), residual = norm(A * c - b))
end

# the Wannier90 amplitude Hᵢⱼ(R); zero if `R` is not among the Wannier90 translations
function _w90_amplitude(w90_hoppings::Dict, R, i::Int, j::Int)
    H = get(w90_hoppings, R, nothing)
    return isnothing(H) ? zero(ComplexF64) : H[i, j]
end

# The terms of `tbm`, grouped by hopping orbit, as ranges of term indices (e.g., `[1:1, 2:3,
# 4:4, …]`): `tb_hamiltonian` lists the terms that span a hopping orbit consecutively (one
# term per basis vector of the orbit's free parameters), so each group is a run of
# consecutive terms with the same block and hopping orbit
function _orbit_term_ranges(tbm::TightBindingModel)
    ranges = UnitRange{Int}[]
    for t in eachindex(tbm.terms)
        if t > 1 && _same_hopping_orbit(tbm.terms[t], tbm.terms[t-1])
            ranges[end] = first(ranges[end]):t # extend the current run
        else
            push!(ranges, t:t) # start a new run
        end
    end
    return ranges
end
function _same_hopping_orbit(tbt::TightBindingTerm, tbt′::TightBindingTerm)
    return tbt.block_ij == tbt′.block_ij && tbt.block.h_orbit === tbt′.block.h_orbit
end
