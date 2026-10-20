# ---------------------------------------------------------------------------------------- #
# Constraints on the coefficient vector `t`
#
# All constraints on `t` (spatial symmetry, time-reversal, and hermiticity) are collected as
# rows of a single sparse, real constraint matrix `A` acting on the re/im-doubled coefficient
# vector `[tᴿ; tᴵ]` (with `t = tᴿ + i tᴵ`). The allowed coefficient vectors span the null
# space of `A`, which we obtain from its Gram matrix `AᵀA` (see `constraint_nullspace!`).
# In this formulation, constraints need not be deduplicated (repeated or collinear rows
# leave the null space unchanged), and the null space of several kinds of constraints is
# simply that of their combined rows.

"""
    ConstraintRows(J::Int)

`(i,j,v)` "triplet"-format of the rows of a sparse, real constraint matrix acting on the
re/im-doubled coefficient vector `[tᴿ; tᴵ]` (`2J` elements, with `t = tᴿ + i tᴵ`), intended
to accumulate constraints, block by block, via `add_constraints!` and
`add_complex_constraints!`; entries added repeatedly to the same element are summed.
The accumulated constraint matrix is obtained via `sparse(::ConstraintRows)` (or, consuming
the `ConstraintRows`, via `sparse!(::ConstraintRows)`).
"""
mutable struct ConstraintRows
    const J::Int              # no. of (complex) elements of `t` (half the no. of columns of realified constraint matrix)
    nrows::Int                # total no. of rows/constraints
    const is::Vector{Int}     # rows i
    const js::Vector{Int}     # columns j
    const vs::Vector{Float64} # nonzero values at [i,j]
end
ConstraintRows(J::Int) = ConstraintRows(J, 0, Int[], Int[], Float64[])

function SparseArrays.sparse(constraints::ConstraintRows)
    (; is, js, vs, nrows, J) = constraints
    return sparse(is, js, vs, nrows, 2J) # sums duplicate entries
end

"""
    sparse!(constraints::ConstraintRows) --> SparseMatrixCSC{Float64, Int}

Return the accumulated constraint matrix, as `sparse(constraints)`, but reusing the buffers
of `constraints` as the row-index and value storage of the returned matrix, rather than
allocating new ones.

!!! warning
    `constraints` is consumed: its buffers `is`, `js`, and `vs` are resized and overwritten
    to hold the returned matrix's data, so `constraints` no longer represents the
    constraints and must not be used afterwards.
"""
function SparseArrays.sparse!(constraints::ConstraintRows)
    (; is, js, vs, nrows, J) = constraints
    return sparse!(is, js, vs, nrows, 2J) # sums duplicate entries
end

"""
    RowPermutation(σ::Vector{Int})

A permutation matrix `Π` with `Π[σ[i], i] = 1`, stored (lazily) via the permutation `σ`.

Its product `Π * M` with a sparse matrix `M` is a `RowPermuted`, rather than a new sparse
matrix: the permutation of the rows, `(Π * M)[σ[i], :] = M[i, :]`, is only applied when the
product is added to a `ConstraintRows` (in a `RealConstraintTerm` or
`ComplexConstraintTerm`). `transpose(Π)` is the
inverse permutation, and `sparse(Π)` and `Matrix(Π)` materialize `Π`.
"""
struct RowPermutation
    σ::Vector{Int} # `Π[σ[i], i] = 1`
    function RowPermutation(σ::Vector{Int})
        isperm(σ) || throw(ArgumentError("input is not a permutation"))
        return new(σ)
    end
end
Base.size(Π::RowPermutation) = (length(Π.σ), length(Π.σ))
Base.size(Π::RowPermutation, d::Integer) = d ≤ 2 ? length(Π.σ) : 1
Base.transpose(Π::RowPermutation) = RowPermutation(invperm(Π.σ))
SparseArrays.sparse(Π::RowPermutation) = sparse(Π.σ, eachindex(Π.σ), one(Int), size(Π)...)
Base.Matrix(Π::RowPermutation) = Matrix(sparse(Π))

struct RowPermuted # `Π * M`, for `Π::RowPermutation`
    Π::RowPermutation
    M::SparseMatrixCSC{Int, Int}
end
function Base.:*(Π::RowPermutation, M::SparseMatrixCSC{Int, Int})
    size(Π, 2) == size(M, 1) || throw(DimensionMismatch("incompatible sizes of `Π` & `M`"))
    return RowPermuted(Π, M)
end
Base.size(A::RowPermuted, d...) = size(A.M, d...)

"""
    RealConstraintTerm(a::Real, b::Real, M)

A term `[a M  b M]` of a block of real constraint rows acting on `[tᴿ; tᴵ]` (see
`add_constraints!`), with `M` a V×J matrix: an element of `Mm` (a `SparseMatrixCSC`), or a
row permutation of one (`Π * M`, a `RowPermuted`).
"""
struct RealConstraintTerm{A <: Union{SparseMatrixCSC{Int, Int}, RowPermuted}}
    a::Float64 # coefficient of `M` acting on `tᴿ`
    b::Float64 # coefficient of `M` acting on `tᴵ`
    M::A
end
RealConstraintTerm(a::Real, b::Real, M::A) where {A} = RealConstraintTerm{A}(a, b, M)

"""
    ComplexConstraintTerm(c::Number, M)

A term `c M` of a block of complex constraint rows acting on the complex coefficient vector
`t` (see `add_complex_constraints!`), with `M` a V×J matrix: an element of `Mm` (a
`SparseMatrixCSC`), or a row permutation of one (`Π * M`, a `RowPermuted`).
"""
struct ComplexConstraintTerm{A <: Union{SparseMatrixCSC{Int, Int}, RowPermuted}}
    c::ComplexF64
    M::A
end
ComplexConstraintTerm(c::Number, M::A) where {A} = ComplexConstraintTerm{A}(c, M)

# push the nonzero entries of `x * M` to `constraints`, with their row & column indices
# offset by `row_offset` & `col_offset` (and the rows permuted by `σ`, if not `nothing`);
# `x = 0` would only add explicit zeros, so is skipped
function _push_scaled!(
    constraints::ConstraintRows,
    row_offset::Int,
    col_offset::Int,
    x::Real,
    M::SparseMatrixCSC{Int, Int},
    σ::Union{Nothing, Vector{Int}} = nothing,
)
    iszero(x) && return constraints
    rv, nzv = rowvals(M), nonzeros(M)
    for j in axes(M, 2), p in nzrange(M, j)
        i = isnothing(σ) ? rv[p] : σ[rv[p]]
        push!(constraints.is, row_offset + i)
        push!(constraints.js, col_offset + j)
        push!(constraints.vs, x * nzv[p])
    end
    return constraints
end
_push_scaled!(constraints, row_offset, col_offset, x, A::RowPermuted) =
    _push_scaled!(constraints, row_offset, col_offset, x, A.M, A.Π.σ)

"""
    add_constraints!(constraints::ConstraintRows, terms, [more_terms]) --> ConstraintRows

Append a block of `V` real constraint rows `Σₖ [aₖ Mₖ  bₖ Mₖ]`, acting on `[tᴿ; tᴵ]`, to
`constraints`. `terms` (nonempty) and `more_terms` are vectors or tuples of
`RealConstraintTerm`s, with coefficients `aₖ` & `bₖ` and `V×J` matrices `Mₖ`, all summed
into the same block. Two collections allow terms of different types (e.g., with and
without a row permutation) to be added together while each collection remains concretely
typed.
"""
function add_constraints!(constraints::ConstraintRows, terms, more_terms = ())
    V = size(first(terms).M, 1)
    _push_terms!(constraints, V, terms)
    _push_terms!(constraints, V, more_terms)
    constraints.nrows += V
    return constraints
end
function _push_terms!(constraints::ConstraintRows, V::Int, terms)
    (; J, nrows) = constraints
    for (; a, b, M) in terms::Union{AbstractVector{<:RealConstraintTerm},
                                     Tuple{Vararg{RealConstraintTerm}}}
        size(M) == (V, J) || throw(DimensionMismatch("terms must be V×J matrices"))
        _push_scaled!(constraints, nrows, 0, a, M) # `aₖ Mₖ`, acting on `tᴿ`
        _push_scaled!(constraints, nrows, J, b, M) # `bₖ Mₖ`, acting on `tᴵ`
    end
    return constraints
end

"""
    add_complex_constraints!(constraints::ConstraintRows, terms, [more_terms]) --> ConstraintRows

Append a block of `V` complex constraint rows `C = Σₖ cₖ Mₖ`, acting on the complex
coefficient vector `t`, to `constraints`. `terms` (nonempty) and `more_terms` are vectors or
tuples of `ComplexConstraintTerm`s, with complex coefficients `cₖ` and `V×J` matrices `Mₖ`,
all summed into the same block (cf. `add_constraints!`).

The rows are added in realified form, i.e., as the `2V` real rows of `[Re(C) -Im(C);
Im(C) Re(C)]`, acting on `[tᴿ; tᴵ]`: `C t = 0` is equivalent to `Re(C) tᴿ - Im(C) tᴵ = 0`
and `Im(C) tᴿ + Re(C) tᴵ = 0`.
"""
function add_complex_constraints!(constraints::ConstraintRows, terms, more_terms = ())
    V = size(first(terms).M, 1)
    _push_complex_terms!(constraints, V, terms)
    _push_complex_terms!(constraints, V, more_terms)
    constraints.nrows += 2V
    return constraints
end
function _push_complex_terms!(constraints::ConstraintRows, V::Int, terms)
    (; J, nrows) = constraints
    for (; c, M) in terms::Union{AbstractVector{<:ComplexConstraintTerm},
                                  Tuple{Vararg{ComplexConstraintTerm}}}
        size(M) == (V, J) || throw(DimensionMismatch("terms must be V×J matrices"))
        # NB: `cₖ` is often purely real or imaginary, so some of these are skipped
        _push_scaled!(constraints, nrows,     0,  real(c), M) # Re(C) tᴿ
        _push_scaled!(constraints, nrows,     J, -imag(c), M) # -Im(C) tᴵ
        _push_scaled!(constraints, nrows + V, 0,  imag(c), M) # Im(C) tᴿ
        _push_scaled!(constraints, nrows + V, J,  real(c), M) # Re(C) tᴵ
    end
    return constraints
end

"""
    constraint_nullspace!(constraints::ConstraintRows; atol::Real) --> Matrix{Float64}

Return an orthonormal basis (as columns) for the null space of the constraint matrix
`A = sparse(constraints)`, i.e., for the coefficient vectors `[tᴿ; tᴵ]` that satisfy every
constraint. `constraints` is consumed in the process (see `sparse!(::ConstraintRows)`).

The null space is obtained from that of the Gram matrix `AᵀA`, via a pivoted Cholesky
factorization (see `psd_nullspace!`), treating pivots below `atol²` as zero: since the
eigenvalues of `AᵀA` are the squared singular values of `A`, this roughly corresponds to
treating singular values of `A` below `atol` as zero, as in `nullspace(A; atol)`. `AᵀA`
is only 2J×2J, whereas `A` typically has many more rows than columns.
"""
function constraint_nullspace!(constraints::ConstraintRows; atol::Real)
    A = sparse!(constraints)
    G = Symmetric(Matrix(A' * A))
    # `A` & `G = AᵀA` share their null space (`Gx = 0` ⇒ `‖Ax‖² = xᵀGx = 0`); `G` is only
    # 2J×2J and positive semi-definite, so its null space is cheaply found by a pivoted
    # Cholesky factorization (cheaper, and more robust for the degenerate null eigenvalues,
    # than partial eigensolvers; and much cheaper than an SVD of the tall `A`). Working with
    # `G` squares `A`'s condition number, which is harmless here: `A`'s nonzero singular
    # values lie in ~[1, 5] (empirically, across all space groups), far from `atol`.
    return psd_nullspace!(G; tol = atol^2)
end

"""
    obtain_basis_free_parameters(
        h_orbit::HoppingOrbit{D},
        brₐ::BandRep{D}, 
        brᵦ::BandRep{D}, 
       [orderingₐ = OrbitalOrdering(brₐ), 
        orderingᵦ = OrbitalOrdering(brᵦ),
        diagonal_block::Bool = true,
        hermiticity::Hermiticity = HERMITIAN]
    )              --> Tuple{Matrix{SparseMatrixCSC{Int, Int}}, Vector{Vector{Float64}}}

Obtain the basis of free parameters for the hopping terms between `brₐ` and `brᵦ` 
associated with the hopping orbit `h_orbit`.

## Note

The presence or absence of time-reversal symmetry is inferred implicitly from `brₐ` and `brᵦ`.
"""
function obtain_basis_free_parameters(
    h_orbit::HoppingOrbit{D},
    brₐ::BandRep{D},
    brᵦ::BandRep{D},
    orderingₐ::OrbitalOrdering{D} = OrbitalOrdering(brₐ),
    orderingᵦ::OrbitalOrdering{D} = OrbitalOrdering(brᵦ),
    diagonal_block::Bool = true,
    hermiticity::Hermiticity = HERMITIAN,
) where {D}
    num(brₐ) == num(brᵦ) || error("input `BandRep`s must come from the same space group")
    brₐ.timereversal == brᵦ.timereversal || error("input `BandRep`s must have identical time-reversal symmetry")
    timereversal = brₐ.timereversal

    # generators in a primitive basis (double group generators, for spinful band
    # representations)
    gens = primitivized_generators(brₐ)

    # encode Hamiltonian as a coefficient matrix sandwiched by exponentials & hopping ampl.
    Mm = construct_M_matrix(h_orbit, brₐ, brᵦ, orderingₐ, orderingᵦ)

    # actual computation
    tₐᵦ_basis_reim = _obtain_basis_free_parameters(
        h_orbit,
        brₐ,
        brᵦ,
        orderingₐ,
        orderingᵦ,
        Mm,
        gens,
        timereversal,
        diagonal_block,
        hermiticity,
    )

    return Mm, tₐᵦ_basis_reim
end

function _obtain_basis_free_parameters(
    h_orbit::HoppingOrbit{D},
    brₐ::BandRep{D},
    brᵦ::BandRep{D},
    orderingₐ::OrbitalOrdering{D},
    orderingᵦ::OrbitalOrdering{D},
    Mm::AbstractMatrix{SparseMatrixCSC{Int, Int}},
    gens::AbstractVector{<:AbstractOperation{D}},
    timereversal::Bool,
    diagonal_block::Bool,
    hermiticity::Hermiticity,
) where {D}
    # collect all constraints on the coefficient vector `t` associated with `h_orbit`, as
    # rows of a real constraint matrix acting on the re/im-doubled vector `[tᴿ; tᴵ]`: we
    # want to work with strictly real coefficients for the basis vectors, so we split the
    # possible products of complex coefficients with a complex basis vector into two
    # vectors, where the lower half of each represents the imaginary part of the product.
    # A benefit of this is that it becomes possible to take the complex conjugate of each
    # element "symbolically" by simply multiplying with [1 0; 0 -1] (block matrix), which
    # is needed for the time-reversal and hermiticity constraints. Crucially, this must be
    # kept in mind everywhere that we use a coefficient vector. In particular, we must
    # always interpret `Mm` as acting like `[Mm iMm]` on a "doubled" basis vector.
    constraints = ConstraintRows(coefficient_count(Mm))

    # spatial symmetry constraints on Hₐᵦ
    add_symmetry_constraints!(constraints, Mm, h_orbit, brₐ, brᵦ, gens, orderingₐ, orderingᵦ)

    # time-reversal constraints
    if timereversal
        add_timereversal_constraints!(constraints, Mm, brₐ, brᵦ, orderingₐ, orderingᵦ)
    end

    # hermiticity constraints, if this is a diagonal block
    if diagonal_block && hermiticity ≠ NONHERMITIAN
        antihermitian = hermiticity === ANTIHERMITIAN
        add_hermiticity_constraints!(constraints, Mm, h_orbit, antihermitian)
    end

    # the allowed coefficient vectors span the null space of the combined constraints
    tₐᵦ_basis_reim_matrix = constraint_nullspace!(constraints; atol = NULLSPACE_ATOL_DEFAULT)
    isempty(tₐᵦ_basis_reim_matrix) && return Vector{Vector{Float64}}()

    # ------------------------------------------------------------------------------------ #
    # Make sure we have a reasonably pretty-looking basis in the end by sparsifying and
    # dropping near-zero elements explicitly

    # convert null-space to a sparse column form
    tₐᵦ_basis_reim_matrix_form′ = _poormans_sparsification(tₐᵦ_basis_reim_matrix)
    tₐᵦ_basis_reim = [collect(v) for v in eachcol(tₐᵦ_basis_reim_matrix_form′)]

    # prune near-zero elements of basis vectors
    _prune_at_threshold!(tₐᵦ_basis_reim)

    return tₐᵦ_basis_reim
end

"""
Poor man's "matrix sparsification" via the reduced row echelon form.

!!! warning
    This function is an internal helper function for `obtain_basis_free_parameters` and
    is not part of the public API.
"""
function _poormans_sparsification(
    A::AbstractMatrix{<:Number};
    rref_tol::Union{Nothing, Float64} = SPARSIFICATION_ATOL_DEFAULT,
)
    # following appendix E of the Qsymm paper (https://arxiv.org/abs/1806.08363) [copied
    # over from Neumann.jl] [copied over from MPBUtils.jl]
    if !isnothing(rref_tol)
        # use a relatively low tolerance in `rref` to avoid explosions of errors
        # NB: this optional tolerance argument of `rref!` is undocumented :(
        return transpose(rref!(copy(transpose(A)), rref_tol))
    end
    return transpose(rref(transpose(A)))
end

"""
Prune near-zero elements of vectors in `vs`.

!!! warning
    This function is an internal helper function for `obtain_basis_free_parameters` and
    is not part of the public API.
"""
function _prune_at_threshold!(
    vs::AbstractVector{<:AbstractVector{T}};
    atol::Real = PRUNE_ATOL_DEFAULT,
) where {T <: Complex}
    for v in vs
        for (j, vⱼ) in enumerate(v)
            rⱼ, iⱼ = reim(vⱼ)
            rⱼ′ = ifelse(abs(rⱼ) < atol, zero(real(T)), rⱼ)
            iⱼ′ = ifelse(abs(iⱼ) < atol, zero(real(T)), iⱼ)
            v[j] = T(rⱼ′, iⱼ′)
        end
    end
    return vs
end

function _prune_at_threshold!(
    vs::AbstractVector{<:AbstractVector{T}};
    atol::Real = PRUNE_ATOL_DEFAULT,
) where {T <: Real} # needed for type consistency in TRS implementation
    for v in vs
        for (j, vⱼ) in enumerate(v)
            v[j] = abs(vⱼ) < atol ? zero(T) : vⱼ
        end
    end
    return vs
end
