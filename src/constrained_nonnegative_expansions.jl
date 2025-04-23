"""
    find_all_admissible_expansions(
        basis::AbstractVector{<:AbstractVector{<:Integer}}
        basis_occupations::AbstractVector{<:Integer},
        occupation::Integer,
        constraints::AbstractVector{<:Integer},
        idxs::AbstractVector{<:Integer};
        basis_idxs = eachindex(basis),
        maxdepth = div(occupation, minimum(basis_occupations), RoundDown)
        )

Given a basis of vectors `basis` ``= [𝐧₁, 𝐧₂, ...]`` with associated non-negative, integer
"occupations" `basis_occupations` ``= [μ₁, μ₂, ...]``, find *all* admissible expansion
coefficients ``{𝐜ⁱ} = {[c₁ⁱ, c₂ⁱ, ...]}`` and associated expansions

``c₁ⁱ𝐧₁ + c₂ⁱ𝐧₂ + ... = 𝐧``

such that ``𝐧`` satisfies the constraints:

1. *occupation constraint*:  each expansion's total occupation ``μ`` is exactly equal to
   `occupation` (i.e., satisfies a linear Diophantine equation):

``c₁ⁱμ₁ + c₂ⁱμ₂ + ... = μ``

2. *symmetry constraint*: each expansion satisfies a set of non-negative, integer
   constraints specified by `constraints`, s.t.:

``(𝐧 = c₁𝐧₁ᴴ + c₂𝐧₂ᴴ + ...)```[idxs][j]` ``≥`` `constraints[j]`

for all `j ∈ eachindex(constraints)`.

# Keyword arguments
- `basis_idxs`: Optionally, if the caller wants to restrict the expansion to a subset of the
  bases in `basis`, the argument `basis_idxs` can provide an indexing into allowable
  elements of `basis`.
- `maxdepth`: include at most `maxdepth` basis vectors, counted with multiplicity (see
  Implementation notes below).

# Implementation
Recursion is used to build a nested set of for loops, of depth `maxdepth`, corresponding 
to the inclusion of at most `maxdepth` basis vectors (this limits the maximum meaningful 
value of `maxdepth` to `div(μ, minimum(μⱼ), RoundDown)`; its default value). 
"""
function find_all_admissible_expansions(
    basis::AbstractVector{<:AbstractVector{<:Integer}},
    basis_occupations::AbstractVector{<:Integer},
    occupation::Integer,
    constraints::AbstractVector{<:Integer},
    idxs::AbstractVector{<:Integer};
    basis_idxs=eachindex(basis),
    maxdepth=div(occupation, minimum(basis_occupations), RoundDown)
)

    occupation > 0 || throw(DomainError(occupation, "must be positive"))

    cⁱs = Vector{Int}[] # solution vector storage
    constraints′ = similar(constraints) # buffer
    _constrained_expansions!(
        cⁱs, constraints′, (), occupation, constraints, basis_occupations, basis, idxs,
        1, length(basis_idxs), 1, maxdepth, basis_idxs)
end
function _constrained_expansions!(
    cⁱs, constraints′, ijks, occupation, constraints, basis_occupations, basis,
    idxs, startidx, stopidx, depth, maxdepth, basis_idxs)
    depth > maxdepth && return cⁱs
    for idxᵢ in startidx:stopidx
        i = basis_idxs[idxᵢ]
        μ = test_expansion_add_if_valid!(cⁱs, constraints′, (ijks..., i), occupation,
            constraints, basis_occupations, basis, idxs)
        μ ≥ occupation && continue # matched/overflowed occupation-constraint; no more to add

        # did not yet match/overflow filling constraint: add more Hilbert basis vectors
        _constrained_expansions!(
            cⁱs, constraints′, (ijks..., i), occupation, constraints, basis_occupations,
            basis, idxs, idxᵢ, stopidx, depth + 1, maxdepth, basis_idxs)
    end
    return cⁱs
end

function test_expansion_add_if_valid!(
    cⁱs, constraints′, # push to cⁱs; use constraints′ as an updating buffer
    ijks::NTuple{N,Int}, occupation, constraints, basis_occupations, basis, idxs
) where {N}

    μ = _sum_fillings(ijks, basis_occupations) # accumulate band fillings
    μ ≠ occupation && return μ                 # return early if μ overflows `occupation`

    # update `constraints′`
    _update_constraints!(constraints′, ijks, constraints, basis, idxs)

    # check if nᵢ+nⱼ+nₖ+... fulfil constraints from `constraints`
    if all(≤(0), constraints′) # check if nᵢ+nⱼ+nₖ+... fulfill `constraints`
        add_solution!(cⁱs, ijks) # push a solution "i+j+k+..." to storage `cⁱs`
    end

    return μ # return occupation associated with `ijks` expansion
end

# equivalent of μ = μs[i] + μs[j] + μs[k] + ... for i,j,k, in ijks, recursively
_sum_fillings(ijks::NTuple{1,Int}, μs) = μs[first(ijks)]
function _sum_fillings(ijks::NTuple{N,Int}, μs) where {N}
    μs[first(ijks)] + _sum_fillings(Base.tail(ijks), μs)
end

# update constraints, assigning to `constraints′`
@inline function _update_constraints!(
    constraints′, ijks::NTuple{N,Int}, constraints, basis, idxs
) where {N}

    b = basis
    c = constraints
    c′ = constraints′
    if N == 1
        i, = ijks
        @views c′ .= c .- b[i][idxs]
    elseif N == 2
        i, j = ijks
        @views c′ .= c .- b[i][idxs] .- b[j][idxs]
    elseif N == 3
        i, j, k = ijks
        @views c′ .= c .- b[i][idxs] .- b[j][idxs] .- b[k][idxs]
    elseif N == 4
        i, j, k, l = ijks
        @views c′ .= c .- b[i][idxs] .- b[j][idxs] .- b[k][idxs] .- b[l][idxs]
    elseif N == 5
        i, j, k, l, o = ijks
        @views c′ .= c .- b[i][idxs] .- b[j][idxs] .- b[k][idxs] .- b[l][idxs] .- b[o][idxs]
    else # fall back to looping
        c′ .= c
        for ijk in ijks
            @views c′ .-= b[ijk][idxs]
        end
    end
    return c′
end

function add_solution!(cⁱs::Vector{Vector{Int}}, ijks::NTuple{N,Int}) where {N}
    # push `ijks` to solution storage `cⁱs` as a vector of indices
    push!(cⁱs, [idx for idx in ijks])
end

# we bother to optimize this, as it can be a bottleneck; much faster than a naive 
# implementation like `sum(basis[idxs])`
function add_basis_vecs!(n, basis::AbstractVector{<:AbstractVector{<:Integer}}, idxs)
    Nⁱʳʳ = length(n)
    copy!(n, basis[first(idxs)]) # peel 1st iter & ensure invariance to n's initialization
    @inbounds for idx in @view idxs[2:end]
        nᴴ = basis[idx]
        for i in 1:Nⁱʳʳ
            n[i] += nᴴ[i]
        end
    end
    return n
end
function add_basis_vecs(basis::AbstractVector{<:AbstractVector{<:Integer}}, idxs)
    n = similar(first(basis))
    return add_basis_vecs!(n, basis, idxs)
end

function coef2idxs(c::AbstractVector{<:Integer})
    N = sum(c)
    cⁱ = Vector{Int}(undef, N)
    pos₁, pos₂, idx = 0, 0, 0
    while true
        idx = findnext(≠(0), c, idx + 1)
        pos₁ = pos₂ + 1
        pos₂ = pos₂ + c[idx]
        cⁱ[pos₁:pos₂] .= idx
        pos₂ == N && break
    end
    return cⁱ
end

function idxs2coef(cⁱ, N_basis) # `N_basis = length(basis)`
    c = zeros(Int, N_basis)
    for i in cⁱ
        c[i] += 1
    end
    return c
end

function isvalid_solution(
    cⁱ::AbstractVector{<:Integer},
    occupation::Integer,
    constraints::AbstractVector{<:Integer},
    basis::AbstractVector{<:AbstractVector{<:Integer}},
    idxs::AbstractVector{<:Integer}
)
    n = add_basis_vecs(basis, cⁱ)
    return all(n[idxs] .≥ constraints) && n[end] == occupation
end