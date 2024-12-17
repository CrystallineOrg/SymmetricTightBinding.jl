using LinearAlgebra: nullspace
using Crystalline: isapproxin
"""
Compute the symmetry related hoppings relative vectors from the WP of `br1` to the WP of `br2`
displaced a set of primitive lattice vectors `Rs`.

How it works: 
1. Take a point `a` in the WP of `br1` and a point `b` in the WP of `br2`. We compute the 
    displacement vector `δ = b + R - a`, where `R ∈ Rs`.
2. If `δ ∈ representatives` then we add `δ => (a, b, R)` to the list of hoppings of that 
    representative and continue. If not then, we search inside of all the representatives
    for the one that `δ => (a, b, R)` in the list of hoppings. If not found, then we add `δ`
    as a new representative and add `δ =>(a, b, R)` to its list of hoppings.
3. Take `g ∈ generators` and compute `δ' = g δ` and `(a', b', R') = (g a,g b, g R)`, and 
    repeat step 2.
4. 
"""
function obtain_symmetry_related_hoppings(
    Rs::AbstractVector{V}, # must be specified in the primitive basis
    br1::NewBandRep{D},
    br2::NewBandRep{D}
) where {V<:Union{AbstractVector{<:Integer},RVec{D}}} where {D}

    sgnum = num(br1)
    num(br2) == sgnum || error("both band representations must be in the same space group")
    # we only want to include the wyckoff positions in the primitive cell - but the default
    # listings from `spacegroup` include operations that are "centering translations";
    # fortunately, the orbit returned for a `NewBandRep` do not include these redundant
    # operations - but is still specified in a conventional basis. So, below, we remove
    # redundant operations from the space group, and also change both the operations and the
    # positions from a conventional to a primitive basis
    cntr = centering(sgnum, D)
    ops = primitivize(spacegroup(sgnum, D))
    wps1 = primitivize.(orbit(group(br1)), cntr)
    wps2 = primitivize.(orbit(group(br2)), cntr)

    # we are going to create a dictionary of dictionaries. The first set of keys will indicate
    # the representative of the hopping distance and inside it, each dictionary will be 
    # associated with a point in the "orbit" of that hopping distance. Addtionally, if you 
    # look at the value of and specific point inside the "orbit" of one representative you 
    # will obtain a set of tuples that will encode the points of the WPs involved in the 
    # hopping and the displacement vector R. (q_i, w_j , R) => q_i -> w_j + R

    δsdd = Dict{RVec{D},Vector{Pair{RVec{D},Vector{NTuple{3,RVec{D}}}}}}()
    for R in Rs
        R = RVec{D}(R) # change the type of R to be type consistent
        for (qₐ, qᵦ) in Iterators.product(wps1, wps2)
            qₐ = parent(qₐ) # work with RVec directly rather than WyckoffPosition
            qᵦ = parent(qᵦ)
            δ = qᵦ + R - qₐ
            if (!isapproxin(δ, keys(δsdd)) &&
                !isapproxin(δ, Iterators.map(first, Iterators.flatten(values(δsdd)))))
                push!(δsdd, δ => [δ => [(qₐ, qᵦ, R)]])
                for g in ops
                    Ρ = SymOperation(rotation(g)) # type consistency for the rotation 
                    # we want to keep the WPs inside the unit cell and put all the posible
                    # translations into the RVec 'R', so we can keep the notation such that
                    # q_i -> w_j + R. For that reason we make the following changes:
                    _qₐ′ = Ρ * qₐ
                    _qᵦ′ = Ρ * qᵦ
                    qₐ′ = RVec(reduce_translation_to_unitrange(constant(_qₐ′)), free(_qₐ′))
                    qᵦ′ = RVec(reduce_translation_to_unitrange(constant(_qᵦ′)), free(_qᵦ′))
                    dₐ = (Ρ * qₐ) - qₐ′
                    dᵦ = (Ρ * qᵦ) - qᵦ′
                    R′ = (Ρ * R) + dᵦ - dₐ
                    δ′ = Ρ * δ
                    idx = findfirst(v -> first(v) ≈ δ′, δsdd[δ])
                    if !isapproxin(δ′, keys(δsdd)) && isnothing(idx)
                        push!(δsdd[δ], δ′ => [(qₐ′, qᵦ′, R′)])
                    end

                    isnothing(idx) && continue
                    idx = something(idx)
                    # evaluate `bool = (qₐ′, qᵦ′, R′) ∈ last(δsdd[δ][idx])`, w/ approximate
                    # equality comparison
                    bool = any(last(δsdd[δ][idx])) do (qₐ′′, qᵦ′′, R′′)
                        isapprox(qₐ′, qₐ′′) && isapprox(qᵦ′, qᵦ′′) && isapprox(R′, R′′)
                    end
                    if !bool
                        push!(last(δsdd[δ][idx]), (qₐ′, qᵦ′, R′))
                    end
                end
            end
        end
    end
    return δsdd
end

# EBRs: (q|A), (w|B)
# Wyckoff positions: q, w
#   q: q1, ..., qN
#   w: w1, ..., wM
# Site symmetry irreps: A, B
#   A: A1, ..., AJ
#   B: B1, ..., BK
# δs = [δ1, δ2, ..., δn]
#   δ1: qi₁¹ -> wj₁¹, qi₁² -> wj₁², ...
#   δ2: qi₂¹ -> wj₂¹, qi₂² -> wj₂², ...
# v = [exp(ik⋅δ1), exp(ik⋅δ2), ..., exp(ik⋅δn)]
# t = [[t(δ1) ...], [t(δ2) ...], ..., [t(δn) ...]]
#   t(δ1): [t(qi₁ᵅ -> wj₁ᵅ, A_f -> B_g) ...]

# Current example: (1a|E), (2c|A)
#   ___w2__
#  |   x   |
#  |q1 x   x w1
#  |_______|
#   δs = [1/2x, -1/2x, 1/2y, -1/2y]
#      δ1: q1 -> w1 + G1
#      δ2: q1 -> w1 + G2
#      δ3: q1 -> w2 + G3
#      δ4: q1 -> w2 + G4
# t = [t(δ1)..., t(δ2)..., t(δ3)..., t(δ4)...]
#   t(δ1): [t(q1 -> w1, G1, E1 -> A1), t(q1 -> w1, G1, E2 -> A1)]
#   t(δ2): [t(q1 -> w1, G2, E1 -> A1), t(q1 -> w1, G2, E2 -> A1)]
#   t(δ3): [t(q1 -> w2, G3, E1 -> A1), t(q1 -> w2, G3, E2 -> A1)]
#   t(δ4): [t(q1 -> w2, G4, E1 -> A1), t(q1 -> w2, G4, E2 -> A1)]

"""
Gives an order for the Hamiltonian's term concerning hopping between `br1` and `br2`.

A matrix with the dimensions of the Hamiltonian is given back where each term indicates the 
hopping term considered in the Hamiltonian at that position.
"""
function hamiltonian_term_order(
    br1::NewBandRep{D},
    br2::NewBandRep{D},
) where {D}

    sgnum = num(br1)
    num(br2) == sgnum || error("Both band representations must be in the same space group")
    # we only want to include the wyckoff positions in the primitive cell - but the default
    # listings from `spacegroup` include operations that are "centering translations";
    # fortunately, the orbit returned for a `NewBandRep` do not include these redundant
    # operations - but is still specified in a conventional basis. So, below, we remove
    # redundant operations from the space group, and also change both the operations and the
    # positions from a conventional to a primitive basis
    cntr = centering(sgnum, D)
    wp1 = primitivize.(orbit(group(br1)), cntr)
    wp2 = primitivize.(orbit(group(br2)), cntr)

    V1, V2, = length(wp1), length(wp2)
    Q1, Q2 = irdim(br1.siteir), irdim(br2.siteir)
    order = Matrix{Pair{Tuple{Int64,WyckoffPosition{D}},
        Tuple{Int64,WyckoffPosition{D}}}}(undef, V1 * Q1, V2 * Q2)
    for i in 1:V1
        for j in 1:V2
            for k in 1:Q1
                for l in 1:Q2
                    order[i*k, j*l] = (k, wp1[i]) => (l, wp2[j])
                end
            end
        end
    end
    return order
end

"""
Construct a set of matrices that encodes a Hamiltonian's term which resembles the hopping
from EBR `br1` to EBR `br2`. 

The Hamiltonian's order which is implicitly used is returned as output, and the matrices 
are stored on a 4D matrix which last two axes indicate the Hamiltonian term position it is 
describing and the first axis refer to the vector `δs` and the second axis to the vector `t`. 
See `devdocs.md` for details.
"""
function construct_M_matrix(
    δs::Vector{Pair{RVec{D},Vector{NTuple{3,RVec{D}}}}},
    br1::NewBandRep{D},
    br2::NewBandRep{D};
    order=hamiltonian_term_order(br1, br2) # an internal order for the Hamiltonian's terms
) where {D}
    V = length(δs)
    E = length(last(first(δs))) # number of elements in each δ_r (constant for all r)
    foreach(δs) do (_, δ_r)
        length(δ_r) == E || error("Unexpectedly had different counts of elements across δ_r")
    end
    Q1, Q2 = irdim(br1.siteir), irdim(br2.siteir)
    Q = Q1 * Q2

    # matrix of matrices that will store the matrices for each
    # Hamiltonian's term
    Mm = zeros(Int, V, V * E * Q, size(order, 1), size(order, 2))

    # fill in the unit-elements of `Mm`
    for α in axes(order, 1), β in axes(order, 2)
        (i, q), (j, w) = order[α, β]
        q = parent(q)
        w = parent(w)

        # assign a 1 to the correct position
        for (r, (_, δ_r)) in enumerate(δs)
            offset0 = (r - 1) * E * Q
            for (x, hop) in enumerate(δ_r)
                if hop[1] ≈ q && hop[2] ≈ w
                    offset1 = (x - 1) * Q
                    c = offset0 + offset1 + (j - 1) * Q1 + i
                    Mm[r, c, α, β] = 1
                end
            end
        end
    end

    return Mm
end

# H_{s,t} = v_i M_{i,j,s,t} t_j

# build the Q matrix for a particular symmetry operation (or, equivalently, a particular
# matrix from the site-symmetry representation), acting on the M matrix. Relative to our
# white-board notes, Q has swapped indices, in the sense we below give Q[i,j,r,f].
function representation_constraint_matrices(
    Mm::Array{Int,4},
    ρ_αα::AbstractMatrix{<:Number},
    ρ_ββ::AbstractMatrix{<:Number}
)

    ρ_αα = Matrix(ρ_αα) # since `/` doesn't extend to BlockArrays currently
    ρ_ββ = Matrix(ρ_ββ)

    Q = zeros(ComplexF64, size(Mm))
    for i in axes(Mm, 1), j in axes(Mm, 2)
        Q[i, j, :, :] .= ρ_αα * Mm[i, j, :, :] / ρ_ββ # = ρ_αα * Mm[i, j, :, :] * inv(ρ_ββ)
    end
    return Q
end

function constraint_matrices(
    br_α::NewBandRep{D},
    br_β::NewBandRep{D},
    δs,
    order=hamiltonian_term_order(br1, br2)
) where {D}
    # We obtain the needed representations over the generators of each bandrep
    gens_α, ρs_αα = sgrep_induced_by_siteir_generators(br_α)
    gens_β, ρs_ββ = sgrep_induced_by_siteir_generators(br_β)
    @assert gens_α == gens_β # must be from same space group and in same sorting

    # cast generators to primitive basis
    cntr = centering(num(br_α), D)
    gens = cntr ∈ ('P', 'p') ? gens_α : primitivize.(gens_α, cntr)

    # compute the tensor M that encodes the Hamiltonian as a numerical matrix
    Mm = construct_M_matrix(δs, br_α, br_β, order)

    # compute the Q tensor, encoding representation constraints on H_αβ
    Qs = Vector{Array{Complex,4}}(undef, length(gens))
    for (i, (ρ_αα, ρ_ββ)) in enumerate(zip(ρs_αα, ρs_ββ))
        Qs[i] = representation_constraint_matrices(Mm, ρ_αα, ρ_ββ)
    end

    # compute the Z tensor, encoding reciprocal-rotation constraints on H_αβ
    Zs = reciprocal_contraints_matrices(Mm, gens, δs)

    # build an aggregate constraint matrix, over all generators, acting on the hopping
    # coefficient vector t_αβ associated with δs
    constraint_ms = Vector{Matrix{ComplexF64}}()
    for (Q, Z) in zip(Qs, Zs)
        for s in axes(Q, 3), t in axes(Q, 4)
            q = Q[:, :, s, t]
            z = Z[:, :, s, t]
            push!(constraint_ms, q - z)
        end
    end
    constraints = reduce(vcat, constraint_ms)
    t_αβ_basis_matrix_form = nullspace(constraints; atol=NULLSPACE_ATOL_DEFAULT)

    # convert null-space to a sparse column form
    t_αβ_basis_matrix_form′ = poormans_sparsification(t_αβ_basis_matrix_form)
    t_αβ_basis = [collect(v) for v in eachcol(t_αβ_basis_matrix_form′)]

    # prune near-zero elements of basis vectors
    prune_at_threshold!(t_αβ_basis)

    return Mm, t_αβ_basis, order
end

function reciprocal_contraints_matrices(
    Mm::Array{Int,4},
    gens::AbstractVector{<:SymOperation},
    δs
)
    Zs = Vector{Array{Int,4}}(undef, length(gens))
    δ_hops = first.(δs)
    for (i, op) in enumerate(gens)
        Z = zeros(ComplexF64, size(Mm))
        P = permute_symmetry_related_hoppings_under_symmetry_operation(δ_hops, op)
        Pᵀ = transpose(P)
        for l in axes(P, 1), j in axes(Mm, 2), s in axes(Mm, 3), t in axes(Mm, 4)
            Z[l, j, s, t] = sum(Pᵀ[l, i] * Mm[i, j, s, t] for i in axes(P, 2))
        end
        Zs[i] = Z
    end
    return Zs
end

function permute_symmetry_related_hoppings_under_symmetry_operation(
    δ_hops::Vector{RVec{D}},
    op::SymOperation
) where {D}
    P = zeros(Int, length(δ_hops), length(δ_hops))
    for (i, δ) in enumerate(δ_hops)
        δ′ = compose(op, δ)
        j = findfirst(≈(δ′), δ_hops)
        isnothing(j) && error(lazy"hopping element $δ not closed under $op in $δ_hops")
        P[i, j] = 1
        # P acts as g on v: g v = P v so (gv)ᵢ = vⱼ = ∑ₖ Pᵢₖ vₖ so Pᵢₖ = δᵢⱼ
    end
    return P
end

"""
Poor man's "matrix sparsification" via the reduced row echelon form.
"""
function poormans_sparsification(
    A;
    rref_tol::Union{Nothing,Float64}=SPARSIFICATION_ATOL_DEFAULT)
    # following appendix E of the Qsymm paper (https://arxiv.org/abs/1806.08363) [copied
    # over from Neumann.jl]
    if !isnothing(rref_tol)
        # use a relatively low tolerance in `rref` to avoid explosions of errors
        # NB: this optional tolerance argument of `rref!` is undocumented :(
        return transpose(rref!(copy(transpose(A)), rref_tol))
    end
    return transpose(rref(transpose(A)))
end

function prune_at_threshold!(
    vs::AbstractVector{<:AbstractVector{T}};
    atol::Real=PRUNE_ATOL_DEFAULT
) where {T<:Complex}

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

# ---------------------------------------------------------------------------------------- #

function tb_hamiltonian(
    cbr::CompositeBandRep{D},
    Rs::AbstractVector{Vector{Int}} # "global" translation-representatives of hoppings to consider
) where {D}
    if any(c -> !isinteger(c) || c < 0, cbr.coefs)
        error("provided composite bandrep is not Wannierizable: contains negative or noninteger coefficients")
    end
    coefs = round.(Int, cbr.coefs)

    brs = Vector{NewBandRep{D}}(undef, sum(coefs))
    idx = 0
    for (i, c) in enumerate(coefs)
        for _ in 1:c
            brs[idx+=1] = cbr.brs[i]
        end
    end
    # find all families of hoppings between involved band representations
    representative_δs = RVec{D}[]
    for br1 in brs
        for br2 in brs
            δss = obtain_symmetry_related_hoppings(Rs, br1, br2)
            for δ in keys(δss)
                if !isapproxin(δ, representative_δs)
                    push!(representative_δs, δ)
                end
            end
        end
    end

    Norbs = count_bandrep_orbitals.(brs)
    tbs = [BlockMatrix{TightBindingElementString,Matrix{TightBindingBlock{D}}}(
        undef_blocks, Norbs, Norbs) for _ in representative_δs]
    c_idx_start = 1
    for (block_i, br1) in enumerate(brs)
        # TODO: maybe only need to go over upper triangular part of loop cf. hermicity
        #       (br1 vs br2 ~ br2 vs. br1)?
        for (block_j, br2) in enumerate(brs)
            δss = obtain_symmetry_related_hoppings(Rs, br1, br2)
            seen_n = Set{Int}()
            order = hamiltonian_term_order(br1, br2)
            for (δ, δhops) in δss
                n = something(findfirst(≈(δ), representative_δs))
                push!(seen_n, n)
                Mm, t_αβ_basis, _ = constraint_matrices(br1, br2, δhops, order)

                A = TightBindingBlock{D}(
                    (block_i, block_j), (Norbs[block_i], Norbs[block_j]), br1, br2,
                    order, Mm, t_αβ_basis, δhops, c_idx_start:c_idx_start+length(t_αβ_basis))
                tbs[n][Block(block_i), Block(block_j)] = A
                c_idx_start += length(t_αβ_basis)
            end
            # blocks for other values of `n` are not featured in `δss` - i.e., vanish, so we
            # put in manually construct zero blocks for those spots:
            for n in eachindex(representative_δs)
                n ∈ seen_n && continue
                tbs[n][Block(block_i), Block(block_j)] = TightBindingBlock{D}(
                    (block_i, block_j), (Norbs[block_i], Norbs[block_j]), br1, br2,
                    order,
                    zeros(Int, 0, 0, size(order, 1), size(order, 2)),    # Mm
                    Vector{Vector{ComplexF64}}(),                        # t_αβ_basis
                    Vector{Pair{RVec{D},Vector{NTuple{3,RVec{D}}}}}(), # δs
                    1:0)                                                 # c_idxs (empty)
            end
        end
    end
    return tbs
end

struct TightBindingElementString
    s::String
end
Base.show(io::IO, tbe_str::TightBindingElementString) = print(io, tbe_str.s)

struct TightBindingBlock{D} <: AbstractMatrix{TightBindingElementString}
    block_ij::Tuple{Int,Int}
    global_ij::Tuple{Int,Int}
    br1::NewBandRep{D}
    br2::NewBandRep{D}
    order::Matrix{Pair{Tuple{Int64,WyckoffPosition{D}},Tuple{Int64,WyckoffPosition{D}}}}
    Mm::Array{Int,4}
    t_αβ_basis::Vector{Vector{ComplexF64}}
    δs::Vector{Pair{RVec{D},Vector{NTuple{3,RVec{D}}}}}
    c_idxs::UnitRange{Int}
    # TODO: figure out how much of this is needed, e.g.:
    # TODO: is it redundant to have both `order` and `δs`? Not completely I guess
    # TODO: find better schema for naming the free coefficients than the arbitrary `c_idxs`
    # TODO: do we need the indexing information in `block_ij` and `global_ij`? Probably not
end
Base.size(tbb::TightBindingBlock) = (size(tbb.Mm, 3), size(tbb.Mm, 4))
function Base.getindex(tbb::TightBindingBlock, i::Int, j::Int)
    exp_strs = Vector{String}(undef, length(tbb.δs))
    for (n, δ) in enumerate(tbb.δs)
        io_kr = IOBuffer()
        first_nonzero = true
        for (l, δₗ) in enumerate(δ[1].cnst)
            abs2δₗ = 2abs(δₗ)
            abs2δₗ < SPARSIFICATION_ATOL_DEFAULT && continue
            if δₗ < 0 || !first_nonzero
                print(io_kr, Crystalline.signaschar(δₗ))
            end
            first_nonzero = false
            str_enum, v_r, v_str = _stringify_characters(abs2δₗ)
            str_enum == REAL_STR || error("unexpected imaginary component in exponential argument $abs2δₗ")
            isone(v_r) || print(io_kr, v_str)
            print(io_kr, "k", Crystalline.subscriptify(string(l)))
        end
        exp_arg = String(take!(io_kr))
        if !isempty(exp_arg)
            exp_strs[n] = "𝕖(" * exp_arg * ")" # short-hand: 𝕖(x) = exp(iπx)
        else
            exp_strs[n] = "" # = 1, but omit for compactness
        end
    end
    exp_strs

    Mm = tbb.Mm

    io = IOBuffer()
    first_t_αβ_basis_vec = true
    for k in eachindex(tbb.t_αβ_basis)
        Mⁱʲtᵏ = Mm[:, :, i, j] * tbb.t_αβ_basis[k]
        nnz_els = count(v -> abs(v) > SPARSIFICATION_ATOL_DEFAULT, Mⁱʲtᵏ)
        nnz_els == 0 && continue
        first_t_αβ_basis_vec || (first_t_αβ_basis_vec = false; print(io, " + "))
        print(io, "c", Crystalline.subscriptify(string(tbb.c_idxs[k])))
        nnz_els > 1 && print(io, "[")
        first = true
        for (n, v) in enumerate(Mⁱʲtᵏ)
            abs(v) < SPARSIFICATION_ATOL_DEFAULT && continue
            str_enum, v_r, v_str = _stringify_characters(v)
            if str_enum == REAL_STR || str_enum == IMAG_STR
                if isone(abs(v_r))
                    v_str = ""
                else
                    v_str = lstrip(v_str, ('-', '+'))
                end
                if first
                    v_str = (v_r < 0 ? "-" : "") * v_str
                else
                    v_str = (v_r < 0 ? " - " : " + ") * v_str
                end
            else
                v_str = (first ? "" : " + ") * "(" * v_str * ")"
            end
            print(io, v_str, exp_strs[n])
            first = false
        end
        nnz_els > 1 && print(io, "]")
    end
    s = String(take!(io))
    return TightBindingElementString(s)
end
Base.setindex!(::TightBindingBlock, v, ij...) = error("setindex! is not supported")


function count_bandrep_orbitals(br::NewBandRep{D}) where {D}
    mult = multiplicity(position(br)) # multiplicity in conventional setting
    # we need the Wyckoff multiplicity, excluding conventional-centering copies, so we
    # divide by the the number of centering-translations
    cntr = centering(num(br), D)
    denom = centering_volume_fraction(cntr, Val(D), Val(D))
    mult = div(mult, denom)

    return mult * irdim(br.siteir)
end

# ---------------------------------------------------------------------------------------- #
# pretty-printing of scalars

@enum StrPrintAs begin
    REAL_STR
    IMAG_STR
    COMPLEX_STR
end
function _stringify_characters(c::Number; digits::Int=3)
    c′ = round(c; digits)
    cr, ci = reim(c′)
    if iszero(ci)     # real
        isinteger(cr) && return (REAL_STR, cr, string(Int(cr)))
        return (REAL_STR, cr, string(cr))

    elseif iszero(cr) # imaginary
        isinteger(ci) && return (IMAG_STR, ci, string(Int(ci)) * "i")
        return (IMAG_STR, ci, string(ci) * "i")

    else              # complex
        if isinteger(cr) && isinteger(ci)
            return COMPLEX_STR, zero(cr), _complex_as_compact_string(Complex{Int}(cr, ci))
        else
            return COMPLEX_STR, zero(cr), _complex_as_compact_string(c′)
        end
    end
end
function _complex_as_compact_string(c::Complex) # usual string(::Complex) has spaces; avoid that
    io = IOBuffer()
    print(io, real(c), Crystalline.signaschar(imag(c)), abs(imag(c)), "i")
    return String(take!(io))
end