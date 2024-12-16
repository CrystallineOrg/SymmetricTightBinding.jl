using LinearAlgebra: nullspace

"""
Compute the symmetry related hoppings relative vectors from the WP of `br1` to the WP of `br2`
displaced a set of primitive lattice vectors `Rs`.
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

    δsdd = Dict{RVec{D},Vector{Pair{RVec{D},Vector{Tuple{RVec{D},RVec{D},RVec{D}}}}}}()
    for R in Rs
        R = RVec(R) # change the type of R to be type consistent
        for (qₐ, qᵦ) in Iterators.product(wps1, wps2)
            qₐ = parent(qₐ) # to work with the RVec type directly
            qᵦ = parent(qᵦ) # same as above
            δ = qᵦ + R - qₐ
            if !in(δ, keys(δsdd)) && !in(δ, Iterators.flatten(first.(values(δsdd))))
                push!(δsdd, δ => [δ => [(qₐ, qᵦ, R)]])
                for g in ops
                    Ρ = SymOperation(rotation(g)) # type consitentency for the rotation 
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
                    idx = findfirst(v -> first(v) == δ′, δsdd[δ])
                    if !in(δ′, keys(δsdd)) && isnothing(idx)
                        push!(δsdd[δ], δ′ => [(qₐ′, qᵦ′, R′)])
                    elseif !isnothing(idx) && !in((qₐ′, qᵦ′, R′), last(δsdd[δ][something(idx)]))
                        push!(last(δsdd[δ][something(idx)]), (qₐ′, qᵦ′, R′))
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
    δs::Vector{Pair{RVec{D},Vector{Tuple{RVec{D},RVec{D},RVec{D}}}}},
    br1::NewBandRep{D},
    br2::NewBandRep{D},
) where {D}

    V = length(δs)
    E = length(last(first(δs))) # number of elements in each δ_r (constant for all r)
    foreach(δs) do (_, δ_r)
        length(δ_r) == E || error("Unexpected had different counts of elements across δ_r")
    end
    Q1, Q2 = irdim(br1.siteir), irdim(br2.siteir)
    Q = Q1 * Q2

    # we need to declare an order that the Hamiltonian's term internally has
    order = hamiltonian_term_order(br1, br2)

    # matrix of matrices that will store the matrices for each
    # Hamiltonian's term
    Mm = zeros(Int, V, V * E * Q, size(order)[1], size(order)[2])

    for α in axes(order, 1), β in axes(order, 2)
        ((i, q), (j, w)) = order[α, β]
        q = parent(q)
        w = parent(w)

        ## Introduce code that assigns a 1 to the correct position
        for (r, (_, δ_r)) in enumerate(δs)
            offset0 = (r - 1) * E * Q
            for (x, hop) in enumerate(δ_r)
                if hop[1] == q && hop[2] == w
                    offset1 = (x - 1) * Q
                    c = offset0 + offset1 + i + j - 1
                    Mm[r, c, α, β] = 1
                end
            end
        end
    end

    return order, Mm
end

# H_{s,t} = v_i M_{i,j,s,t} t_j

"""
Build the Q matrix for a particular symmetry operation (or, equivalently, a particular
matrix from the site-symmetry representation), acting on the M matrix. Relative to our
white-board notes, Q has swapped indices, in the sense we below give Q[i,j,r,f].
"""
function representation_constraint_matrices(
    Mm::AbstractArray{<:Number,4},
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

"""
Find the null-space of the constraints impose by symmetries on the Hamiltonian, which 
translate to Zs = Qs, and returns a basis for the null-space.
"""
function symmetry_constraint_solution(
    br_α::NewBandRep,
    br_β::NewBandRep,
    δs::Vector{Pair{RVec{D},Vector{Tuple{RVec{D},RVec{D},RVec{D}}}}},
) where {D}
    # We obtain the needed representations over the generators of each bandrep
    gens_and_ρs_αα = sgrep_induced_by_siteir_generators(br_α)
    gens_and_ρs_ββ = sgrep_induced_by_siteir_generators(br_β)

    # compute the tensor M that encodes the Hamiltonian as a numerical matrix
    or, Mm = construct_M_matrix(δs, br_α, br_β)

    # compute the Q tensor, encoding representation constraints on H_αβ
    Qs = Vector{Array{Complex,4}}(undef, length(gens_and_ρs_αα))
    gens = Vector{SymOperation}(undef, length(gens_and_ρs_αα))
    for (i, (gen_and_ρ_αα, gen_and_ρ_ββ)) in enumerate(zip(gens_and_ρs_αα, gens_and_ρs_ββ))
        gen_α, gen_β = first(gen_and_ρ_αα), first(gen_and_ρ_ββ)
        gen_α == gen_β || error(lazy"bandrep generators differ; they should not")
        gens[i] = gen_α

        ρ_αα, ρ_ββ = last(gen_and_ρ_αα), last(gen_and_ρ_ββ)
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

    return or, Mm, t_αβ_basis
end

"""
Build the Z matrix for a particular symmetry operation acting on k-space over the M matrix.
"""
function reciprocal_contraints_matrices(
    Mm::AbstractArray{<:Number,4},
    gens::Vector{SymOperation},
    δs::Vector{Pair{RVec{D},Vector{Tuple{RVec{D},RVec{D},RVec{D}}}}}
) where {D}
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

"""
Build the P matrix for a particular symmetry operation acting on k-space, which permutes 
the rows of the M matrix. 
"""
function permute_symmetry_related_hoppings_under_symmetry_operation(
    δ_hops::Vector{RVec{D}}, op::SymOperation
) where {D}
    P = zeros(Int, length(δ_hops), length(δ_hops))
    for (i, δ) in enumerate(δ_hops)
        δ′ = compose(op, δ)
        j = findfirst(==(δ′), δ_hops)
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
    A::AbstractMatrix{<:Number};
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

"""
Prune near-zero elements of vectors in `vs`.
"""
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