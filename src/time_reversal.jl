# TODO: add all of this as if in the previous functions
# TODO: improve the code so it works for COMPLEX and PSEUDOREAL irreps

# we are going to assume for now that the transform as real irreps of the site-symmetry 
# group. TRS can be understood as a spacial symmetry when acting on the Hamiltonian:
# D(𝒯)H(k)D(𝒯)⁻¹ = H(𝒯k) -> Γ(𝒯)H*(k)Γ(𝒯)⁻¹ = H(-k), where D is the whole 
# operator and Γ is only the unitary part, so D(𝒯) = Γ(𝒯)𝒯
# If the site-symmetry irrep is real, Γ(𝒯) = I -> H*(k) = H(-k).


"""
    obtain_basis_free_parameters_TRS(
                                    brₐ::NewBandRep{D}, 
                                    brᵦ::NewBandRep{D}, 
                                    h_orbit::HoppingOrbit{D}, 
                                    order=hamiltonian_term_order(brₐ, brᵦ)
                                    ) --> Tuple{Array{Int,4}, Vector{Vector{ComplexF64}}, 
                                                Matrix{Pair{Tuple{Int,WyckoffPosition{D}},
                                                            Tuple{Int,WyckoffPosition{D}}}}

Obtain the basis of free parameters for the hopping terms between `brₐ` and `brᵦ` associated
with the hopping orbit `h_orbit`. The Hamiltonian's default order is given by `order`.
The constraints assume also time-reversal symmetry. A differentiation between 
real and imaginary components is performed.
"""
function obtain_basis_free_parameters_TRS(
    brₐ::NewBandRep{D},
    brᵦ::NewBandRep{D},
    h_orbit::HoppingOrbit{D},
    order=hamiltonian_term_order(brₐ, brᵦ)
) where {D}
    # We obtain the needed representations over the generators of each bandrep
    gensₐ = generators(num(brₐ), SpaceGroup{D})
    gensᵦ = generators(num(brᵦ), SpaceGroup{D})
    @assert gensₐ == gensᵦ # must be from same space group and in same sorting

    # cast generators to primitive basis
    cntr = centering(num(brₐ), D)
    gens = cntr ∈ ('P', 'p') ? gensₐ : primitivize.(gensₐ, cntr)

    # compute the tensor M that encodes the Hamiltonian as a numerical matrix
    Mm = construct_M_matrix(h_orbit, brₐ, brᵦ, order)

    #compute the Q tensor, encoding representation constraints on Hₐᵦ
    Qs = representation_constraint_matrices(Mm, brₐ, brᵦ)

    # compute the Z tensor, encoding reciprocal-rotation constraints on Hₐᵦ
    Zs = reciprocal_constraints_matrices(Mm, gens, h_orbit)

    # now we need to add time-reversal constraints. For this we first duplicate 
    # the M matrix since H = vᵢᵀ (Mᵢⱼ|Mᵢⱼ) (tᴿⱼ,tᴵⱼ)
    # WARNING : we assume that we have ±δ in v. Is this true? why? it is ok physically
    #           and in here we will see if they are equal or not.

    brₐ.siteir.iscorep == brᵦ.siteir.iscorep == false || error("Not implemented for COMPLEX or PSEUDOREAL irreps")

    # compute the Z tensor, encoding time-reversal constraints on H for the k-space
    # part. This is done by H(-k) = (Ρv)ᵢᵀ (Mᵢⱼ|Mᵢⱼ) (tᴿⱼ,tᴵⱼ) = vᵢᵀ Ρᵀ (Mᵢⱼ|Mᵢⱼ) (tᴿⱼ,tᴵⱼ)
    # vᵢᵀ (Ρᵀ Mᵢⱼ|Ρᵀ Mᵢⱼ) (tᴿⱼ,tᴵⱼ)

    Z_trs = reciprocal_constraints_trs(Mm, h_orbit)

    # compute the Q tensor, encoding time-reversal constraints on H for the free-
    # parameter part. This is done by H*(k) = vᵢᵀ* (Mᵢⱼ|-Mᵢⱼ) (tᴿⱼ,tᴵⱼ) = 
    # (Pvᵢ)ᵀ (Mᵢⱼ|-Mᵢⱼ) (tᴿⱼ,tᴵⱼ) = vᵢᵀ (Pᵀ Mᵢⱼ|-Pᵀ Mᵢⱼ) (tᴿⱼ,tᴵⱼ)

    Q_trs = representation_constraint_trs(Mm, h_orbit)

    Zs_trs = Vector{Array{ComplexF64,8}}(undef, length(Zs) + 1)
    Qs_trs = Vector{Array{ComplexF64,8}}(undef, length(Qs) + 1)

    for i in axes(Zs, 2)
        Zs_trs[i] = [Zs[i] Zs[i]]
        Qs_trs[i] = [Qs[i] Qs[i]]
    end

    Zs_trs[end] = Z_trs
    Qs_trs[end] = Q_trs

    # build an aggregate constraint matrix, over all generators, acting on the hopping
    # coefficient vector tₐᵦ associated with h_orbit
    constraint_ms = Vector{Matrix{ComplexF64}}()
    for (Q, Z) in zip(Qs_trs, Zs_trs)
        for s in axes(Q, 3), t in axes(Q, 4)
            q = Q[:, :, s, t]
            z = Z[:, :, s, t]
            push!(constraint_ms, q - z)
        end
    end
    constraints = reduce(vcat, constraint_ms)
    tₐᵦ_basis_matrix_form = nullspace(constraints; atol=NULLSPACE_ATOL_DEFAULT)

    # convert null-space to a sparse column form
    tₐᵦ_basis_matrix_form′ = _poormans_sparsification(tₐᵦ_basis_matrix_form)
    tₐᵦ_basis = [collect(v) for v in eachcol(tₐᵦ_basis_matrix_form′)]

    # prune near-zero elements of basis vectors
    _prune_at_threshold!(tₐᵦ_basis)

    return [Mm Mm], tₐᵦ_basis, order
end

"""
If time reversal symmetry is present, we need to add the constraint. It is given by
the association k -> -k.
"""
function reciprocal_constraints_trs(
    Mm::Array{Int,4},
    h_orbit::HoppingOrbit{D}
) where {D}
    Z_trs = zeros(ComplexF64, size(Mm))
    op = SymOperation([-I(D) zeros(D)])
    P = _permute_symmetry_related_hoppings_under_symmetry_operation(h_orbit, op)
    Pᵀ = transpose(P)
    for l in axes(P, 1), j in axes(Mm, 2), s in axes(Mm, 3), t in axes(Mm, 4)
        Z_trs[l, j, s, t] = sum(Pᵀ[l, i] * Mm[i, j, s, t] for i in axes(P, 2))
    end
    return [Z_trs Z_trs]
end

"""
If time reversal symmetry is present, we need to add the constraint. It is given by
the association δ -> -δ and the complex conjugation in the free-parameter part
tⱼ -> tⱼ* = (tⱼᴿ|tⱼᴵ) -> (tⱼᴿ|-tⱼᴵ).
"""
function representation_constraint_trs(
    Mm::AbstractArray{<:Number,4},
    h_orbit::HoppingOrbit{D},
) where {D}
    Q_trs = zeros(ComplexF64, size(Mm))
    op = SymOperation([-I(D) zeros(D)])
    P = _permute_symmetry_related_hoppings_under_symmetry_operation(h_orbit, op)
    Pᵀ = transpose(P)
    for l in axes(P, 1), j in axes(Mm, 2), s in axes(Mm, 3), t in axes(Mm, 4)
        Q_trs[l, j, s, t] = sum(Pᵀ[l, i] * Mm[i, j, s, t] for i in axes(P, 2))
    end
    return [Q_trs -Q_trs]
end