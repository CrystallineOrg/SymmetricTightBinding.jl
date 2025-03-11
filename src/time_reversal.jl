# TODO: improve the code so it works for COMPLEX and PSEUDOREAL irreps

# we are going to assume for now that the transform as real irreps of the site-symmetry 
# group. TRS can be understood as a spacial symmetry when acting on the Hamiltonian:
# D(𝒯)H(k)D(𝒯)⁻¹ = H(𝒯k) -> Γ(𝒯)H*(k)Γ(𝒯)⁻¹ = H(-k), where D is the whole 
# operator and Γ is only the unitary part, so D(𝒯) = Γ(𝒯)𝒯
# If the site-symmetry irrep is real, Γ(𝒯) = I -> H*(k) = H(-k).

# below i have assumed that the unitary part Γ doesn't have any complex entry.
# if a complex entry is present, we need to perform an extra permutation on 
# (Mᵢⱼ| -Mᵢⱼ). Imagine we have an operation: 
# (ΓH*(k)Γ⁻¹)ₛₜ = vᵀ Pᵀ [Γₛₙ (Mₙᵣ| -Mₙᵣ) Γᵣₜ⁻¹] (tᴿ,tᴵ)
#
# for simplicity i am going to pick Γ = [0 i; -i 0] (complex entries in the 
# representation). Then, i can build an auxiliary matrix 𝒢 such that 
# 𝒢(Γ₁₂) = [0 1; -1 0] = -𝒢(Γ₂₁⁻¹) and the other being zero. Then, we have that 
# the permutation is given by 
# vᵀ Pᵀ [Γ₁₂ (Mₙᵣ| -Mₙᵣ) Γ₁₂⁻¹] (tᴿ,tᴵ) = vᵀ Pᵀ [(Mₙᵣ| -Mₙᵣ)𝒢(Γ₁₂)𝒢(Γ₁₂⁻¹)] (tᴿ,tᴵ)
# This is based on the fact that matrices acting on the right can permute columns.
# TODO: is it better to expand by rows instead of columns? 
#       (Mₙᵣ|Mₙᵣ)(tᴿ,tᴵ) -> (Mₙᵣ,Mₙᵣ)(tᴿ|tᴵ)
# WARNING: does something similar apply also to crystalline symmetries?
#          check example SG213 EBR (4a|E) -> yes I included function `split_complex`
#          (maybe unnecessary as standalone function)

"""
    obtain_basis_free_parameters_TRS(
                                    brₐ::NewBandRep{D}, 
                                    brᵦ::NewBandRep{D}, 
                                    h_orbit::HoppingOrbit{D}, 
                                    order=hamiltonian_term_order(brₐ, brᵦ)
                                    ) --> Tuple{Array{Int,4}, Vector{Vector{ComplexF64}}, 
                                                Matrix{Pair{Tuple{Int,WyckoffPosition{D}},
                                                            Tuple{Int,WyckoffPosition{D}}}}

Obtain the basis of free parameters for the hopping terms between `brₐ` and `brᵦ` 
associated with the hopping orbit `h_orbit`. The Hamiltonian's default order is 
given by `order`. The constraints assume also time-reversal symmetry. A 
differentiation between real and imaginary components is performed.
"""
function obtain_basis_free_parameters_TRS(
    brₐ::NewBandRep{D},
    brᵦ::NewBandRep{D},
    h_orbit::HoppingOrbit{D},
    order=hamiltonian_term_order(brₐ, brᵦ)
) where {D}
    # compute the tensor M that encodes the Hamiltonian as a numerical matrix
    Mm = construct_M_matrix(h_orbit, brₐ, brᵦ, order)

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
    # parameter part. This is done by H*(k) = vᵢᵀ* (Mᵢⱼ| -Mᵢⱼ) (tᴿⱼ,tᴵⱼ) = 
    # (Pvᵢ)ᵀ (Mᵢⱼ| -Mᵢⱼ) (tᴿⱼ,tᴵⱼ) = vᵢᵀ (Pᵀ Mᵢⱼ| -Pᵀ Mᵢⱼ) (tᴿⱼ,tᴵⱼ)

    Q_trs = representation_constraint_trs(Mm, h_orbit)

    # build an constraint matrix acting on the hopping coefficient vector tₐᵦ 
    #associated with h_orbit
    constraint_vs = Vector{Matrix{ComplexF64}}()
    for s in axes(Q_trs, 3), t in axes(Q_trs, 4)
        q = @view Q_trs[:, :, s, t]
        z = @view Z_trs[:, :, s, t]
        c = q - z
        filtered_rows = filter(r -> norm(r) > 1e-10, eachrow(c))
        isempty(filtered_rows) && continue # don't add empty constraints
        append!(constraint_vs, filtered_rows)
    end

    constraints = stack(constraint_vs, dims=1)
    tₐᵦ_basis_matrix_form = nullspace(constraints; atol=NULLSPACE_ATOL_DEFAULT)

    # convert null-space to a sparse column form
    tₐᵦ_basis_matrix_form′ = _poormans_sparsification(tₐᵦ_basis_matrix_form)
    tₐᵦ_basis = [collect(v) for v in eachcol(tₐᵦ_basis_matrix_form′)]
    # TODO: maybe we can keep it as matrix

    # prune near-zero elements of basis vectors
    _prune_at_threshold!(tₐᵦ_basis)

    return [Mm Mm], tₐᵦ_basis, order
end

"""
    reciprocal_constraints_trs(Mm::Array{Int,4}, h_orbit::HoppingOrbit{D}) 
    --> Array{ComplexF64,4}

If time reversal symmetry is present, we need to add the constraint. It is given
by the association k -> -k.
"""
function reciprocal_constraints_trs(
    Mm::Array{Int,4},
    h_orbit::HoppingOrbit{D}
) where {D}
    Z_trs = zeros(ComplexF64, size(Mm))
    op = SymOperation([-I(D) zeros(D)])
    P = _permute_symmetry_related_hoppings_under_symmetry_operation(h_orbit, op)
    for l in axes(P, 2), j in axes(Mm, 2), s in axes(Mm, 3), t in axes(Mm, 4)
        Z_trs[l, j, s, t] = sum(P[i, l] * Mm[i, j, s, t] for i in axes(P, 1))
    end
    return [Z_trs Z_trs]
end

"""
    representation_constraint_trs(Mm::Array{Int,4}, h_orbit::HoppingOrbit{D})
    --> Array{ComplexF64,4}

If time reversal symmetry is present, we need to add the constraint. It is given
by the association δ -> -δ and the complex conjugation in the free-parameter part
tⱼ -> tⱼ* = (tⱼᴿ|tⱼᴵ) -> (tⱼᴿ| -tⱼᴵ).
"""
function representation_constraint_trs(
    Mm::AbstractArray{<:Number,4},
    h_orbit::HoppingOrbit{D},
) where {D}
    Q_trs = zeros(ComplexF64, size(Mm))
    op = SymOperation([-I(D) zeros(D)])
    P = _permute_symmetry_related_hoppings_under_symmetry_operation(h_orbit, op)
    for l in axes(P, 2), j in axes(Mm, 2), s in axes(Mm, 3), t in axes(Mm, 4)
        Q_trs[l, j, s, t] = sum(P[i, l] * Mm[i, j, s, t] for i in axes(P, 1))
    end
    return [Q_trs -Q_trs]
end