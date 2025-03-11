# TODO: add all of this to the previous functions as an `if` statement

# Similar approach to TRS but simpler
# The condition that we want to impose now is H(k) = (Hᵀ(k))*, then we have that
# (H*(k))ₛₜ = (Ρv)ᵀ (Mₛₜ| -Mₛₜ) (tᴿ,tᴵ) => ((H*(k))ᵀ)ₛₜ = (Ρv)ᵀ (Mₜₛ| -Mₜₛ) (tᴿ,tᴵ)
# where Ρ accounts for the complex conjugation on v.


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
function obtain_basis_free_parameters_hermiticity(
    brₐ::NewBandRep{D},
    brᵦ::NewBandRep{D},
    h_orbit::HoppingOrbit{D},
    order=hamiltonian_term_order(brₐ, brᵦ)
) where {D}
    size(order, 1) == size(order, 2) || error("The Hamiltonian must be a square matrix")

    # compute the tensor M that encodes the Hamiltonian as a numerical matrix
    Mm = construct_M_matrix(h_orbit, brₐ, brᵦ, order)

    # now we need to add hermiticity constraints. For this we first duplicate 
    # the M matrix since H = vᵢᵀ (Mᵢⱼ|Mᵢⱼ) (tᴿⱼ,tᴵⱼ)
    # WARNING : we assume that we have ±δ in v. Is this true? why? it is ok physically
    #           and in here we will see if they are equal or not.

    # compute the Z tensor, encoding the hermiticity constraint on H for the k-space
    # part. Since it is just H(k), we just need to double Mm to account for (tᴿ,tᴵ)

    Z_hermiticity = [Mm Mm]

    # compute the Q tensor, encoding hermiticity constraints on H for the free-
    # parameter part. This is done by ((H*(k))ᵀ)ₛₜ = (Ρv)ᵀ (Mₜₛ| -Mₜₛ) (tᴿ,tᴵ) =
    # where Ρ accounts for the complex conjugation on v.

    Q_hermiticity = representation_constraint_hermiticity(Mm, h_orbit)

    # build an constraint matrix acting on the hopping coefficient vector tₐᵦ 
    #associated with h_orbit
    constraint_vs = Vector{Matrix{ComplexF64}}()
    for s in axes(Q_hermiticity, 3), t in axes(Q_hermiticity, 4)
        q = @view Q_hermiticity[:, :, s, t]
        z = @view Z_hermiticity[:, :, s, t]
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
    representation_constraint_hermiticity(Mm::Array{Int,4}, h_orbit::HoppingOrbit{D})
    --> Array{ComplexF64,4}

If hermiticity is present, we need to add the constraint. It is given
by the association δ -> -δ and the complex conjugation in the free-parameter part
tⱼ -> tⱼ* = (tⱼᴿ|tⱼᴵ) -> (tⱼᴿ| -tⱼᴵ), and a transposition in the `Mm` matrix.
"""
function representation_constraint_hermiticity(
    Mm::AbstractArray{<:Number,4},
    h_orbit::HoppingOrbit{D},
) where {D}
    Q_hermiticity = zeros(ComplexF64, size(Mm))
    op = SymOperation([-I(D) zeros(D)])
    P = _permute_symmetry_related_hoppings_under_symmetry_operation(h_orbit, op)
    for l in axes(P, 2), j in axes(Mm, 2), s in axes(Mm, 3), t in axes(Mm, 4)
        Q_hermiticity[l, j, t, s] = sum(P[i, l] * Mm[i, j, s, t] for i in axes(P, 1))
    end
    return [Q_hermiticity -Q_hermiticity]
end