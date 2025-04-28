# Now we need to take care of what is real or complex to properly apply TRS.

# the vector `t` could be multiplied by a complex number. We want now to make 
# everything real. For this, we do the following extra step to the one above:
# `αt = (a + im*b)t = at + b (im*t) = [a; b] [t im*t]`
# where `α∈ℂ` and `a,b∈ℝ`. This way we duplicate the dimensionality of the problem,
# but we keep everything real.

# NOTE: this is implemented in the code just by expanding the basis of the nullspace 
# from `{tᵢ}ᵢ -> {tᵢ, im*tᵢ}ᵢ`, so for each complex variable `α` we now have two 
# real variables `a` and `b`.

# but in that way, `t∈ℂ`. If we want it to be real also, we can separated it into 
# its real and imaginary part in the following way: `t -> [real(t); imag(t)] = 
# [tᴿ; tᴵ].`

# so joining everything together, we have: αt -> [a; b] [t im*t] -> 
# [a; b] [real(t) real(im*t); imag(t) imag(im*t)]

# again note that we are just expanding the basis of the nullspace in the following 
# way `{tᵢ}ᵢ -> {tᵢ, im*tᵢ}ᵢ`, so for each complex variable `α` we now have two 
# real variables `a` and `b`.

# now, what happens with the irreps assuming that we are working with real numbers 
# only?

# we are going to assume for now that the transform as real irreps of the site-symmetry 
# group. TRS can be understood as a spacial symmetry when acting on the Hamiltonian:
# `D(𝒯)H(k)D(𝒯)⁻¹ = H(𝒯k) -> Γ(𝒯)H*(k)Γ(𝒯)⁻¹ = H(-k)`, where `D` is the whole 
# operator and `Γ` is only the unitary part, so `D(𝒯) = Γ(𝒯)𝒯`
# If the site-symmetry irrep is real, `Γ(𝒯) = I -> H*(k) = H(-k)`.

# we have made the following decomposition: `t -> [real(t); imag(t)]`.
# how does this affect the Hamiltonian representation?

# we have that `Hᵢⱼ(k) = vᵀ(k) Mᵢⱼ t`. Then, we need to make the following change:
# `Hᵢⱼ(k) = vᵀ(k) [Mᵢⱼ Mᵢⱼ] [real(t); imag(t)]`

# Then applying TRS will be:

# 1. `H(-k) = vᵀ(-k) [Mᵢⱼ Mᵢⱼ] [real(t); imag(t)] = (Pv)ᵀ(k) [Mᵢⱼ Mᵢⱼ] 
# [real(t); imag(t)] = vᵀ(k) [Pᵀ Mᵢⱼ Pᵀ Mᵢⱼ] [real(t); imag(t)]`

# 2. `H*(k) = (v*)ᵀ(k) [Mᵢⱼ Mᵢⱼ] [real(t); -imag(t) -imag(im*t)] = (Pv)ᵀ(k) 
# [Mᵢⱼ -Mᵢⱼ] [real(t); imag(t)] = vᵀ(k) [Pᵀ Mᵢⱼ -Pᵀ Mᵢⱼ] [real(t); imag(t)]`

"""
    obtain_basis_free_parameters_TRS(
        h_orbit::HoppingOrbit{D}, 
        brₐ::NewBandRep{D}, 
        brᵦ::NewBandRep{D}, 
        orderingₐ::OrbitalOrdering{D} = OrbitalOrdering(brₐ),
        orderingᵦ::OrbitalOrdering{D} = OrbitalOrdering(brᵦ),
        Mm::Array{4, Int} = construct_M_matrix(h_orbit, brₐ, brᵦ, orderingₐ, orderingᵦ)
        )                             --> Tuple{Array{Int,4}, Vector{Vector{ComplexF64}}}}

Obtain the basis of free parameters for the hopping terms between `brₐ` and `brᵦ` 
associated with the hopping orbit `h_orbit` under time-reversal symmetry.

Real and imaginary parts of the basis vectors are differentiated explicitly: internally,
we consider only variables.
"""
function obtain_basis_free_parameters_TRS(
    h_orbit::HoppingOrbit{D},
    brₐ::NewBandRep{D},
    brᵦ::NewBandRep{D},
    orderingₐ::OrbitalOrdering{D} = OrbitalOrdering(brₐ),
    orderingᵦ::OrbitalOrdering{D} = OrbitalOrdering(brᵦ),
    Mm::Array{Int, 4} = construct_M_matrix(h_orbit, brₐ, brᵦ, orderingₐ, orderingᵦ)
) where {D}
    # To add time-reversal constraints, we duplicate the `M` matrix so that we can transfer
    # the conjugation from `t` to `M`. We exploit the following splitting:
    #    `Hᵢⱼ = vᵀ Mᵢⱼ t = vᵀ Mᵢⱼ (tᴿ + itᴵ) = vᵀ [Mᵢⱼ Mᵢⱼ] [tᴿ; itᴵ]`

    # Step 1: compute the Z tensor, encoding time-reversal constraints on H for the k-space
    # part. This is done by `Hᵢⱼ(-k) = vᵀ(-k) [Mᵢⱼ Mᵢⱼ] [tᴿ; tᴵ]`
    # `= (Pv)ᵀ(k) [Mᵢⱼ Mᵢⱼ] [tᴿ; tᴵ]`
    # `= vᵀ(k) [Pᵀ Mᵢⱼ Pᵀ Mᵢⱼ] [tᴿ; tᴵ]`
    Z_trs = reciprocal_constraints_trs(Mm, h_orbit)
    # QUESTION: in the above assume that we have ±δ in v. Is this true? why? it is ok
    #           physically and in here we will see if they are equal or not.

    # Step 2: compute the Q tensor, encoding time-reversal constraints on H for the free-
    # parameter part. This is done by `H*(k) = (v*)ᵀ(k) [Mᵢⱼ Mᵢⱼ] [tᴿ; -itᴵ]`
    # `= (Pv)ᵀ(k) [Mᵢⱼ -Mᵢⱼ] [tᴿ; tᴵ]`
    # `= vᵀ(k) [Pᵀ Mᵢⱼ -Pᵀ Mᵢⱼ] [tᴿ; tᴵ]`
    Q_trs = representation_constraint_trs(Mm, h_orbit)

    # Step 3: build an constraint matrix acting on the doubled hopping coefficient vector
    # `[tᴿ; itᴵ]` associated with `h_orbit`
    constraint_vs = Vector{Vector{Float64}}()
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
    # TODO: are this vectors real? -> I think so but let's check it for now
    all(isreal.(tₐᵦ_basis)) || error("TRS constraint returns COMPLEX basis. Case not considered.")

    # prune near-zero elements of basis vectors
    _prune_at_threshold!(tₐᵦ_basis)

    return [Mm Mm], tₐᵦ_basis
end

"""
    reciprocal_constraints_trs(Mm::Array{Int,4}, h_orbit::HoppingOrbit{D}) 
    --> Array{ComplexF64,4}

Time reversal symmetry action on reciprocal space. It is given by the association 
`k -> -k => H(k) -> H(-k)`.
"""
function reciprocal_constraints_trs(
    Mm::Array{Int,4},
    h_orbit::HoppingOrbit{D}
) where {D}
    Z_trs = zeros(Int, size(Mm))
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

Time reversal symmetry action on the Hamiltonian. It is given by the association δ -> -δ and 
the complex conjugation in the free-parameter part:
``tⱼ -> tⱼ* ⇒ [tⱼᴿ, itⱼᴵ] -> [tⱼᴿ, -itⱼᴵ] ⇒ H(k) -> H*(k)``.
"""
function representation_constraint_trs(
    Mm::AbstractArray{<:Number,4},
    h_orbit::HoppingOrbit{D},
) where {D}
    Q_trs = zeros(Int, size(Mm))
    op = SymOperation([-I(D) zeros(D)])
    P = _permute_symmetry_related_hoppings_under_symmetry_operation(h_orbit, op)
    for l in axes(P, 2), j in axes(Mm, 2), s in axes(Mm, 3), t in axes(Mm, 4)
        Q_trs[l, j, s, t] = sum(P[i, l] * Mm[i, j, s, t] for i in axes(P, 1))
    end
    return [Q_trs -Q_trs]
end
