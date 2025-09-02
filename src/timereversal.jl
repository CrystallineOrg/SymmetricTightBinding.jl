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
# `D(𝒯)H(k)D(𝒯)† = H(𝒯k) -> Γ(𝒯)H*(k)Γ(𝒯)† = H(-k)`, where `D` is the whole 
# operator and `Γ` is only the unitary part, so `D(𝒯) = Γ(𝒯)𝒯`
# If the site-symmetry irrep is real, `Γ(𝒯) = I -> H*(k) = H(-k)`.

# we have made the following decomposition: `t -> [real(t); imag(t)]`.
# how does this affect the Hamiltonian representation?

# we have that `Hᵢⱼ(k) = vᵀ(k) Mᵢⱼ t`. Then, we need to make the following change:
# `Hᵢⱼ(k) = vᵀ(k) [Mᵢⱼ Mᵢⱼ] [real(t); imag(t)]`

# Then, applying TRS will be:

# 1. `H(-k) = vᵀ(-k) [Mᵢⱼ Mᵢⱼ] [real(t); imag(t)] = (Pv)ᵀ(k) [Mᵢⱼ Mᵢⱼ] 
# [real(t); imag(t)]`

# 2. `H*(k) = (v*)ᵀ(k) [Mᵢⱼ Mᵢⱼ] [real(t); -imag(t) -imag(im*t)] = (Pv)ᵀ(k) 
# [Mᵢⱼ -Mᵢⱼ] [real(t); imag(t)]`

# Imposing the condition `H(-k) = H*(k)` we get:
# `(Pv)ᵀ(k) [Mᵢⱼ Mᵢⱼ] [real(t); imag(t)] = (Pv)ᵀ(k) [Mᵢⱼ -Mᵢⱼ] [real(t); imag(t)]`
# which can be rewritten as:
# `(Pv)ᵀ(k) [0 2Mᵢⱼ] [real(t); imag(t)] = 0`

# NOTE: This way of casting the problem allows us to not use the "extended" 𝐯-vector with 
# the reversed hopping in non-diagonal blocks, which could potentially enforce extra 
# symmetries in the system.

"""
    obtain_basis_free_parameters_TRS(
        h_orbit::HoppingOrbit{D}, 
        brₐ::NewBandRep{D}, 
        brᵦ::NewBandRep{D}, 
        orderingₐ::OrbitalOrdering{D} = OrbitalOrdering(brₐ),
        orderingᵦ::OrbitalOrdering{D} = OrbitalOrdering(brᵦ),
        Mm::AbstractArray{4, Int} = construct_M_matrix(h_orbit, brₐ, brᵦ, orderingₐ, orderingᵦ)
        )                             --> Tuple{Array{Int,4}, Vector{Vector{ComplexF64}}}}

Obtain the basis of free parameters for the hopping terms between `brₐ` and `brᵦ` associated
with the hopping orbit `h_orbit` under time-reversal symmetry.

Real and imaginary parts of the basis vectors are differentiated explicitly: internally,
we consider only variables.
"""
function obtain_basis_free_parameters_TRS(
    h_orbit::HoppingOrbit{D},
    brₐ::NewBandRep{D},
    brᵦ::NewBandRep{D},
    orderingₐ::OrbitalOrdering{D} = OrbitalOrdering(brₐ),
    orderingᵦ::OrbitalOrdering{D} = OrbitalOrdering(brᵦ),
    Mm::AbstractArray{Int, 4} = construct_M_matrix(h_orbit, brₐ, brᵦ, orderingₐ, orderingᵦ),
) where {D}
    # To add time-reversal constraints, we duplicate the `M` matrix so that we can transfer
    # the conjugation from `t` to `M`. We exploit the following splitting:
    #    `Hᵢⱼ = vᵀ Mᵢⱼ t = vᵀ Mᵢⱼ (tᴿ + itᴵ) = vᵀ [Mᵢⱼ Mᵢⱼ] [tᴿ; itᴵ]`

    # Step 1: compute the Z tensor, encoding time-reversal constraints on H for the k-space
    # part. This is done by `Hᵢⱼ(-k) = vᵀ(-k) [Mᵢⱼ Mᵢⱼ] [tᴿ; tᴵ]`
    # `= (Pv)ᵀ(k) [Mᵢⱼ Mᵢⱼ] [tᴿ; tᴵ]`
    Z = [Mm Mm]

    # Step 2: compute the Q tensor, encoding time-reversal constraints on H for the free-
    # parameter part. This is done by `H*(k) = (v*)ᵀ(k) [Mᵢⱼ Mᵢⱼ] [tᴿ; -itᴵ]`
    # `= (Pv)ᵀ(k) [Mᵢⱼ -Mᵢⱼ] [tᴿ; tᴵ]`
    Q = [Mm -Mm]

    # Step 3: build an constraint matrix acting on the doubled hopping coefficient vector
    # `[tᴿ; itᴵ]` associated with `h_orbit`; each row is a constraint
    constraints = _aggregate_constraints(Q, Z)
    tₐᵦ_basis_matrix = nullspace(constraints; atol = NULLSPACE_ATOL_DEFAULT)

    return tₐᵦ_basis_matrix
end
