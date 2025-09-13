
# TRS can be understood as a spatial symmetry when acting on the Hamiltonian:
#   `D(𝒯)H(k)D(𝒯)† = H(𝒯k) -> Γ(𝒯)H*(k)Γ(𝒯)† = H(-k)`, where `D` is the whole 
# operator and `Γ` is only the unitary part, so `D(𝒯) = Γ(𝒯)𝒯`
# If the basis of the representation is real, `Γ(𝒯) = I -> H*(k) = H(-k)`. We will
# assume this is the case from now on.

# NB: The "realification" of the representation is performed in Crystalline.jl when 
# time-reversal symmetry is present.

# We have that `Hᵢⱼ(k) = vᵀ(k) Mᵢⱼ t`. Remember that we split `t` into real and imaginary
# parts: `t = [real(t); i imag(t)]`, and `Mm` acts as `[Mm Mm]`. Then, we can rewrite
# the Hamiltonian as:
#   `Hᵢⱼ(k) = vᵀ(k) [Mᵢⱼ Mᵢⱼ] [real(t); i*imag(t)]`

# Then, applying TRS will be:

# 1. `H(-k) = vᵀ(-k) [Mᵢⱼ Mᵢⱼ] [real(t); i*imag(t)]`

# 2. `H*(k) = (v*)ᵀ(k) [Mᵢⱼ Mᵢⱼ] [real(t); -i*imag(t)] = (v*)ᵀ(k) [Mᵢⱼ -Mᵢⱼ] [real(t); i*imag(t)]
# = `vᵀ(-k) [Mᵢⱼ -Mᵢⱼ] [real(t); i*imag(t)]`

# where we make use of the fact that `v*(k) = v(-k)`.

# Imposing the condition `H(-k) = H*(k)` we get:
#   `vᵀ(-k) [Mᵢⱼ Mᵢⱼ] [real(t); i*imag(t)] = vᵀ(-k) [Mᵢⱼ -Mᵢⱼ] [real(t); i*imag(t)]`
# which can be rewritten as:
#   `vᵀ(-k) [0 2Mᵢⱼ] [real(t); i*imag(t)] = 0`

# This way of casting the problem is very convenient (and doesn't require any modification of 
#  the `v` vectors to include "reversed" hoppings, which is not in general a necessary feature).
# We just need to ensure that we intersect the `[real(t); i*imag(t)]` basis with the null space
# of  `[0 2Mᵢⱼ]`.

"""
    obtain_basis_free_parameters_TRS(
        h_orbit::HoppingOrbit{D}, 
        brₐ::NewBandRep{D}, 
        brᵦ::NewBandRep{D}, 
        orderingₐ::OrbitalOrdering{D} = OrbitalOrdering(brₐ),
        orderingᵦ::OrbitalOrdering{D} = OrbitalOrdering(brᵦ),
        Mm::AbstractArray{4, Int} = construct_M_matrix(h_orbit, brₐ, brᵦ, orderingₐ, orderingᵦ)
    )                             --> Tuple{Array{Int,4}, Vector{Vector{ComplexF64}}}

Obtain the basis of free parameters for the hopping terms between `brₐ` and `brᵦ` associated
with the hopping orbit `h_orbit` under time-reversal symmetry.

Real and imaginary parts of the basis vectors are differentiated explicitly.
"""
function obtain_basis_free_parameters_TRS(
    h_orbit::HoppingOrbit{D},
    brₐ::NewBandRep{D},
    brᵦ::NewBandRep{D},
    orderingₐ::OrbitalOrdering{D} = OrbitalOrdering(brₐ),
    orderingᵦ::OrbitalOrdering{D} = OrbitalOrdering(brᵦ),
    Mm::AbstractArray{Int, 4} = construct_M_matrix(h_orbit, brₐ, brᵦ, orderingₐ, orderingᵦ),
) where {D}
    # NB: we want to keep `_aggregate_constraints` due to its efficiency in building the
    # constraint matrix. That's why we are going to define to artificial tensors `Z` and `Q`
    # whose subtraction will result in `[0 2Mᵢⱼ]`.

    # Step 1: compute the Z tensor, which will just be `[0 2Mᵢⱼ]`
    Z = [zeros(Int, size(Mm)) 2 * Mm]

    # Step 2: compute the Q tensor, which will be just a matrix of zeros 
    Q = zeros(Int, size(Z))

    # Step 3: make use of `_aggregate_constraints` to build the constraint matrix
    constraints = _aggregate_constraints(Q, Z)
    tₐᵦ_basis_matrix = nullspace(constraints; atol = NULLSPACE_ATOL_DEFAULT)

    return tₐᵦ_basis_matrix
end
