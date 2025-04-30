using Pkg
Pkg.activate(@__DIR__)

using Crystalline, TETB

# This example shows a non-inversion symmetric case with 2D complex site-symmetry representations
pgnum, D = 13, 2
timereversal = false
brs = calc_bandreps(pgnum, Val(D); timereversal)

# The band representations are: (1b|A) and (1c|A), they will depend on the choice of timereversal
coefs = zeros(length(brs))

if timereversal
    br₁ = brs[1]
    br₂ = brs[3]
    coefs[[1, 3]] .= 1
else
    br₁ = brs[1]
    br₂ = brs[4]
    coefs[[1, 4]] .= 1
end

cbr = CompositeBandRep(coefs, brs)

# I will study here the term involving br₁ → br₂, which will be a non-diagonal term

Rs = [[0, 0]]

hops = obtain_symmetry_related_hoppings(Rs, br₁, br₂; timereversal)

tb_model = tb_hamiltonian(cbr, Rs; timereversal)