# This example shows a particular plane group (P3 #13) which doesn't have inversion and 
# do contain a complex representation. Good example for testing TRS.

using Pkg
Pkg.activate(@__DIR__)
using Crystalline, TETB

sgnum, D = 13, 2
timereversal = false
brs = calc_bandreps(pgnum, Val(D); timereversal)

# we want to pick EBRs (1c|A) and (1b|A), which position will depend on timereversal
if timereversal
    br₁ = brs[1]
    br₂ = brs[3]
    cbr = @composite brs[1] + brs[3]
else
    br₁ = brs[1]
    br₂ = brs[4]
    cbr = @composite brs[1] + brs[4]
end

# chose a representative jump until where we will consider hopping terms.
Rs = [[0, 0]]

hops = obtain_symmetry_related_hoppings(Rs, br₁, br₂)

# this results in a vector v bigger than expected probably the sort of error.
