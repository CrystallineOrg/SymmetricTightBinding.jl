using Pkg
Pkg.activate(@__DIR__)

using Crystalline, SymmetricTightBinding

sgnum = 10
brs = calc_bandreps(sgnum, Val(2))
cbr = @composite brs[1] + brs[end]

Rs = [[0, 0]]

## -------------------------------------------------------------------------------------- ##

sym_hops = obtain_symmetry_related_hoppings(Rs, cbr.brs[1], cbr.brs[5])

## -------------------------------------------------------------------------------------- ##

tb = tb_hamiltonian(cbr, Rs)
