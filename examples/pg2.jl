using Pkg
Pkg.activate(@__DIR__)

using Crystalline, TETB

sgnum, D = 2, 1
brs = calc_bandreps(2, Val(D))
cbr = @composite brs[1] + brs[3]

Rs = [[0], [1]]

tbs = tb_hamiltonian(cbr, Rs, false)
