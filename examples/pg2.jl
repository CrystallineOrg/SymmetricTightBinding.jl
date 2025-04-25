using Pkg
Pkg.activate(@__DIR__)

using Crystalline, TETB

sgnum, D = 2, 1
brs = calc_bandreps(2, Val(D))
coefs = zeros(Int, length(brs))
coefs[[1, 3]] .= 1
cbr = CompositeBandRep(coefs, brs)

Rs = [[0,], [1,]]

tbs = tb_hamiltonian(cbr, Rs, false)