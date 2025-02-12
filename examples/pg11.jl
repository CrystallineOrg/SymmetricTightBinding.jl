using Pkg
Pkg.activate(@__DIR__)

using Crystalline, TETB

D, sgnum = 2, 11
brs = calc_bandreps(sgnum, Val(D))
coefs = zeros(Int, length(brs));
coefs[end] = 1; # pick (1a|E) EBR
cbr = CompositeBandRep(coefs, brs)
Rs = [[1, 0], [1, 1]]
tbs = TETB.tb_hamiltonian(cbr, Rs)
