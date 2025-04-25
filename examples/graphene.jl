using Pkg
Pkg.activate(@__DIR__)
using Crystalline, TETB

pg_num, D = 17, 2
brs = calc_bandreps(pg_num, Val(D))
coefs = zeros(Int, length(brs))
coefs[5] = 1
cbr = CompositeBandRep(coefs, brs)

tb_model = TETB.tb_hamiltonian(cbr, [[0, 0]], timereversal=true)