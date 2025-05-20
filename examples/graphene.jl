using Pkg
Pkg.activate(@__DIR__)
using Crystalline, TETB

sgnum, D = 17, 2
brs = calc_bandreps(sgnum, Val(D))
cbr = @composite brs[5]

tb_model = TETB.tb_hamiltonian(cbr; timereversal = true)
