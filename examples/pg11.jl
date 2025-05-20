using Pkg
Pkg.activate(@__DIR__)

using Crystalline, TETB

D, sgnum = 2, 11
brs = calc_bandreps(sgnum, Val(D))
cbr = @omposite brs[14] # pick (1a|E) EBR
Rs = [[1, 0], [1, 1]]
tbs = TETB.tb_hamiltonian(cbr, Rs)
