using Pkg
Pkg.activate(@__DIR__)

using Crystalline, SymmetricTightBinding

sgnum, D = 224, 3
brs = calc_bandreps(sgnum, Val(D))
cbr = @composite brs[13] + brs[19]
Rs = [[0, 0, 0]]

tbm = tb_hamiltonian(cbr, Rs)

ptbm = tbm(randn(length(tbm)))

println(SymmetryVector(cbr))
println(collect_compatible(ptbm))
