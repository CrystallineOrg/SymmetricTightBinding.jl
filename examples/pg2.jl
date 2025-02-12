using Pkg
Pkg.activate(@__DIR__)

using Crystalline, TETB

##- Compute the necessary things for obtaining the hoppings

pg_num, D = 2, 1
brs = calc_bandreps(2, Val(D))
coefs = zeros(Int, length(brs))
coefs[[1, 3]] .= 1
cbr = CompositeBandRep(coefs, brs)

Rs = [[0,]]

tbs = tb_hamiltonian(cbr, Rs)