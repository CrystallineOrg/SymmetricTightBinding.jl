using TETB, Test

@testset "TETB" begin
    include("pg_tb_hamiltonian.jl")    # plane groups
    include("sg_tb_hamiltonian.jl")    # space groups
    include("site_representations.jl") # site representations
end
