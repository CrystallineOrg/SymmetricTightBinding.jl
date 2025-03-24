using TETB, Test

@testset "TETB" begin
    # Plane groups
    include("pg_tb_hamiltonian.jl")

    # Space groups
    include("sg_tb_hamiltonian.jl")

    # Site representations
    include("site_representations.jl")

    # decomposition into EBR
    include("mpb_ebr_decomposition.jl")
end
