using TETB, Test

@testset "TETB" begin
    # decompositions into EBRs
    include("ebr_decomposition.jl")
    include("mpb_ebr_decomposition.jl")
end
