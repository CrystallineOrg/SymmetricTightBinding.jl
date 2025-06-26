using SymmetricTightBinding, Test

include("pg_tb_hamiltonian.jl")    # plane groups
include("sg_tb_hamiltonian.jl")    # space groups
include("site_representations.jl") # site representations

#include("symmetry_analysis.jl")    # check that each tb model is symmetry compatible with
#                                    the constituents EBRs 