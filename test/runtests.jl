using SymmetricTightBinding, Test

# The test files, in run order. Passing no arguments (the default, and what CI does) runs all
# of them; passing arguments runs only the matching subset, e.g.
#
#     julia --project -e 'using Pkg; Pkg.test(; test_args = ["show", "gradients"])'
#
# Arguments are matched as substrings of the names below, so e.g. "symmetry_analysis" selects
# both symmetry-analysis files. An argument matching nothing is an error, so that a typo'ed
# selector fails loudly instead of silently running an empty test suite.
#
# This exists because `symmetry_analysis` sweeps every EBR of every space group in 1D–3D and
# dominates the suite's runtime and memory use: verifying a limited change is much better
# served by running only the files that cover it. See also CLAUDE.md's "Running tests".
const TESTFILES = [
    "pg_tb_hamiltonian",        # plane groups
    "sg_tb_hamiltonian",        # space groups
    "site_representations",     # site representations
    "symmetry-breaking",        # symmetry breaking
    "berry",                    # berry curvature and chern numbers
    "dos",                      # density of states (Gilat–Raubenheimer)
    "spectrum",                 # spectrum evaluation
    "show",                     # show/display methods
    "nonhermitian",             # NONHERMITIAN models
    "gradients",                # hopping and momentum gradients
    "symmetry_analysis",        # each tb model is symmetry compatible w/ its constituent EBRs
    "symmetry_analysis_manual", # paired-down, manual version of above, individual cases
    "misc",                     # AbstractArray interface & assorted issue regressions
]

function select_testfiles(args)
    isempty(args) && return TESTFILES
    selected = String[]
    for arg in args
        matches = filter(testfile -> occursin(arg, testfile), TESTFILES)
        isempty(matches) && error(
            "test selector \"$arg\" matches none of the test files: $(join(TESTFILES, ", "))")
        append!(selected, matches)
    end
    return filter(in(selected), TESTFILES) # de-duplicate & restore canonical run order
end

for testfile in select_testfiles(ARGS)
    include(testfile * ".jl")
end
