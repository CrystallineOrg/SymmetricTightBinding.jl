using SymmetricTightBinding, Test

# The test files, in run order. `Pkg.test` forwards its `test_args` here as `ARGS`: with no
# arguments every file runs, while each argument selects a file, and each "-"-prefixed
# argument excludes one. Exclusions alone run everything else; alongside inclusions, they are
# subtracted from those. A selector must name a file below exactly, with or without its ".jl"
# extension; anything else is an error.
#
#     Pkg.test(; test_args = ["show", "gradients"])   # just these two
#     Pkg.test(; test_args = ["-symmetry_analysis"])  # all but the slow sweep
#
# CLAUDE.md's "Running tests" gives per-file timings and the recommended local invocation.
const TESTFILES = [
    "pg_tb_hamiltonian.jl",        # plane groups
    "sg_tb_hamiltonian.jl",        # space groups
    "site_representations.jl",     # site representations
    "symmetry-breaking.jl",        # symmetry breaking
    "berry.jl",                    # berry curvature and chern numbers
    "dos.jl",                      # density of states (Gilat–Raubenheimer)
    "spectrum.jl",                 # spectrum evaluation
    "show.jl",                     # show/display methods
    "nonhermitian.jl",             # NONHERMITIAN models
    "gradients.jl",                # hopping and momentum gradients
    "fitting.jl",                  # fitting (Optim extension) + TightBindingCache
    "symmetry_analysis.jl",        # ⚠️ every EBR of every SG in 1D-3D; minutes to hours
    "symmetry_analysis_manual.jl", # paired-down, manual version of above, individual cases
    "misc.jl",                     # AbstractArray interface & assorted issue regressions
]

"Resolve a selector to the test file it names, erroring if it names none."
function match_testfile(selector)
    testfile = endswith(selector, ".jl") ? selector : selector * ".jl"
    testfile ∈ TESTFILES || error(
        "unknown test selector \"$selector\"; must name one of: $(join(TESTFILES, ", "))")
    return testfile
end

function select_testfiles(args)
    isempty(args) && return TESTFILES
    included, excluded = String[], String[]
    for arg in args
        if startswith(arg, '-')
            push!(excluded, match_testfile(arg[nextind(arg, 1):end]))
        else
            push!(included, match_testfile(arg))
        end
    end
    isempty(included) && (included = TESTFILES) # exclusions alone ⇒ run everything else
    # filtering over `TESTFILES` de-duplicates and restores the canonical run order
    return filter(testfile -> testfile ∈ included && testfile ∉ excluded, TESTFILES)
end

for testfile in select_testfiles(ARGS)
    include(testfile)
end
