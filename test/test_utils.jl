# Shared helpers for the test suite. Included unconditionally by `runtests.jl`, ahead of the
# selected test files; individual test files that use these helpers also include this file
# themselves (guarded), so that they remain runnable on their own.

using Test
using DeepDiffs: deepdiff

"""
    test_show(expected::AbstractString, observed::AbstractString)

Compare two strings, printing a colored diff of the two if they differ.

Used for `show`-regression tests, which compare complete output rather than substrings.
Adapted from Crystalline.jl's `test/show.jl`.
"""
function test_show(expected::AbstractString, observed::AbstractString)
    if expected == observed
        @test true
    else
        old = Base.have_color
        @eval Base have_color = true
        try
            println(deepdiff(expected, observed))
        finally
            @eval Base have_color = $old
        end
        @test :expected == :observed
    end
end

"""
    test_tp_show(v, expected::AbstractString)

Compare the `MIME"text/plain"` representation of `v` against `expected`, via [`test_show`](@ref).
"""
test_tp_show(v, expected::AbstractString) = test_show(repr(MIME"text/plain"(), v), expected)
