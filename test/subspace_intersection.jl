using Test
using LinearAlgebra
using Random
using SymmetricTightBinding: subspace_intersection

# random `n×k` matrix with orthonormal columns
randorth(T, n, k) = Matrix(qr(randn(T, n, k)).Q)

# random orthonormal `U` (`n×k`) and `W` (`n×l`) whose spans share exactly the span of `C`
# (`n×r`), provided `k + l - r ≤ n`
function predetermined_overlap(T, n, k, l, r)
    C = randorth(T, n, r)
    U = Matrix(qr([C randn(T, n, k - r)]).Q)
    W = Matrix(qr([C randn(T, n, l - r)]).Q)
    return U, W, C
end

# whether the columns of `A` and `B` span the same subspace
samespan(A, B) = size(A, 2) == size(B, 2) && rank([A B]; atol = 1e-8) == size(A, 2)

@testset "subspace_intersection" begin
    Random.seed!(141)
    atol = 1e-5

    @testset "Predetermined overlap ($T, n=$n)" for T in (Float64, ComplexF64),
                                              (n, k, l, r) in ((12, 6, 5, 3), (60, 30, 35, 10))
        U, W, C = predetermined_overlap(T, n, k, l, r)
        Z = subspace_intersection(U, W; atol)
        @test size(Z) == (n, r)
        @test Z' * Z ≈ I        # orthonormal
        @test samespan(Z, C)    # spans exactly the predetermined overlap
    end

    # regression: a relative tolerance scaled by the largest singular value of `(I - UU†)W`
    # collapses when span(W) ⊆ span(U), since all singular values are then rounding noise
    @testset "Nested subspaces" begin
        U = randorth(Float64, 20, 8)
        @test samespan(subspace_intersection(U, U; atol), U)
        W = randorth(Float64, 8, 5)
        W = U * W # rotated 5-dimensional subspace of span(U)
        @test samespan(subspace_intersection(U, W; atol), W)
        @test samespan(subspace_intersection(W, U; atol), W) # symmetric in arguments
    end

    @testset "Generic subspaces do not intersect" begin
        U, W = randorth(Float64, 20, 8), randorth(Float64, 20, 12) # k + l ≤ n
        @test size(subspace_intersection(U, W; atol)) == (20, 0)
    end

    # `atol` bounds the sine of the angle between a direction in span(W) and span(U)
    @testset "Tolerance" begin
        e = Matrix{Float64}(I, 4, 4)
        U = e[:, 1:2]
        for (θ, dim) in ((1e-7, 1), (1e-3, 0))
            W = reshape(cos(θ) * e[:, 1] + sin(θ) * e[:, 3], :, 1)
            @test size(subspace_intersection(U, W; atol), 2) == dim
        end
    end

    @testset "Empty input and dimension mismatch" begin
        U = randorth(ComplexF64, 6, 3)
        Z = subspace_intersection(U, zeros(6, 0); atol)
        @test size(Z) == (6, 0) && eltype(Z) == ComplexF64
        @test size(subspace_intersection(zeros(6, 0), U; atol)) == (6, 0)
        R = randorth(Float64, 5, 2)
        @test_throws DimensionMismatch subspace_intersection(U, R; atol)
    end
end
