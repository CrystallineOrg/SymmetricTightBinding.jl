# implementation of the Zassenhaus algorithm.

using RowEchelon: rref

# my idea is to have U and W as the matrices {aᵢⱼ} and {bᵢⱼ} from the Zassenhaus algorithm
# U is a matrix of n×m and W is a matrix of k×m.

"""
    zassenhaus_intersection(
                            U::AbstractArray{<:Number},
                            W::AbstractArray{<:Number}) -> AbstractArray{<:Number}

Finds the intersection of two bases `U` and `W` using the Zassenhaus algorithm.
It assumes that the basis are given by columns.

# References: 
- https://en.wikipedia.org/wiki/Zassenhaus_algorithm
"""
function zassenhaus_intersection(
    U::AbstractArray{T},
    W::AbstractArray{T}
) where {T<:Number}
    U = transpose(U) # the algorithm assumes that the basis are given by rows
    W = transpose(W) # for me is more natural to give it by columns, so I just
    # transpose them   transpose them
    size(U, 2) == size(W, 2) || error("The matrices must have the same number of columns")
    # Combine the bases of U and W into a single matrix
    A = [U U; W zeros(size(W, 1), size(U, 2))]

    # Perform row reduction (Gaussian elimination)
    R = rref(A)

    inter_basis = Vector{Vector{T}}()
    for v in eachrow(R)
        if all(vᵢ -> abs(vᵢ) < ZASSENHAUS_ATOL_DEFAULT, @view v[1:size(U, 2)])
            append!(inter_basis, [v[size(U, 2)+1:end]])
        end
    end

    return inter_basis
end