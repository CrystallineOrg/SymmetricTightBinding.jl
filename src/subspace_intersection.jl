# Intersection of subspaces. Conventionally, this could be computed by the Zassenhaus
# algorithm; but this is quite slow compared to using more traditional SVD-based approaches.
# Below, we exploit that the subspaces we want to intersect all come with orthonormalized
# bases, which allows some simplifications

"""
    subspace_intersection(U, W; atol) --> Matrix

Return an orthonormal basis for span(U) ∩ span(W), for matrices `U` and `W` with orthonormal
columns (as returned, e.g., by `nullspace`).

The implementation uses the SVD and the following facts:
- A vector `w` lies in span(U) iff its component orthogonal to span(U), `(I - UU†)w`,
  vanishes. Here, we use that `U` is orthonormal.
- The intersection span(U) ∩ span(W) can therefore be constructed by identifying all linear
  combinations of the columns of `W` that vanish under `(I - UU†)`, i.e., by forming
  `(I - UU†)W` and identifying its nullspace, `Q`. The intersection is then `WQ`.
- We identify the nullspace `Q` from the right singular vectors of `(I - UU†)W`, including
  only those whose singular values are below `atol`.
"""
function subspace_intersection(
    U::AbstractMatrix,
    W::AbstractMatrix;
    atol::Real # required
    # TODO: rtol; maybe some day
)
    # check input
    size(U,1) == size(W,1) || throw(DimensionMismatch("incompatible number of rows"))
    
    # early return for empty
    if (isempty(U) || isempty(W)) 
        T = float(promote_type(eltype(U), eltype(W)))
        return Matrix{T}(undef, size(U, 1), 0)
    end

    # project `W` onto orthogonal complement of `U`
    tmp = U * (U' * W) # = UU†W
    tmp .= W .- tmp    # = W - UU†W

    # compute `Q = nullspace(tmp; atol)`, but reusing `tmp` for `svd!`
    F = svd!(tmp) # thin SVD; size(W,2) ≤ size(W,1) by orthonormality (NB: `F.S` is sorted)
    indstart = count(>(atol), F.S) + 1 # get columns of `F.V` with singular values `<(atol)`
    Q = (@view F.Vt[indstart:end,:])'  # but take from `Vt` cf. strided memory layout

    # form the span by multiplying with `W` again; result is orthonormal, since W & Q are
    return W * Q
end