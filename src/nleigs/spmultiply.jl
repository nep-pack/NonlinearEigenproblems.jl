# Faster way to multiply a sparse matrix by a sparse vector
# See https://github.com/JuliaLang/julia/issues/20447
function spmultiply(A::SparseMatrixCSC{TA,IA}, v::SparseVector{TV,IV}) where {TA,IA,TV,IV}
    A.n == v.n || throw(DimensionMismatch("$(A.n) != $(v.n)"))

    ind_out = promote_type(IA,IV)[]
    val_out = promote_type(TA,TV)[]

    for (vi,vv) in zip(v.nzind, v.nzval)
        lo_ptr = A.colptr[vi]
        hi_ptr = A.colptr[vi+1] - 1
        append!(ind_out, A.rowval[lo_ptr:hi_ptr])
        append!(val_out, A.nzval[lo_ptr:hi_ptr] * vv)
    end

    sparsevec(ind_out, val_out, A.m)
end

function spmultiply(A::SparseMatrixCSC{TA,IA}, v::AbstractVector{TV}) where {TA,IA,TV}
    spmultiply(A, sparse(v))
end

#test
if false
    n = 2_000_000
    p = 2e-5
    A = sprand(n, n, p)
    v = sprand(n, 8/n)

    @time A*v
    @time spmultiply(A, v) # should be faster
end
