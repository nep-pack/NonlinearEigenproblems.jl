module Serialization
    using ZipFile

    export write_sparse_matrices
    export read_sparse_matrices

    # Writes the sparse matrices in text format inside a zip file
    function write_sparse_matrices{S<:AbstractString}(zipname::AbstractString, matrices::Dict{S,SparseMatrixCSC{Float64,Int64}})
        zipfile = ZipFile.Writer(zipname)

        for (name, M) in matrices
            f = ZipFile.addfile(zipfile, name, method = ZipFile.Deflate)
            write(f, "$(size(M, 1))\n$(size(M, 2))\n")
            nz = findnz(M)
            for n = 1:3
                foreach(x -> write(f, "$x\n"), nz[n])
            end
        end

        close(zipfile)
    end

    # Reads and returns the sparse matrices written by write_sparse_matrices
    function read_sparse_matrices(zipname::AbstractString)
        matrices = Dict{String,SparseMatrixCSC{Float64,Int64}}()

        zipfile = ZipFile.Reader(zipname)

        for f in zipfile.files
            data = readlines(f)
            m = parse(Int, data[1])
            n = parse(Int, data[2])
            c = Integer((length(data) - 2) / 3)
            I = map(x -> parse(Int, x), data[3:3+c-1])
            J = map(x -> parse(Int, x), data[3+c:3+2*c-1])
            V = map(x -> parse(Float64, x), data[3+2*c:3+3*c-1])
            matrices[f.name] = sparse(I, J, V, m, n)
        end

        close(zipfile)

        return matrices
    end
end
