module Serialization
    using SparseArrays

    export write_sparse_matrix
    export read_sparse_matrix

    # Writes the sparse matrix M in text format to a file with the given name
    function write_sparse_matrix(filename, M, T=Int64)
        open(filename, "w") do f
            write(f, "$(T(size(M, 1)))\n$(T(size(M, 2)))\n")
            nz = findnz(M)
            for n = 1:2
                foreach(x -> write(f, "$(T(x))\n"), nz[n])
            end
            foreach(x -> write(f, "$(Float64(x))\n"), nz[3])
        end
    end

    # Reads and returns a sparse matrix written by write_sparse_matrix
    function read_sparse_matrix(filename, T=Int64)
        data = open(filename) do f
            readlines(f)
        end
        m = parse(T, data[1])
        n = parse(T, data[2])
        c = Integer((length(data) - 2) / 3)
        I = map(x -> parse(T, x), data[3:3+c-1])
        J = map(x -> parse(T, x), data[3+c:3+2*c-1])
        V = map(x -> parse(Float64, x), data[3+2*c:3+3*c-1])
        return sparse(I, J, V, m, n)
    end
end
