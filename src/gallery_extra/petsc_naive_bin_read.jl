using SparseArrays;

"""
  function naive_petsc_read(filename, int_type=Int32, float_type=ComplexF64)

Loading petsc data into julia arrays and vector
The function does not support all Petsc class types.
It only supports sparse arrays and vectors.

File format documentation:
http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MatLoad.html
https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Vec/VecLoad.html
"""

function naive_petsc_read(filename, int_type=Int32, float_type=ComplexF64)
    open(filename) do io
        class_id = ntoh(read(io, int_type))

        MAT_FILE_CLASSID=1211216;
        VEC_FILE_CLASSID=1211214;

        if (class_id != MAT_FILE_CLASSID  && class_id != VEC_FILE_CLASSID)
            @show class_id
            error("Unsupported class_id. This function can only load sparse arrays and vectors: ")
        end

        if (class_id == MAT_FILE_CLASSID)
            # It's a sparse array

            # Get size
            rows = ntoh(read(io, int_type)) # From documentation: number of rows
            cols = ntoh(read(io, int_type)) # From documentation: number of columns

            # Number of nnz elements
            nnz = ntoh(read(io, int_type))  # From documentation: total number of nonzeros

            # Prep for reading
            row_ptr = Vector{int_type}(undef, rows+1)
            row_ptr[1] = 1

            ## Read data

            # row lengths
            # From documentation: number nonzeros in each row
            row_ptr[2:end] = ntoh.(read!(io, Vector{int_type}(undef, rows)))
            cumsum!(row_ptr, row_ptr)

            # columns
            # From documentation: column indices of all nonzeros
            colvals = ntoh.(read!(io, Vector{int_type}(undef, nnz)))
            colvals .+= int_type(1) # Start index is 0 in PETsc and 1 in Julia

            # From documentation: values of all nonzeros
            vals = ntoh.(read!(io, Vector{float_type}(undef, nnz)))

            mat = SparseMatrixCSC(cols, rows, row_ptr, colvals, vals)
            return SparseMatrixCSC(transpose(mat))
        else
            # It's a vec

            # From documentation: number of rows
            rows = ntoh(read(io, int_type))

            # From documentation: values of all entries
            vec=ntoh.(read!(io, Vector{float_type}(undef, rows)))
        end

    end

end
