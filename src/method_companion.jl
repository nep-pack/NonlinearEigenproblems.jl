#Methods to transform a PEP to a companion linearized form and solve the corresponding linearized pencil 
export companion
export polyeig #Wrapper around the solver for a linearized PEP pencil

"""
    Return the most commonly used companion linearization(as in "Non-linear eigenvalue problems, a challenge for modern eigenvalue methods", by Mehrmann and Voss) of a PEP. For a k-th degree PEP with n-by-n coefficient matrices, this returns E and A, both kn-by-kn, the linearized problem is Ax = λEx
"""
    function companion(pep::PEP)

        n = size(pep,1);#Size of monomial coefficient matrices

        d = size(pep.A,1)-1;#Degree of pep

        T = eltype(pep.A[1]);#Deduce the type of elements in the PEP matrices
        
        

        ##### Check sparsity of the problem and allocate memory accordingly #####
        if (issparse(pep))
            E = spzeros(T,d*n,d*n);
            A = spzeros(T,n*d,n*d);
            
            #(d-1)n-by-(d-1)n matrix (Used to construct both E and A)
            Iblock = sparse(kron(eye(T,d-1),eye(T,n)));
        else
            E = zeros(T,d*n,d*n);
            A = zeros(T,n*d,n*d);
            
            Iblock = kron(eye(T,d-1),eye(T,n));
        end
        
        E[1:n,1:n] = pep.A[d+1];#Fill block (1,1)
        E[n+1:d*n,n+1:d*n] = Iblock;#Fill all blocks on the diagonal with eye(n)

        #####Construct A #####
        
        #First row block of A
        for i=1:d
           A[1:n,(i-1)*n+1:i*n] = pep.A[d-i+1];
        end
        #Lower part of A
        A[n+1:d*n,1:(d-1)*n] = T(-1.0)*Iblock 

        return E,-A


    end

    #############################################################################
    #Solve the linearized companion of a PEP
    function polyeig(pep::PEP,eigsolvertype::DataType=DefaultEigSolver)

        #Linearize to Ax = λEx
        E,A = companion(pep);

        solver::EigSolver = eigsolvertype(A,E);
        D,V = eig_solve(solver,target=1.0,nev=size(A)[1]);

        return D,V
    end
