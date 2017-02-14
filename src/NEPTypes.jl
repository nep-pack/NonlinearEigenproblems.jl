module NEPTypes
    # Specializalized NEPs 
    export DEP
    export PEP
    export REP
    export SPMF_NEP    

    export interpolate
    export interpolate_cheb
    
    using NEPCore

    # We overload these
    import NEPCore.compute_Mder
    import NEPCore.compute_Mlincomb
    import NEPCore.compute_MM
    import NEPCore.compute_resnorm
    import NEPCore.compute_rf
    
    import Base.size
    import Base.issparse

    export compute_Mder
    export compute_Mlincomb
    export compute_MM
    export compute_resnorm
    export compute_rf    
    export size
    export companion



    
    ###########################################################
    # Delay eigenvalue problems - DEP
    #

    """
### Delay eigenvalue problem
  A DEP is defined by the sum the sum  ``-λI + Σ_i A_i exp(-tau_i λ)``\\
  where all of the matrices are of size n times n\\
  Constructor: DEP(AA,tauv) where AA is an array of the\\
  matrices A_i, and tauv is a vector of the values tau_i
"""
    type DEP <: NEP
        n::Integer
        A::Array     # An array of matrices (full or sparse matrices)
        tauv::Array{Float64,1} # the delays
        function DEP(AA,tauv=[0,1.0])
            n=size(AA[1],1)
            this=new(n,AA,tauv);
            return this;
        end
    end

    """
    compute_Mder(nep::DEP,λ::Number,i::Integer=0)
 Compute the ith derivative of a DEP
"""
     # TODO: this function compute only the first 2 derivatives. Extend.
    function compute_Mder(nep::DEP,λ::Number,i::Integer=0)
        local M,I;
        if issparse(nep)
            M=spzeros(nep.n,nep.n)
            I=speye(nep.n,nep.n)
        else
            M=zeros(nep.n,nep.n)
            I=eye(nep.n,nep.n)        
        end
        if i==0; M=-λ*I;  end
        if i==1; M=-I; end
        for j=1:size(nep.A,1)
            M+=nep.A[j]*(exp(-nep.tauv[j]*λ)*(-nep.tauv[j])^i)
        end
        return M
    end



    """
    compute_MM(nep::DEP,S,V)
 Computes the sum ``Σ_i M_i V f_i(S)`` for a DEP
"""
    function compute_MM(nep::DEP,S,V)
        Z=-V*S;
        for j=1:size(nep.A,1)
            Z+=nep.A[j]*V*expm(full(-nep.tauv[j]*S))
        end
        return Z
    end

    ###########################################################
    # Polynomial eigenvalue problem - PEP
    #

    """
### Polynomial eigenvalue problem
  A PEP is defined by the sum the sum ``Σ_i A_i λ^i``,\\
  where i = 0,1,2,..., and  all of the matrices are of size n times n\\
  Constructor: PEP(AA) where AA is an array of the matrices A_i
"""

    type PEP <: NEP
        n::Integer
        A::Array   # Monomial coefficients of PEP
        function PEP(AA)
            n=size(AA[1],1)
            return new(n,AA)
        end
    end

    """
    compute_MM(nep::DEP,S,V)
 Computes the sum ``Σ_i M_i V f_i(S)`` for a DEP
"""
    function compute_MM(nep::PEP,S,V)
        if (issparse(nep))
            Z=spzeros(size(V,1),size(V,2))
            Si=speye(size(S,1))
        else
            Z=zeros(size(V))
            Si=eye(size(S,1))
        end
        for i=1:size(nep.A,1)
            Z+=nep.A[i]*V*Si;
            Si=Si*S;
        end
        return Z
    end

    """
    compute_Mder(nep::PEP,λ::Number,i::Integer=0)
 Compute the ith derivative of a PEP
"""
    function compute_Mder(nep::PEP,λ::Number,i::Integer=0)
        if (issparse(nep))
            Z=spzeros(size(nep,1),size(nep,1));
        else
            Z=zeros(size(nep,1),size(nep,1));
        end
        for j=(i+1):size(nep.A,1)
            # Derivatives of monimials
            Z+= nep.A[j]*(λ^(j-i-1)*factorial(j-1)/factorial(j-i-1))
        end
        return Z
    end



"""
    interpolate([T::DataType=Complex64,] nep::NEP, intpoints::Array)
 Interpolates a NEP in the points intpoints and returns a PEP.\\
 `T` is the DataType in which the PEP should be defined.
"""
    function interpolate(T::DataType, nep::NEP, intpoints::Array)

        n = size(nep, 1)
        d = length(intpoints)
        
        V = Array{T}(d,d) #Vandermonde matrix
        pwr = ones(d,1)
        for i = 1:d
            V[:,i] = pwr
            pwr = pwr.*intpoints
        end

        if (issparse(nep)) #If Sparse, do elementwise interpolation
            b = Array{SparseMatrixCSC{T},1}(d)
            AA = Array{SparseMatrixCSC{T},1}(d)
            V = factorize(V) # Will be used multiple times, factorize

            for i=1:d
                b[i] = compute_Mder(nep, intpoints[i])
            end
            
            # OBS: The following lines and hence the  following method assumes that Sparsity-structure is the same!
            nnz_AA = nnz(b[1])
            for i=1:d
                AA[i] = spones(b[1])
            end
            
            f = zeros(d,1)
            for i = 1:nnz_AA
                for j = 1:d
                    f[j] = b[j].nzval[i]
                end
                a = \(V,f)
                for j = 1:d
                    AA[j].nzval[i] = a[j]
                end
            end

        else # If dense, use Vandermonde
            b = Array{T}(n*d,n)
            AA = Array{Array{T,2}}(d)
            (L, U, p) = lu(V)

            I = speye(n,n)
            LL = kron(L,I)
            UU = kron(U,I)

            for i = 1:d
                b[(1:n)+(i-1)*n,:] =  compute_Mder(nep,intpoints[p[i]])
            end

            A = \(UU, \(LL,b))

            for i = 1:d
                AA[i] = A[(1:n)+(i-1)*n,:]
            end
        end
        
        return PEP(AA)
    end


    interpolate(nep::NEP, intpoints::Array) = interpolate(Complex128, nep, intpoints)


    """
     interpolate_cheb(nep::NEP,a::Real,b::Real)
  Interpolation in an interval using Chebyshev distribution. Returns a PEP.
  Following Effenberger, Cedric, and Daniel Kressner. "Chebyshev interpolation for nonlinear eigenvalue problems." BIT Numerical Mathematics 52.4 (2012): 933-951.
"""
    function interpolate_cheb(nep::NEP,a::Real,b::Real)
        # Not yet implemented
        # Note: PEP should probably be separated into Mono_PEP and
        # Cheb_PEP depending which inherit from PEP.
    end


"""
    Return the most commonly used companion linearization(as in "Non-linear eigenvalue problems, a challenge for modern eigenvalue methods", by Mehrmann and Voss) of a PEP. For a k-th degree PEP with n-by-n coefficient matrices, this returns E and A, both kn-by-kn, the linearized problem is Ax = λEx
"""
    function companion(pep::PEP)

        n = size(pep,1);#Size of monomial coefficient matrices

        d = size(pep.A,1)-1;#Degree of pep

        T = eltype(pep.A[1]);#Deduce the tyle of elements in the PEP matrices
        
        

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

    ###########################################################
    # Rational eigenvalue problem - REP
        
    """
### Rational eigenvalue problem
  A Rep is defined by the sum the sum ``Σ_i A_i s_i(λ)/q_i(λ)``,\\
  where i = 0,1,2,..., all of the matrices are of size n times n
  and s_i and q_i are polynomials\\
  Constructor: REP(AA) where AA is an array of the matrices A_i
"""

    type REP <: NEP
        n::Integer
        A::Array   # Monomial coefficients of REP
        si::Array  # numerator polynomials
        qi::Array  # demonimator polynomials
        # Initiate with order zero numerators and order one denominators
        # with poles given by poles[]
        function REP(AA,poles::Array)

            n=size(AA[1],1)
            # numerators
            si=Array{Array{Number,1},1}(length(poles))
            for i =1:size(poles,1)
                si[i]=[1];
            end
            # denominators
            qi=Array{Array{Number,1}}(length(poles))            
            for i =1:size(poles,1)
                if poles[i]!=0
                    qi[i]=[1,-poles[i]];
                else
                    qi[i]=[1];                    
                end
            end
            return new(n,AA,si,qi)
        end
    end

    function compute_MM(nep::REP,S,V)
        if (issparse(nep))
            Z=spzeros(size(V,1),size(V,2))
            Si=speye(size(S,1))
        else
            Z=zeros(size(V))
            Si=eye(size(S,1))
        end
        # Sum all the elements            
        for i=1:size(nep.A,1)
            # compute numerator
            Snum=copy(Si);
            Spowj=copy(Si);
            for j=1:length(nep.si[i])
                Snum+=Spowj*nep.si[i][j]
                Spowj=Spowj*S;
            end

            # compute denominator
            Sden=copy(Si);
            Spowj=copy(Si);
            for j=1:length(nep.qi[i])
                Sden+=Spowj*nep.qi[i][j]
                Spowj=Spowj*S;
            end
            
            # Sum it up 
            Z+=nep.A[i]*V*(Sden\Snum)
        end
        return Z
    end
    function compute_Mder(rep::REP,λ::Number,i::Integer=0)
        if (i!=0) # Todo
            error("Higher order derivatives of REP's not implemented")
        end
        S=speye(rep.n)*λ 
        V=speye(rep.n);
        return compute_MM(rep,S,V)  # This call can be slow
    end


    ####################################################### 
    ### Sum of products matrices and functions
    type SPMF_NEP <: NEP
         n::Integer
         A::Array   # Array of Array of matrices
         fi::Array  # Array of functions

         # Sparse zero matrix to be used for sparse matrix creation
         Zero::SparseMatrixCSC  
         function SPMF_NEP(AA,fii::Array)
             n=size(AA[1],1);
             if (issparse(AA[1]))
                 Zero=spones(AA[1]);                 
                 for i=2:length(AA)
                     Zero=Zero+spones(AA[i]);
                 end
                 Zero=(Zero*1im)*0
             else
                 Zero=zeros(n,n);
             end
             this=new(n,AA,fii,Zero);
             return this
         end       
    end
    function compute_MM(nep::SPMF_NEP,S,V)
        if (issparse(V))
            if (size(V)==size(nep))
                Z=copy(nep.Zero)
            else
                Z=spzeros(size(V,1),size(V,2))
            end
        else
            Z=zeros(size(V))
        end
        # Sum together all the terms in the SPMF:
        for i=1:length(nep.A)
            # Compute Fi=f_i(S)

            if (isdiag(S)) # optimize if S is diagonal
                Sd=diag(S);
                if (norm(Sd-Sd[1])==0) # Optimize further if S is a
                                       # multiple of identity
                    Fid=nep.fi[i](reshape([Sd[1]],1,1))[1]*ones(size(Sd,1))
                else  # Diagonal but not constant
                    Fid=zeros(Complex128,size(S,1))
                    for j=1:size(S,1)
                        Fid[j]=nep.fi[i](reshape([Sd[j]],1,1))[1]
                    end                    
                end
                Fi=spdiagm(Fid);
            else  # Otherwise just compute the matrix function operation
                Fi=nep.fi[i](S)
            end
            # Sum it up
            Z=Z+nep.A[i]*(V*Fi);
        end
        return Z
    end
    function compute_Mder(nep::SPMF_NEP,λ::Number,i::Integer=0)
        # This is typically slow for i>1 (can be optimized)
        return compute_Mder_from_MM(nep,λ,i) 
    end

  
         
   #######################################################
   ### Functions in common for many NEPs in NEPTypes

   #
"""
    size(nep::NEP,dim=-1)
 Overloads the size functions for NEPs storing size in nep.n
"""
    function size(nep::Union{DEP,PEP,REP,SPMF_NEP},dim=-1)
        if (dim==-1)
            return (nep.n,nep.n)
        else
            return nep.n
        end
    end

"""
    issparse(nep)
Returns true/false if the NEP is sparse (if compute_Mder() returns sparse)
"""
    function issparse(nep::Union{DEP,PEP,REP,SPMF_NEP})
        return issparse(nep.A[1])   
    end
        
end

