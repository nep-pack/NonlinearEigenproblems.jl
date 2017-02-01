module Gallery
    using NEPCore
    using NEPTypes
    using MATLAB
    using PolynomialRoots
    
    export nep_gallery
    export nlevp_gallery_import
    # We have to explicitly specify functions that we want "overload"
    import NEPCore.compute_Mder
    import NEPCore.size
    
  """
  Returns a NEP object from a gallery of examples of nonlinear eigenvalue problems. name decides which NEP. \\
  Usage:\\
     nep=nep_gallery(name)\\
     or
     nep=nep_gallery(name,params)
  """
  function nep_gallery(name,params...)
      if (name == "dep0")
          # A delay eigenvalue problem
          n=100;
          srand(0) # reset the random seed
          A0=randn(n,n);
          A1=randn(n,n);
          tau=1;
          nep=DEP([A0,A1],[0,tau])
          return nep

      elseif (name == "dep0_sparse")
          # A delay eigenvalue problem with sparse matrices
          n=5;
          srand(0) # reset the random seed
          A0=sparse(randn(n,n));
          A1=sparse(randn(n,n));
          tau=1;
          nep=DEP([A0,A1],[0,tau])
          return nep

      elseif (name == "dep_double")
          # A delay eigenvalue problem with a double non-semisimple eigenvalue in λ=3πi
          # Examle from E. Jarlebring, Convergence factors of Newton methods for nonlinear eigenvalue problems, LAA, 2012
          n=3;

          denom = 8+5*pi;
          a1 = 2/5 *(65*pi + 32)/(denom);
          a2 = 9*pi^2*(13+5*pi)/(denom);
          a3 = 324/5 *pi^2*(5*pi+4)/(denom);
          b1 = (260*pi + 128 + 225*pi^2)/(10*denom);
          b2 = 45*pi^2/denom;
          b3 = 81*pi^2*(40*pi + 32 + 25*pi^2)/(10*denom);
          A0 = [ 0    1    0;  0    0    1;  -a3  -a2  -a1];
          A1 = [ 0    0    0;  0    0    0;  -b3  -b2  -b1];

          tau=1;
          nep=DEP([A0,A1],[0,tau])
          return nep

      elseif (name== "pep0")
          # A polynomial eigenvalue problem
          n=200; # mat size

          srand(0)
          A0=randn(n,n)
          A1=randn(n,n)
          A2=randn(n,n)
          A=[A0,A1,A2]
          nep=PEP(A)
          return nep

       elseif (name== "pep0_sparse_003")
          # A sparse polynomial eigenvalue problem
          n=200; # mat size

          srand(0)
          A0=sprandn(n,n,0.03)
          A1=sprandn(n,n,0.03)
          A2=sprandn(n,n, 0.03)
          A=[A0,A1,A2]
          nep=PEP(A)
          return nep

      end    
  end

    """
    nlevp_gallery(name)
Loads a NEP from the Berlin-Manchester collection of nonlinear
eigenvalue problems
"""
    function nlevp_gallery_import(name::String,nlevp_path::String="../../nlevp3")
        nep=NLEVP_NEP(name,nlevp_path)
        return nep
    end
  
    """
         NLEVP_NEP represents a NEP in the NLEVP-toolbox
    Example usage: nep=NLEVP_NEP("gun")
    """
    type NLEVP_NEP <: NEP
        n::Integer
        name::String
        Ai::Array
        function NLEVP_NEP(name,nlevp_path)
            if (~isfile(joinpath(nlevp_path,"nlevp.m")))
                error("nlevp.m not found. You need to install the Berlin-Manchester collection (http://www.maths.manchester.ac.uk/our-research/research-groups/numerical-analysis-and-scientific-computing/numerical-analysis/software/nlevp/) and specify a nlevp_path.")
            end
            @mput name nlevp_path
            @matlab begin
                addpath(nlevp_path)
                Ai,funs=nlevp(name)
            @matlab end
            @mget Ai # fetch and store the matrices
            this=new(size(Ai[1],1),name,Ai);
        end
    end
    
    function compute_Mder(nep::NLEVP_NEP,λ::Number,i::Integer=0)
        lambda=Complex{Float64}(λ)  # avoid type conversion problems
        #println("type",typeof(lambda))
        ## The following commented code is calling nlevp("eval",...)
        ## directly and does not work. We use functions instead
    #        nep_name::String=nep.name
    #        @mput lambda nep_name
    #        if (i==0)
    #            println(λ)
    #            @matlab begin
    #                ll=1+0.1i
    #                M=nlevp("eval",nep_name,lambda
    #            @matlab end
    #            @mget M
    #            return M
    #        else
    #            @matlab begin
    #                (M,Mp)=nlevp("eval",nep_name,lambda)
    #            @matlab end
    #            @mget Mp
    #            return Mp
    #        end
    #    return f,fp
        D=call_current_fun(lambda,i)
        f=D[i+1,:]
        M=zeros(nep.Ai[1]);
        for i=1:length(nep.Ai)
            M=M+nep.Ai[i]*f[i]
        end
        return M
    end
    
    # Return function values and derivatives of the current matlab session "funs"
    # stemming from a previous call to [Ai,funs]=nlevp(nepname).
    # The returned matrix containing derivatives has (maxder+1) rows 
    function call_current_fun(lambda,maxder::Integer=0)        
        l::Complex64=Complex64(lambda)  # avoid type problems
        @mput l maxder
        eval_string("C=cell(maxder+1,1); [C{:}]=funs(l); D=cell2mat(C);")
        @mget D
        return D
    end


    # size for NLEVP_NEPs
    function size(nep::NLEVP_NEP,dim=-1)
        if (dim==-1)
            return (nep.n,nep.n)
        else
            return nep.n
        end
    end
    
    
end
