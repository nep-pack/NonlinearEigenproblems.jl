  """
  Module containing a gallery of examples of nonlinear eigenvalue problems.\\
  Look at the function 'nep_gallery()' for further instructions.
  """
module Gallery
    using NonlinearEigenproblems.Serialization
    using ..NEPCore
    using ..NEPTypes
    using Random
    using LinearAlgebra
    using SparseArrays
    using PolynomialRoots
    using SpecialFunctions

    export nep_gallery

    function __init__()
        # Add the search-path to the extra galleries - so that unloaded modules can be easily loaded from top level
        this_path = string(@__DIR__, "/gallery_extra")
        if !(this_path in LOAD_PATH)
            push!(LOAD_PATH, this_path)
        end
    end

    # Add the search-path to the extra galleries - temporarily so so needed files can be included
    push!(LOAD_PATH, string(@__DIR__, "/gallery_extra"))
    include("gallery_extra/distributed_example.jl")
    include("gallery_extra/periodic_dde.jl")
    include("gallery_extra/GalleryNLEVP_native.jl")

  """
     nep=nep_gallery(name)\\
     nep=nep_gallery(name,params)\\
**Returns a NEP object from a gallery of examples of nonlinear eigenvalue problems.**\\
    The parameter 'name' decides which NEP.\\
\\

\\
  **Supported "name" and 'params':**\\
\\
     `dep0`\\
Create a random delay eiganvalue problem with one delay tau = 1\\
      * one optional parameter determining the size (default = 5)\\
\\
     `dep0_sparse`\\
Create a random delay eiganvalue problem with sparse matrices and one delay tau = 1\\
      * two optional parameter determining the size (default = 5) and the fill (default = 0.25)\\
\\
      `dep0_tridiag`\\
Create a random delay eiganvalue problem with sparse tridiaognal matrices and one delay tau = 1\\
      * one optional parameter determining the size (default = 100)\\
\\
      `dep_symm_double`\\
Create delay eiganvalue problem with double eigenvalues and sparse symmetric matrices and one delay tau = 1\\
      * one optional parameter determining the size (default = 100)\\
       Examle from H. Voss and M. M. Betcke, Restarting iterative projection methods for Hermitian nonlinear eigenvalue problems with minmax property, Numer. Math., 2017\\
\\
     `dep_double`\\
Create problem with a double non-semisimple eigenvalue in λ=3πi\\
      Example from E. Jarlebring, Convergence factors of Newton methods for nonlinear eigenvalue problems, LAA, 2012\\
\\
     `dep1`\\
      A delay eigenvalue problem with one eigenvalue equal to one.\\
\\
     `pep0`\\
Create a random polynomial eigenvalue problem\\
      * one optional parameter determining the size (default = 200)\\
 \\
      `pep0_sym`\\
Create a random symmetric polynomial eigenvalue problem\\
      * one optional parameter determining the size (default = 200)\\
\\
     `pep0_sparse`\\
Create a random polynomial eigenvalue problem with sparse matrices
      * two optional parameter determining the size (default = 200) and the fill (default = 0.03)
\\
     `real_quadratic`\\
Create a quadratic problem with real eigenvalues\\
          Four smallest eigenvalues of the problem:\\
          -2051.741417993845\\
          -182.101627437811\\
          -39.344930222838\\
          -4.039879577113\\
\\
     `dep_distributed`\\
Creates the NEP associated with example in E. Jarlebring and W. Michiels and K. Meerbergen,  The infinite Arnoldi method and an application to time-delay systems with distributed delays,\\
     Delay Systems - Methods, Applications and New Trends, 2012\\
         Some correct eigenvalues:\\
         -0.400236388049641 + 0.970633098237807i\\
         -0.400236388049641 - 0.970633098237807i\\
          2.726146249832675 + 0.000000000000000i\\
         -1.955643591177653 + 3.364550574688863i\\
         -1.955643591177653 - 3.364550574688863i\\
          4.493937056300693 + 0.000000000000000i\\
         -1.631513006819252 + 4.555484848248613i\\
         -1.631513006819252 - 4.555484848248613i\\
         -1.677320660400946 + 7.496870451838560i\\
         -1.677320660400946 - 7.496870451838560i\\
\\
     `qdep0` \\
Quadratic delay eigenvalue problem in "The infinite Bi-Lanczos method for nonlinear eigenvalue problems",  Sarah W. Gaaf and Elias Jarlebring \\
\\

     `qdep1` \\
Quadratic delay eigenvalue problem in "A linear eigenvalue algorithm for the  nonlinear eigenvalue problem",      Elias Jarlebring, Wim Michiels, Karl Meerbergen \\
\\
     `qep_fixed_eig`\\
Create a quadratic eigenvalue problem with chosen eigenvalues \\
     * two optional parameters determining the size (default = 5)
       and a vector containing the eigenvalues (default = randn)       \\
\\
     `neuron0`\\
This problem stems from "Stability, bifurcation and multistability in a system of two coupled neurons with multiple time delays",
    L. P. Shayer and S. A. Campbell.
    SIAM J. Applied Mathematics, 2000. It is also a benchmark example in DDE-BIFTOOL\\
\\
     `beam`\\
The DEP modelling a beam with delayed stabilizing feedback described in "A rank-exploiting infinite Arnoldi algorithm for nonlinear eigenvalue problems", R. Van Beeumen, E. Jarlebring and W. Michiels, 2016. The A1-term has rank one.\\
     * one optional parameter which is the size of the matrix       \\
\\
     `sine` \\
The NEP formed by the sum of a polynomial and a sine-function in "A rank-exploiting infinite Arnoldi algorithm for nonlinear eigenvalue problems", R. Van Beeumen, E. Jarlebring and W. Michiels, 2016. The sine-term has rank one.\\
\\
     `nlevp_native_gun` \\
The benchmark problem from the NLEVP-collection called "gun", represented in the native NEP-PACK format. B.-S. Liao, Z. Bai, L.-Q. Lee, and K. Ko. Nonlinear Rayleigh-Ritz iterative method for solving large scale nonlinear eigenvalue pro blems.  Taiwan. Journal of Mathematics, 14(3):869–883, 2010\\
\\
     `nlevp_native_fiber` \\
The benchmark problem from the NLEVP-collection called "fiber", represented in the native NEP-PACK format. One of terms in this problem is approximated by interpolation, and may not always coincide with the benchmark. Kaufman, L. 2006. Eigenvalue problems in fiber optic design. SIAM J. Matrix Anal. Appl. 28, 1, 105–117.  and Huang, X., Bai, Z., and Su, Y. 2010. Nonlinear rank-one modification of the symmetric eigenvalue problem. J. Comput. Math. 28, 2, 218–234.\\
\\
   **See also the following galleries:**\\
      * GalleryNLEVP\\
      * GalleryWaveguide\\
  """
    nep_gallery(name::String,params...;kwargs...)=nep_gallery(NEP,name,params...;kwargs...)
    function nep_gallery(::Type{T},name::String,params...;kwargs...) where T<:NEP
        local n
        if (name == "dep0")
            return dep0(params...; kwargs...)

        elseif (name == "dep1")
            return dep1(params...; kwargs...)

        elseif (name == "dep0_sparse")
            return dep0_sparse(params...; kwargs...)

        elseif (name == "dep0_tridiag")
            return dep0_tridiag(params...; kwargs...)

        elseif (name == "dep_symm_double")
            return dep_symm_double(params...; kwargs...)

        elseif (name == "dep_double")
            return dep_double(params...; kwargs...)

        elseif (name == "pep0")
            return pep0(params...; kwargs...)

        elseif (name == "pep0_sym")
            return pep0_sym(params...; kwargs...)

        elseif (name == "pep0_sparse")
            return pep0_sparse(params...; kwargs...)

        elseif (name == "real_quadratic")
            return real_quadratic(params...; kwargs...)

        elseif (name== "dep_distributed")
            return gallery_dep_distributed(params...; kwargs...)

        elseif (name == "qdep0")
            return qdep0(params...; kwargs...)

        elseif (name == "qdep1")
            return qdep1(params...; kwargs...)

        elseif (name == "qep_fixed_eig")
            return qep_fixed_eig(params...; kwargs...)

        elseif (name == "periodicdde")
            return periodic_dde_gallery(PeriodicDDE_NEP, params...; kwargs...)

        elseif (name == "neuron0")
            return neuron0(params...; kwargs...)

        elseif (name == "nlevp_native_gun")
            return nlevp_native_gun(params...; kwargs...)

        elseif (name == "nlevp_native_fiber")
            return nlevp_native_fiber(params...; kwargs...)

      elseif (name == "beam")
          n::Int=100
          if (length(params)>0)
             n=params[1]
          end

          h=1/n;
          ee = ones(n);
          A0 = spdiagm(-1 => ee[1:n-1], 0 => -2*ee, 1 => ee[1:n-1]);
          A0[end,end]=1/h;
          A0[end,end-1]=-1/h;
          A1=sparse([n],[n],[1.0]); # A1=en*en'
          tau=1.0;
          return DEP([A0,A1],[0,tau]);

       elseif (name == "sine")
          data_dir=joinpath(dirname(@__FILE__()), "gallery_extra",   "converted_sine")


          A0=read_sparse_matrix(joinpath(data_dir,"sine_A0.txt"));
          A1=read_sparse_matrix(joinpath(data_dir,"sine_A1.txt"));
          A2=read_sparse_matrix(joinpath(data_dir,"sine_A2.txt"));
          V=Matrix(read_sparse_matrix(joinpath(data_dir,"sine_V.txt")));
          Q=Matrix(read_sparse_matrix(joinpath(data_dir,"sine_Q.txt")));

          n=size(A0,1);
          Z=spzeros(n,n);
          pep=PEP([A0,A1,Z,Z,A2]);
          # Matrix  sine function. Note that the Term is rank two which is not exploited here
          sin_nep=SPMF_NEP([V*Q'], [S-> sin(S)]);

          nep=SPMFSumNEP(pep,sin_nep) # Note: nep has a low-rank term
          return nep;
      else

          error("The name $name is not supported in NEP-Gallery.")
      end
    end


    # A delay eigenvalue problem
    function dep0(n::Int=5)
        Random.seed!(0) # reset the random seed
        A0 = randn(n,n)
        A1 = randn(n,n)
        tau = 1.0
        nep = DEP([A0,A1],[0,tau])
        return nep
    end


    #  A delay eigenvalue problem with one eigenvalue equal to one.
    function dep1()
        A0=([1 2 3 ; 4 5 6; 1 -1 3]);
        A1=((-A0+[1 0 3;0 0 -1;0 0 10])*exp(1));
        Q=[1 0 3; 1 1 -4; 2 3 1];
        A0=Q\(A0*Q);
        A1=Q\(A1*Q);
        nep=DEP([A0,A1],[0,1])
        return nep;
    end


    # A delay eigenvalue problem with sparse matrices
    function dep0_sparse(n::Integer=100, p::Real=0.25)
        Random.seed!(0) # reset the random seed
        A0 = sparse(1:n,1:n,rand(n))+sprand(n,n,p)
        A1 = sparse(1:n,1:n,rand(n))+sprand(n,n,p)
        tau = 1.0
        nep = DEP([A0,A1],[0,tau])
        return nep
    end


    # A delay eigenvalue problem with sparse tridiagonal matrices
    function dep0_tridiag(n::Integer=100)
        Random.seed!(1) # reset the random seed
        K = [1:n;2:n;1:n-1]; J=[1:n;1:n-1;2:n] # sparsity pattern of tridiag matrix
        A0 = sparse(K, J, rand(3*n-2))
        A1 = sparse(K, J, rand(3*n-2))
        tau = 1.0
        nep = DEP([A0,A1],[0,tau])
        return nep
    end


    # A symmetric delay eigenvalue problem with double eigenvalues
    # Examle from H. Voss and M. M. Betcke, Restarting iterative projection methods for Hermitian nonlinear eigenvalue problems with minmax property, Numer. Math., 2017
    function dep_symm_double(n::Integer=100)
        LL=-sparse(1:n,1:n,2*ones(n))+sparse(2:n,1:n-1,ones(n-1),n,n)+sparse(1:n-1,2:n,ones(n-1),n,n)

        x = range(0, stop = pi, length = n)
        h=x[2]-x[1];
        h=pi
        LL=LL/(h^2)
        LL=kron(LL,LL)

        b=broadcast((x,y)->100*abs(sin(x+y)),x,transpose(x))
        a=broadcast((x,y)->-8*sin(x)*sin(y),x,transpose(x))
        B=sparse(1:n^2,1:n^2,b[:])
        A=LL+sparse(1:n^2,1:n^2,a[:])

        nep=DEP([A,B],[0,2.0])
        return nep
    end


    # A delay eigenvalue problem with a double non-semisimple eigenvalue in λ=3πi
    # Examle from E. Jarlebring, Convergence factors of Newton methods for nonlinear eigenvalue problems, LAA, 2012
    function dep_double()
        denom = 8+5*pi;
        a1 = 2/5 *(65*pi + 32)/(denom);
        a2 = 9*pi^2*(13+5*pi)/(denom);
        a3 = 324/5 *pi^2*(5*pi+4)/(denom);
        b1 = (260*pi + 128 + 225*pi^2)/(10*denom);
        b2 = 45*pi^2/denom;
        b3 = 81*pi^2*(40*pi + 32 + 25*pi^2)/(10*denom);
        A0 = [ 0    1    0;  0    0    1;  -a3  -a2  -a1];
        A1 = [ 0    0    0;  0    0    0;  -b3  -b2  -b1];

        tau=1.0;
        nep=DEP([A0,A1],[0,tau])
        return nep
    end


    # A polynomial eigenvalue problem
    function pep0(n::Integer=200)
        Random.seed!(0)
        A0=randn(n,n)
        A1=randn(n,n)
        A2=randn(n,n)
        A=[A0,A1,A2]
        nep=PEP(A)
        return nep
    end


    # A polynomial eigenvalue problem with symmetric matrices
    function pep0_sym(n::Integer=200)
        Random.seed!(0)
        A0 = Symmetric(randn(n,n))
        A1 = Symmetric(randn(n,n))
        A2 = Symmetric(randn(n,n))
        A = [A0]#, A1, A2]
        nep = PEP(A)
        return nep
    end


    # A sparse polynomial eigenvalue problem
    function pep0_sparse(n::Integer=200, p::Real=0.03)
        Random.seed!(0)
        A0=sprandn(n,n,p)
        A1=sprandn(n,n,p)
        A2=sprandn(n,n,p)
        A=[A0,A1,A2]
        nep=PEP(A)
        return nep
    end


    #  A quadratic problem with real eigenvalues
    # Four smallest eigenvalues of the problem:
    # 1.0e+03 * [ -2.051741417993845  -0.182101627437811  -0.039344930222838  -0.004039879577113]
    function real_quadratic()
        A0 = [ 4     0     1     1;
               0     2     1     1;
               1     1     6    -2;
               1     1    -2     3];
        A1 = [ 167  -140    95  -131;
              -140   327    54    85;
                95    54   235   -81;
              -131    85   -81   181];
        A2 = [ 2     1    -1    -1;
               1     5    -3     2;
              -1    -3     3     0;
              -1     2     0     3];
        A=[A0,A1,A2]
        nep=PEP(A)
        return nep
    end


    # A quadratic delay eigenvalue problem in "The infinite Bi-Lanczos method for nonlinear eigenvalue problems",  Sarah W. Gaaf and Elias Jarlebring
    function qdep0()
        qdepbase=joinpath(dirname(@__FILE__()),
                          "gallery_extra", "qdep_infbilanczos_")
        A0=read_sparse_matrix(qdepbase * "A0.txt")
        A1=read_sparse_matrix(qdepbase * "A1.txt")
        tau::Float64 = 1
        quadfun = S -> S^2
        constfun = S -> one(S)
        expfun = S -> exp(-tau*S)

        AA = [-one(A0), A0, A1]
        fi = [quadfun, constfun, expfun]
        return SPMF_NEP(AA, fi)
    end


    # A quadratic delay eigenvalue problem in "A linear eigenvalue algorithm for the  nonlinear eigenvalue problem",  Elias Jarlebring, Wim Michiels, Karl Meerbergen \\
    function qdep1()
        A0=[ 0.3000   -0.6000         0    0.4000
            -0.3000    0.4000   -0.8000    1.9000
             0.1000   -1.6000   -1.3000         0
            -1.4000   -0.9000    0.2000    0.9000];
        A1=[ 0.8000    0.2000   -1.3000   -0.3000
            -1.1000    0.9000    1.2000    0.5000
             0.5000    0.2000   -1.6000   -1.3000
             0.7000    0.4000   -0.4000         0];
        return SPMF_NEP([one(A0), A0, A1], [λ -> -λ^2, λ -> one(λ), λ -> exp(-λ)])
    end


    # A quadratic eigenvalue problem with chosen eigenvalues
    function qep_fixed_eig(n::Integer=5)
        Random.seed!(0) # reset the random seed
        E=randn(2*n)
        return qep_fixed_eig(n,E)
    end
    function qep_fixed_eig(n::Integer=5, E::AbstractVecOrMat=randn(2*n))
        A1 = diagm(0 => E[1:n])
        A2 = diagm(0 => E[n+1:2*n])
        K = one(A1)
        nep=PEP([A1*A2,-A1-A2,K])
        return nep
    end

    # This problem stems from
    # L. P. Shayer and S. A. Campbell.  Stability, bifurcation and multistability in a system of two coupled neurons with multiple time delays. SIAM J. Applied Mathematics , 61(2):673–700, 2000
    # It is also a benchmark example in DDE-BIFTOOL
    function neuron0()
        pars = [1/2; -1; 1; 2.34; 0.2; 0.2 ; 1.5] .+ 0im
        kappa = pars[1]
        beta = pars[2]
        A = [0 pars[3]; pars[4] 0]

        x = [0; 0]  # The zero (trivial) stationary solution
        # A non-trivial stationary solution
        #x=[3.201081590416643561697725111745656884148241428177442574927999582405266342752249e-01
        #   5.096324796647208606096018689631125587762848405395086474417800152349531876959548e-01]

        tauv = [0;0.2;0.2;1.5]
        A0 = -kappa * Matrix(1.0I, 2, 2)
        A1 = A[2,1] * [0 0; (1-tanh(x[2])^2) 0]
        A2 = A[1,2] * [0 (1-tanh(x[1])^2); 0 0]
        A3 = beta * diagm(0 => [(1-tanh(x[1])^2), (1-tanh(x[2])^2)])
        return DEP([A0, A1, A2, A3], tauv)
    end

end
