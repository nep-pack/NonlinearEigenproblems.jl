module GalleryNLEVP
    # A module which handles interaction with the berlin
    # manchester NLEVP collection

    using MATLAB
    using NonlinearEigenproblems.NEPCore
    using NonlinearEigenproblems.NEPTypes
    using NonlinearEigenproblems.Gallery

    # We have to explicitly specify functions that we want "overload"
    import .NEPCore.compute_Mder
    import .NEPCore.size
    import .NEPCore.compute_Mlincomb

    export NLEVP_NEP

    import .Gallery.nep_gallery
    export nep_gallery


    function fetch_nlevp_path()
        # Try to find NLEVP
        nlevp_path=joinpath(ENV["HOME"],"src","nlevp"); # default ~/src/nlevp
        try
            nlevp_path=ENV["NLEVP_PATH"]
        catch
            # Environment variables was not set
        end

        if (!isfile(joinpath(nlevp_path,"nlevp.m")))
            error("Unable to find NLEVP when looking in path=",nlevp_path," Try setting environment variable NLEVP_PATH to the directory containing nlevp.m in the zip-file in http://www.maths.manchester.ac.uk/our-research/research-groups/numerical-analysis-and-scientific-computing/numerical-analysis/software/nlevp/")
        end
        return nlevp_path;

    end


    """
         NLEVP_NEP represents a NEP in the NLEVP-toolbox
    Example usage: nep=NLEVP_NEP("gun")
    """
    struct NLEVP_NEP <: NEP
        n::Integer
        name::String
        Ai::Array
        function NLEVP_NEP(name)
            nlevp_path=fetch_nlevp_path();
            println("path:",nlevp_path);
            mat"""
            addpath($nlevp_path);
            [$Ai,funs] = nlevp($name);
        """

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
        M=zero(nep.Ai[1]);
        for i=1:length(nep.Ai)
            M=M+nep.Ai[i]*f[i]
        end
        return M
    end

    compute_Mlincomb(nep::NLEVP_NEP,λ::Number,V::Union{AbstractMatrix,AbstractVector}, a::Vector) = compute_Mlincomb_from_Mder(nep,λ,V,a)
    compute_Mlincomb(nep::NLEVP_NEP,λ::Number,V::Union{AbstractMatrix,AbstractVector}) = compute_Mlincomb(nep,λ,V, ones(eltype(V),size(V,2)))

    # Return function values and derivatives of the current matlab session "funs"
    # stemming from a previous call to [Ai,funs]=nlevp(nepname).
    # The returned matrix containing derivatives has (maxder+1) rows
    function call_current_fun(lambda,maxder::Integer=0)
        l::ComplexF64=ComplexF64(lambda)  # avoid type problems
        mat"""
    C=cell($maxder+1,1);
    [C{:}]=funs($l);
    $D=cell2mat(C);
    """
        return D
    end




    # size for NLEVP_NEPs
    function size(nep::NLEVP_NEP)
        return (nep.n,nep.n)
    end
    function size(nep::NLEVP_NEP,dim)
        return nep.n
    end

    """
    nep_gallery(NLEVP_NEP, name)
    nep_gallery(NLEVP_NEP, name, nlevp_path)
Loads a NEP from the Berlin-Manchester collection of nonlinear
eigenvalue problems.
"""
    function nep_gallery(::Type{T},name::String) where {T<:NLEVP_NEP}
        nep=NLEVP_NEP(name)
        return nep
    end

end
