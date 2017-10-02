module GalleryNLEVP
    # A module which handles interaction with the berlin
    # manchester NLEVP collection

    using MATLAB
    using NEPCore
    using NEPTypes    
    using Gallery    
    
    export nlevp_make_native
    
    # We have to explicitly specify functions that we want "overload"
    import NEPCore.compute_Mder
    import NEPCore.size
    
    export NLEVP_NEP

    import Gallery.nep_gallery
    export nep_gallery

    

    """
         NLEVP_NEP represents a NEP in the NLEVP-toolbox
    Example usage: nep=NLEVP_NEP("gun")
    """
    type NLEVP_NEP <: NEP
        n::Integer
        name::String
        Ai::Array
        function NLEVP_NEP(name, nlevp_path)
            if (~isfile(joinpath(nlevp_path,"nlevp.m")))
                error("nlevp.m not found. You need to install the Berlin-Manchester collection (http://www.maths.manchester.ac.uk/our-research/research-groups/numerical-analysis-and-scientific-computing/numerical-analysis/software/nlevp/) and specify a nlevp_path. Path:"*nlevp_path)
            end

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
        l::Complex128=Complex128(lambda)  # avoid type problems
        mat"""
    C=cell($maxder+1,1);
    [C{:}]=funs($l);
    $D=cell2mat(C);
    """
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

    """
    nep_gallery(NLEVP_NEP, name)
    nep_gallery(NLEVP_NEP, name, nlevp_path)
Loads a NEP from the Berlin-Manchester collection of nonlinear
eigenvalue problems.
"""
    nep_gallery{T<:NLEVP_NEP}(::Type{T},name::String) = nep_gallery(T,name,joinpath(@__DIR__, "..","..","..","nlevp3"))
    function nep_gallery{T<:NLEVP_NEP}(::Type{T},name::String,nlevp_path::String)
        if(!isabspath(nlevp_path)) #Check if path is relative, then make absoulte
            nlevp_path = abspath(nlevp_path)
        end
        nep=NLEVP_NEP(name,nlevp_path)
        return nep
    end


    """
   nlevp_make_native(nep::NLEVP_NEP)

Tries to convert the NLEVP_NEP a NEP of NEP-PACK types
"""
    function nlevp_make_native(nep::NLEVP_NEP)
        if (nep.name == "gun")
            minusop= S-> -S
            oneop= S -> eye(size(S,1),size(S,2))
            sqrt1op= S -> 1im*sqrtm(full(S))
            sqrt2op= S -> 1im*sqrtm(full(S)-108.8774^2*eye(S))
            nep2=SPMF_NEP(nep.Ai,[oneop,minusop,sqrt1op,sqrt2op])
            return nep2
        elseif (nep.name == "cd_player")
            return PEP(nep.Ai);
        elseif (nep.name == "fiber")
            error("Native implementation of fiber not finished")
        else
            error("Unable to make NEP native")
        end

    end

    # # In progress
    #function fiber2(S)
    #    # poor-mans version of
    #    f0=z -> z.*(- besselk.(1, z)./z - besselk.(0, z))./besselk.(1,z)
    #    f1=z -> 
    #    zsamples=eps()+1+cos(linspace(0,pi,100))
    #
    #    a=eps();
    #    b=3;
    #    α=-(a+b)/(b-a)
    #    β=2/(b-a);
    #
    #    
    #end
    #
end
