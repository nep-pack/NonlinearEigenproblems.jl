"""
    Module with an implementation of CORK (CORK version 0.2)
     http://twr.cs.kuleuven.be/research/software/nleps/cork.php
"""
module NEPSolver_CORK # Will change later

    using LinearAlgebra

    export cork
    export TargetType
    export CORKInputType

    struct TargetType
        Sigma::ComplexF64
        shift::ComplexF64
    end

    struct CORKInputType
        A # Length d array with n x n matrices
        B # Length d array with n x n matrices
        M # (d-1) x d matrix
        N # (d-1) x d matrix

    end
    """

"""
#    function cork(L)
#
#    end
#
#    function cork(L,k::Integer)
#
#    end
#
#    function cork(L,k::Integer,sigma::AbstractFloat)
#
#    end
    function cork(L::CORKInputType,k::Integer,target::TargetType)

        #
        n=size(L.A[1],1)



        ##  Check inputs
        #nb = ceil( (m + maxrest*(m-p) + 1) / floor(length(shifts)) );

        nb=10; # hard coded
        shifts=repmat([target.shift],nb)

        ## Initialize

        # Default values:

        m = max(3k,20);
        maxrest = 50;
        p = max(2k,10);

        d=size(L.A,1)

        ## Initialization
        v0=ones(n)

        # Compute the Q-matrix
        Q,U0 = qr(reshape(v0/opnorm(v0),n,1));
        Q = Matrix(Q)
        r = size(Q,2);

        # Init tensor U
        #U = zeros(r,m+1,d);
        U = zeros(r+1,m+1,d); #Don't use size(U) since it is different from CORK
        U[1:r,1,1:size(U0,1)] = U0;
        j = 1;

        # Init Hessenberg matrices H and K
        H = zeros(m+1,m);
        K = zeros(m+1,m);

        # Storage of eigenpairs and residuals
        X = NaN*zeros(n,m);
        lambda = NaN*zeros(m);
        res = NaN*zeros(m);

        # History storing
        Lam = NaN*zeros(m,m + maxrest*(m-p));
        Res = NaN*zeros(m,m + maxrest*(m-p));

        # What's this for ??
        J = [1;zeros(m + maxrest*(m-p),1)];
        R = [r;zeros(m + maxrest*(m-p),1)];

        # For LU
        LU = Array{Any}(1);
        SHIFTS = [];

        # Don't know how to do cellfun in julia

        # Which matrices in L.A and L.B have non-zero elements
        NNZ_A=zeros(Bool,d)
        NNZ_B=zeros(Bool,d)
        for i=1:d
            NNZ_A[i]=countnz(L.A[i])>0
            NNZ_B[i]=countnz(L.B[i])>0
        end

        i=1
        nbrest = 0;
        l = 0;
        count = 0;
        flag = true;

        println("shifts:",shifts)
        while i <= m + maxrest*(m-p)


            println("size(Q):",size(Q))
            println("r:",r," i=",i," j=",j, " Q=",Q)
            println("U=",U)
            r = corkstep!(shifts[i],r,j,i,Q,U,L,NNZ_A,NNZ_B)
            println("r=",r)
            return
            # r = corkstep(shifts(i),r,j,i);

            #            flag = ritzpairs(r,j,i,l);
            #            if j==m
            #                if nbrest < maxrest
            #                    nbrest = nbrest + 1;
            #                    [r,j,l] = implicitrestart(r,nbrest);
            #                else
            #                    break
            #                end
            #            end

            i+=1
            j+=1
        end



    end

    function corkstep!(shift,r,j,i,Q,U,L,NNZ_A,NNZ_B)
        println("Iteration ",i)
        eta=1/sqrt(2);
        d=size(L.A,1)


        # Next vector
        v=Q[:,1:r]*reshape(U[1:r,j,1:d],r,d);


        (v1,MsN1inv0,MsN1invN)=backslash(shift,v,L,NNZ_A,NNZ_B)
        # Not yet implemented

        # level 1 orthogonalization
        q = copy(v1);
        delta = Inf; nborth = 0;
        maxreorth=2
        while (opnorm(q) < eta*delta) && (nborth <= maxreorth)
            delta = opnorm(q);
            if nborth == 0
                u1 = Q[:,1:r]'*v1;
                q = q - Q[:,1:r]*u1;
            else
                q = q - Q[:,1:r]*(Q[:,1:r]'*q);
            end
            nborth = nborth + 1;
            println("Reorth Q ",nborth)

        end
        delta = opnorm(q)+eps();
        #println("XXXXXXXXXXX")

        println("delta=",delta)
        # Update Q
        if delta > eps()
            rnew = r + 1;
            q = q/delta;
            #Q[:,rnew] = q;
            Q=hcat(Q,q)
            println("size(Q)=",size(Q))
        else
            rnew = r;
        end
        println("rnew=",rnew)
        # level 2 orthogonalization
        if rnew > r
            ## How to expand tensors? (Workaround: init with larger)
            U[rnew,:,:] = 0;
            u1 = [u1;q'*v1];
        end
        tmp=broadcast(*,transpose(MsN1inv0),u1)

        println("size(u1)=",size(u1))
        u2=reshape(U[1:rnew,j,1:d],rnew,d)*transpose(MsN1invN)-tmp

        println("*************** tmp=", tmp)
        println(" u2 =",u2)
        ##
        u = [u1;reshape(u2,size(u2,1)*size(u2,2),1)];
        println(" opnorm(u) =",opnorm(u))
        while (opnorm(u) < eta*delta) && (nborth <= maxreorth)
            delta = opnorm(u);
            h = reshape(conj(permute(U(1:rnew,1:j,1:d),[2,1,3])),[],d*rnew)*u;

        end
        rnew=1
        return rnew,Q

    end
    function backslash(shift,y,L,NNZ_A,NNZ_B)

        reuselu=false
        newlu=true


        # initialize
        d=size(L.A,1)
        m0 = L.M[1:d-1,1];
        M1 = L.M[1:d-1,2:d];
        n0 = L.N[1:d-1,1];
        N1 = L.N[1:d-1,2:d];

        # Kronecker coefficients
        MsN1 = M1 - shift*N1;
        msn0 = m0 - shift*n0;
        MsN1inv0 = MsN1\msn0;
        MsN1invN = MsN1\L.N[1:d-1,1:d];

        ## Build At
        At = L.A[1] - shift*L.B[1];
        for ii = 2:d
            At = At - MsN1inv0[ii-1]*(L.A[ii] - shift*L.B[ii]);
            println("MsN1inv0[ii-1]=",MsN1inv0[ii-1])
        end

        println("opnorm(At)=",opnorm(At))

        println(d)

        ## intermediate x used in yt
        x = y*transpose(MsN1invN);

        ## build yt
        yt = L.B[1]*y[:,1];
        for ii = 2:d
            if NNZ_A[ii] && NNZ_B[ii]
                yt = yt - L.A[ii]*x[:,ii-1] + L.B[ii]*(shift*x[:,ii-1]+y[:,ii]);
            elseif NNZ.A[ii]
                yt = yt - L.A[ii]*x[:,ii-1];
            else
                yt = yt + L.B[ii]*(shift*x[:,ii-1]+y[:,ii]);
            end
        end

        ## build x1
        println("At",size(At))
        x1 = At\yt;


        return x1,MsN1inv0,MsN1invN



    end
end
