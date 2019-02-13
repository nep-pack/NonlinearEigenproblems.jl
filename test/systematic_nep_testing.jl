# Tests
using Random;
using LinearAlgebra


struct NEPTestProblemSingleVec
    nep::NEP;
    λ::Number
    v::Vector;
    λapprox::Number;
    vapprox::Vector;
    algtol::Float64
    testtol::Float64;
    name::String;
end

macro include_neptest_file(filename)
    include(filename)
    neptest=NEPTestProblemSingleVec(nep,λ,v,
                                    λapprox,vapprox,
                                    algtol,testtol,name);
    push!(tests,neptest)
end



tests=Vector();

rng = MersenneTwister(1234);
#
include("systematic_testing/dep0.jl")
neptest=NEPTestProblemSingleVec(nep,λ,v,
                                λapprox,vapprox,
                                algtol,testtol,name);
push!(tests,neptest)


include("systematic_testing/sine.jl")
neptest=NEPTestProblemSingleVec(nep,λ,v,
                                λapprox,vapprox,
                                algtol,testtol,name);
push!(tests,neptest)


include("systematic_testing/neuron0.jl")
neptest=NEPTestProblemSingleVec(nep,λ,v,
                                λapprox,vapprox,
                                algtol,testtol,name);
push!(tests,neptest)

include("systematic_testing/beam.jl")
neptest=NEPTestProblemSingleVec(nep,λ,v,
                                λapprox,vapprox,
                                algtol,testtol,name);
push!(tests,neptest)

include("systematic_testing/dep0_sparse.jl")
neptest=NEPTestProblemSingleVec(nep,λ,v,
                                λapprox,vapprox,
                                algtol,testtol,name);
push!(tests,neptest)


@include_neptest_file("systematic_testing/qdep1.jl")
@include_neptest_file("systematic_testing/schrodinger_movebc.jl")
@include_neptest_file("systematic_testing/nlevp_native_cd_player.jl")
@include_neptest_file("systematic_testing/hadeler.jl")

#
#
#
#
#nep=nep_gallery("schrodinger_movebc",10);
#(λ,v)=augnewton(Float64,nep,λ=1,v=ones(size(nep,1)));
#v=compute_Mder(nep,λ)\v; # One extra step
#v=v/norm(v);
#vapprox=v+0.01*rand!(rng, zeros(size(nep,1)));
#vapprox=vapprox/norm(vapprox);
#λapprox=λ+0.01*λ;
#neptest=NEPTestProblemSingleVec(nep,λ,v,
#                                λapprox,vapprox,
#                                1e-9,1e-10,"schrodinger_movebc");
#push!(tests,neptest)
#
#nep=nep_gallery("nlevp_native_cd_player");
#(λ,v)=augnewton(Float64,nep,λ=-0.027,v=ones(size(nep,1)),tol=1e-11)
#v=compute_Mder(nep,λ)\v; # One extra step
#v=v/norm(v);
#vapprox=v+0.01*rand!(rng, zeros(size(nep,1)));
#vapprox=vapprox/norm(vapprox);
#λapprox=λ+0.005*λ;
#neptest=NEPTestProblemSingleVec(nep,λ,v,
#                                λapprox,vapprox,
#                                1e-9,1e-10,"cd_player");
#push!(tests,neptest)
#
#
#nep=nep_gallery("dep_distributed");
#(λ,v)=augnewton(nep,λ=-0.4+1im,v=ones(size(nep,1)))
#v=compute_Mder(nep,λ)\v; # One extra step
#v=v/norm(v);
#vapprox=v+0.01*rand!(rng, zeros(size(nep,1)));
#vapprox=vapprox/norm(vapprox);
#λapprox=λ+0.01*λ;
#neptest=NEPTestProblemSingleVec(nep,λ,v,
#                                λapprox,vapprox,
#                                eps()*100,1e-10,"dep_distributed");
#push!(tests,neptest)
#
#
#nep=nep_gallery(WEP,nx=3,nz=3);
#(λ,v)=augnewton(nep,λ=-2.8-3.1im,v=ones(size(nep,1)),tol=1e-14)
#v=v/norm(v);
#vapprox=v+0.03*rand!(rng, zeros(size(nep,1)));
#vapprox=vapprox/norm(vapprox);
#λapprox=λ+0.1*λ;
#neptest=NEPTestProblemSingleVec(nep,λ,v,
#                                λapprox,vapprox,
#                                eps()*500,1e-10,"WEP");
#push!(tests,neptest)
#

struct MethodTest
    evalstr::String;
    skipnep::Vector;
end
function MethodTest(s::String)
    return MethodTest(s,Vector());
end
function MethodTest(s::String,skiptest::String)
    return MethodTest(s,Vector([skiptest]));
end
function MethodTest(s::String,skiptest1::String,skiptest2::String)
    return MethodTest(s,Vector([skiptest1,skiptest2]));
end



methodlist=Vector();
push!(methodlist,MethodTest("augnewton(nep,v=v0,λ=λ0,tol=algtol,displaylevel=d)",Vector()));
push!(methodlist,MethodTest("resinv(nep,v=v0,λ=λ0,tol=algtol,displaylevel=d)",Vector(["sine"])));
push!(methodlist,MethodTest("implicitdet(nep,v=v0,λ=λ0,tol=algtol,displaylevel=d)",Vector()));
push!(methodlist,MethodTest("quasinewton(nep,v=v0,λ=λ0,tol=algtol,displaylevel=d)",Vector()));
push!(methodlist,MethodTest("broyden(nep,σ=λ0,pmax=1,tol=algtol)"));
push!(methodlist,MethodTest("newtonqr(nep,v=v0,λ=λ0,tol=algtol)"));
push!(methodlist,MethodTest("blocknewton(nep,S=reshape([λ0],1,1), X=reshape(v0,size(v0,1),1),tol=algtol)"))
push!(methodlist,MethodTest("iar(nep,v=v0,σ=λ0,Neig=1,tol=algtol)"));
push!(methodlist,MethodTest("nleigs(nep,region,tol=algtol)"));



for test in tests
    global nep=test.nep;
    println("test name:",test.name, "(n=",size(nep,1),")")
    for m in methodlist
        skipneplist=m.skipnep;
        II=findall(k->skipneplist[k]==test.name,1:size(skipneplist,1))

        if (size(II,1)==0)
            method=m.evalstr;
            print("  $method");
            methodSym=Meta.parse(method);
            global λ0=test.λapprox;
            global v0=test.vapprox;
            global algtol=test.algtol;
            global d=0;

            global region=test.λ .+ abs(test.λapprox-test.λ)*[-1-1.0im; -1+1im; 1+1im; 1-1im]

            local λ,v;
            local error=false;
            try
                (λ,v)=eval(:($methodSym))
            catch
                error=true;
            end

            if (error || size(λ,1)==0)
                println("(E)");
            else
                λ=λ[1]; v=v[:,1];
                v=v/norm(v);
                λerr=abs(λ-test.λ)/abs(test.λ);
                verr=norm(vec(v)-test.v);
                println((λerr < test.testtol) ? "(+)" : "(- $λerr)")
            end
        end


#        println("vdiff:", verr < eps()*100 )
    end
    println("");
end
