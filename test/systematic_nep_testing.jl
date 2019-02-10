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


tests=Vector();

rng = MersenneTwister(1234);
#
nep=nep_gallery("dep0");
(λ,v)=augnewton(Float64,nep,λ=0,v=v=ones(size(nep,1)));
v=v/norm(v);
vapprox=v+0.5*rand!(rng, zeros(size(nep,1)));
vapprox=vapprox/norm(vapprox);
λapprox=λ+0.5*λ;
neptest=NEPTestProblemSingleVec(nep,λ,v,
                                λapprox,vapprox,
                                eps()*100,eps()*100,"dep0");

push!(tests,neptest)


nep=nep_gallery("sine");
(λ,v)=augnewton(Float64,nep,λ=-4,v=ones(size(nep,1)),tol=1e-10);
v=compute_Mder(nep,λ)\v; # One extra step
v=v/norm(v);
vapprox=v+0.001*rand!(rng, zeros(size(nep,1)));
vapprox=vapprox/norm(vapprox);
λapprox=λ-0.05*λ;
neptest=NEPTestProblemSingleVec(nep,λ,v,
                                λapprox,vapprox,
                                1e-10,eps()*100,"sine");
push!(tests,neptest)


nep=nep_gallery("neuron0");
(λ,v)=augnewton(Float64,nep,λ=-4,v=ones(size(nep,1)));
v=compute_Mder(nep,λ)\v; # One extra step
v=v/norm(v);
vapprox=v+0.01*rand!(rng, zeros(size(nep,1)));
vapprox=vapprox/norm(vapprox);
λapprox=λ+0.9*λ;
neptest=NEPTestProblemSingleVec(nep,λ,v,
                                λapprox,vapprox,
                                eps()*100,eps()*100,"neuron0");
push!(tests,neptest)


nep=nep_gallery("beam");
(λ,v)=augnewton(Float64,nep,λ=1,v=ones(size(nep,1)));
v=compute_Mder(nep,λ)\v; # One extra step
v=v/norm(v);
vapprox=v+0.01*rand!(rng, zeros(size(nep,1)));
vapprox=vapprox/norm(vapprox);
λapprox=λ+0.9*λ;
neptest=NEPTestProblemSingleVec(nep,λ,v,
                                λapprox,vapprox,
                                eps()*100,eps()*100,"beam");
push!(tests,neptest)



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
push!(methodlist,MethodTest("newtonqr(nep,v=v0,λ=λ0,tol=algtol)"));
push!(methodlist,MethodTest("blocknewton(nep,S=reshape([λ0],1,1), X=reshape(v0,size(v0,1),1),tol=algtol)"))
push!(methodlist,MethodTest("iar(nep,v=v0,σ=λ0,Neig=1,tol=algtol)"));
push!(methodlist,MethodTest("nleigs(nep,region,tol=algtol)"));



for test in tests
    global nep=test.nep;
    print("test name:",test.name, "(n=",size(nep,1),")")
    for m in methodlist
        skipneplist=m.skipnep;
        II=findall(k->skipneplist[k]==test.name,1:size(skipneplist,1))

        if (size(II,1)==0)
            method=m.evalstr;
            print(" $method");
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
                print("(E)");
            else
                λ=λ[1]; v=v[:,1];
                v=v/norm(v);
                λerr=abs(λ-test.λ)/(test.λ);
                verr=norm(vec(v)-test.v);
                print((λerr < eps()*100) ? "(+)" : "(-)")
            end
        end


#        println("vdiff:", verr < eps()*100 )
    end
    println("");
end
