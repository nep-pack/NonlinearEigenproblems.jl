# Tests
using Random;
using LinearAlgebra
using NonlinearEigenproblems;

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






@include_neptest_file("systematic_testing/dep0.jl")
@include_neptest_file("systematic_testing/sine.jl")
@include_neptest_file("systematic_testing/neuron0.jl")
@include_neptest_file("systematic_testing/beam.jl")
@include_neptest_file("systematic_testing/dep0_sparse.jl")
@include_neptest_file("systematic_testing/qdep1.jl")
@include_neptest_file("systematic_testing/schrodinger_movebc.jl")
@include_neptest_file("systematic_testing/nlevp_native_cd_player.jl")
@include_neptest_file("systematic_testing/hadeler.jl")
@include_neptest_file("systematic_testing/dep_distributed.jl")
@include_neptest_file("systematic_testing/WEP.jl")



algtol=-1;
testtol=-1;
tol0=-1;
name="acoustic_wave_1d";
λ0=-0.18577070376166677 + 0.41962694087391367im;
@include_neptest_file("systematic_testing/nlevp_problem.jl")
name="acoustic_wave_2d";
λ0=-0.6771810313836967 + 0.0897217725561532im;
@include_neptest_file("systematic_testing/nlevp_problem.jl")
name="bicycle"; λ0= -0.3228664290041067;
@include_neptest_file("systematic_testing/nlevp_problem.jl")
name="bilby"; λ0=0.23381573761990057;
@include_neptest_file("systematic_testing/nlevp_problem.jl")
name="butterfly"; λ0=-0.9703704498578221 - 1.0017769654495365im;
@include_neptest_file("systematic_testing/nlevp_problem.jl")
name="closed_loop"; λ0=-0.33107672343097816 ;
@include_neptest_file("systematic_testing/nlevp_problem.jl")
name="concrete"; λ0=44im; tol0=1e-7;
@include_neptest_file("systematic_testing/nlevp_problem.jl")
name="damped_beam"; λ0=7.42298-72.2im; tol0=1e-7
@include_neptest_file("systematic_testing/nlevp_problem.jl")
name="dirac"; λ0=-1.145+0.06im;
@include_neptest_file("systematic_testing/nlevp_problem.jl")




struct MethodTest
    evalstr::String;
    skipnep::Vector;
    results::Vector;
end
function MethodTest(s::String)
    return MethodTest(s,Vector(),Vector());
end
function MethodTest(s::String,skiptest::String)
    return MethodTest(s,Vector([skiptest]),Vector());
end
function MethodTest(s::String,skiptest1::String,skiptest2::String)
    return MethodTest(s,Vector([skiptest1,skiptest2]),Vector());
end



methodlist=Vector();
push!(methodlist,MethodTest("augnewton(nep,v=v0,λ=λ0,tol=algtol,displaylevel=d)"));
#push!(methodlist,MethodTest("resinv(nep,v=v0,λ=λ0,tol=algtol,displaylevel=d)","sine"));
#push!(methodlist,MethodTest("implicitdet(nep,v=v0,λ=λ0,tol=algtol,displaylevel=d)"));
#push!(methodlist,MethodTest("quasinewton(nep,v=v0,λ=λ0,tol=algtol,displaylevel=d)"));
#push!(methodlist,MethodTest("broyden(nep,σ=λ0,pmax=1,tol=algtol)"));
#push!(methodlist,MethodTest("newtonqr(nep,v=v0,λ=λ0,tol=algtol)"));
#push!(methodlist,MethodTest("blocknewton(nep,S=reshape([λ0],1,1), X=reshape(v0,size(v0,1),1),tol=algtol)","WEP"))
#push!(methodlist,MethodTest("iar(nep,v=v0,σ=λ0,Neig=1,tol=algtol)"));
#push!(methodlist,MethodTest("nleigs(nep,region,tol=algtol)","WEP"));



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
            catch e
                error=true;
                push!(m.results,e);
            end

            if (error || size(λ,1)==0)
                println("(E)");
            else
                λ=λ[1]; v=v[:,1];
                v=v/norm(v);
                λerr=abs(λ-test.λ)/abs(test.λ);
                verr=norm(vec(v)-test.v);
                if (λerr < test.testtol)

                    println("(+)")
                    push!(m.results,"+");
                else
                    println("(- $λerr)")
                    push!(m.results,"(- $λerr)");
                end


            end
        end


#        println("vdiff:", verr < eps()*100 )
    end
    println("");
end
