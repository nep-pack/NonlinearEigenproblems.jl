workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using Gallery
using PyPlot
using PyCall
76i8yuiytui
n=4;
A0=[0.3000   -0.6000         0    0.4000
   -0.3000    0.4000   -0.8000    1.9000
    0.1000   -1.6000   -1.3000         0
   -1.4000   -0.9000    0.2000    0.9000];

A1=[0.8000    0.2000   -1.3000   -0.3000
   -1.1000    0.9000    1.2000    0.5000
    0.5000    0.2000   -1.6000   -1.3000
    0.7000    0.4000   -0.4000         0];


# function myexpm(A::Array{T,2}) where {T<:Number}
#     println("call expm with",typeof(A),"\n")    A=Array{Complex128,2}(A);
#     F=zeros(T,size(A,1),size(A,2))
#     if (size(A)==(1,1))
#         F[:]=exp(A[1,1]);
#         return F
#     end
#     Bi=eye(T,size(A,1),size(A,2))
#     for k=0:10
#         F=F+Bi/factorial(real(T(k)));
#         Bi=Bi*A;
#     end
#     #F=Array{Complex128,2}(F);
#     err=norm(expm(A)-F,1)/norm(F,1);
#     if(err>eps()*100)
#         println("Warning: error large:",err, " size:",size(A), " norm(A):",norm(A));
#
#     end
#
#     return F
# end

function myexpm(A)
    F=expm(A);
    return F
end

nep=SPMF_NEP([eye(4), A0, A1],[λ->-λ^2,λ->eye(λ),λ->myexpm(-λ)])


k=3; n=4;
S=rand(k,k);    V=rand(n,k);
MM=compute_MM(nep,S,V);

λ=rand()+rand()*im;
z=compute_Mlincomb(nep,λ,V)
