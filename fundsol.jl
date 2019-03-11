
using PyPlot;
clf()
phiv=range(-2*pi,stop=2*pi,length=100)
sx=sin.(phiv);
sy=cos.(phiv);
sz=phiv/10;
S=[sx sy sz];
plot3D(sx,sy,sz,"*")
n=length(phiv);

# Points on unit sphere
X=randn(n,3);
map(i->normalize!(view(X,i,1:3)),1:size(X,1));
X=X/2;

plot3D(X[:,1],X[:,2],X[:,3],"ro");

G=(λ,x,s) ->  exp(-1im*λ*norm(x-s))/(4*pi*norm(x-s));
Gp=(λ,x,s) -> (-1im*λ*exp(-1im*λ*norm(x-s))*4*pi*norm(x-s)-4*pi*exp(-1im*λ*norm(x-s)))/(4*pi*norm(x-s))^2;
K=(λ,x,s) -> Gp(λ,x,s)*norm(x)/norm(x-s);
#Kp=(λ,x,s) ->  -1im*exp(-1im*λ*norm(x-s))/(4*pi);


Kold=(λ,x,s) ->  exp(-1im*λ*norm(x-s))/(4*pi*norm(x-s));
Kpold=(λ,x,s) ->  -1im*exp(-1im*λ*norm(x-s))/(4*pi);

λ=3

function Mder(λ,der)
    M=zeros(ComplexF64,n,n);
    for i=1:n
        for j=1:n
            if (der==0)
                M[i,j]=K(λ,X[i,:],S[j,:]);
            else
#                M[i,j]=Kp(λ,X[i,:],S[j,:]);
            end
        end
    end
    return M
end
nep=Mder_NEP(n,Mder);
(λ,v)=broyden(nep,pmax=2,displaylevel=1)

#quasinewton(nep,λ=-3,displaylevel=1,armijo_factor=0.1,maxit=500)
