function [ nep ] = load_pep( MM )
%LOAD_DEP Summary of this function goes here
%   Detailed explanation goes here

n=length(MM{1});
nep.n=n;

nep.M =@(l,v)   Meval(MM, l, v);
nep.Md =@(j)    Md_eval(MM, j)*factorial(j);

if issparse(MM{1})
    [L,U,P,Q] = lu(MM{1});
    nep.M_inv=@(v) Q*(U\(L\(P*v)));
else
    nep.M_inv=@(v) MM{1}\v;
end

nep.resnorm=@(l,v) norm(nep.M(l,v))/(n*norm(v));

end

function Md = Md_eval(MM, j)

md=length(MM)+1;
n=length(MM{1});

if j+1<md
    Md=MM{j+1};
else
    Md=sparse(n,n);
end


end

function Mv = Meval(MM, l, v)

md=length(MM);

Mv=0;
for i=1:md
    Mv=Mv+l^(i-1)*MM{i}*v;
end


end