function Z=compute_derivative_k(s,k)
     randn('seed',0);
     n=10;
     A0=randn(n,n); A1=randn(n,n);
     Z=zeros(n,n);
     if (k==0)
         Z=A0+s*A1;
     end
     if (k==1)
         Z=A1;
     end
     Z=Z+(A1^k)*expm(s*A1);
end
