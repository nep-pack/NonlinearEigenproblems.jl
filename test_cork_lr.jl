using NonlinearEigenproblems

A1=rand(2,2); A2=rand(2,2);
B1=rand(2,2); B2=rand(2,2);
M=rand(2,2);  N=rand(2,2);
v1=rand(2,1); v2=rand(2,1);
w1=rand(2,1); w2=rand(2,1);
Z=rand(1,2)


Av=[A1, A2]; AvLR=[v1, v2];
Bv=[A1, A2]; BvLR=[w1, w2];
c=CORKPencilLR(M,N,Av,AvLR,Bv,BvLR,Z)
