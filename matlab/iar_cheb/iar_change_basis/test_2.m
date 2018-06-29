close all
clear all
clc

a=-10;   b=10;    m=6;
P=P_mat(m+1,2/(b-a),(a+b)/(a-b));
Pinv=Pinv_mat(m+1,2/(b-a),(a+b)/(a-b));
P


D=diag(diag(P));
PP=D\P;
PP
P-PP
[V,D] = eig(P);

norm(P-V*D*inv(V))