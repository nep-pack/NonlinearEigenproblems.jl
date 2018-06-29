close all
clear all
clc

a=-1;   b=0;    m=20;
P=P_mat(m+1,2/(b-a),(a+b)/(a-b));
Pinv=Pinv_mat(m+1,2/(b-a),(a+b)/(a-b));

P