close all
clear all
clc

A0=randn(3); A1=randn(3); A2=randn(3);
tds.A={A0,A1,A2};
tds.hA=[0,1,1.5];
[evps,VV]=tds_arnoldi_pub(tds);
s=min(evps);
M=-s*eye(3)+A0+A1*exp(-tds.hA(2)*s)+A2*exp(-tds.hA(3)*s);
shouldbezero=min(svd(M))