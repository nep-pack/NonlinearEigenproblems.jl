close all
clear all
clc

n=10;
a=rand(n,1)+1;
norm(a-cheb2mon(mon2cheb(a)))/norm(a)
norm(a-mon2cheb(cheb2mon(a)))/norm(a)

