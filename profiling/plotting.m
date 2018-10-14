close all
clear all
clc

err=importdata("error_ilan.txt");

[m,p]=size(err);

% fix error plots

for j=1:p
    err(1:m,j)=sort(err(1:m,j),'descend');
end


for j=1:p
    semilogy(1:m,err(1:m,j),'-k');
    hold on
end
ylim([1e-20 1e1])