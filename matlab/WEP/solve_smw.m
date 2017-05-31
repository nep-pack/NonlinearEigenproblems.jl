function [ X ] = solve_smw( M, C, Linv, dd1, dd2, Pm, Pp, K )
%SOLVE_SMW solve the SMW problem
%   M is the SMW matrix - it can either be given as an explicit matrix or
%      as a struct with M.L and M.U being the LU-factorization
%   C is the right hand side
%   Linv is the inverse-operator of the linear matrix equation

if(isstruct(M))
  mm = length(M.L);
else
  mm = length(M);
end
N = sqrt(mm+4)-2;      %OBS: N^2 + 4N = length(M)

nz = size(C,1);
nx = size(C,2);
L = nz/N;
LL = L*L;

% block extract index
II=@(i) (i-1)*L+1:i*L;
JJ=@(j) ((j-3)*L+1:(j-2)*L) + 2;
JJ_2=@(j) (j==1)*1 + (j==2)*2 + (j==N+3)*(nx-1) + (j==N+4)*nx;

% compute the right hand side
LinvC = Linv(C);
b=zeros(mm,1);
for k=1:mm
    [i,j]=k2ij(k,N);
    % evaluate the linear functional
    if((j==1)||(j==2)||(j==N+3)||(j==N+4))
        b(k) = sum(sum(LinvC(II(i),JJ_2(j))))/L;
    else
        b(k) = sum(sum(LinvC(II(i),JJ(j))))/LL;
    end
end
clear LinvC;


% solve for the coefficients
if(isstruct(M))
  alpha = M.U\(M.L\b);
else
  alpha = M\b;
end


% build the solution
Y = zeros(nz,nx);
ek = zeros(nz,1);
for k=1:mm
    
    [i,j] = k2ij(k,N);
    
    % Ek tilde implicitly created
    if (j==1)
        Y(II(i),1) = Y(II(i),1) + alpha(k)*K(II(i),1);
        ek = 0*ek;
        ek(II(i)) = dd1;            
        Y(:,1) = Y(:,1) + alpha(k)*Pm(ek);
    elseif (j==2)
        Y(II(i),2) = Y(II(i),2) + alpha(k)*K(II(i),2);
        ek = 0*ek;
        ek(II(i)) = dd2;            
        Y(:,1) = Y(:,1) + alpha(k)*Pm(ek);
    elseif (j==N+4)
        Y(II(i),nx) = Y(II(i),nx) + alpha(k)*K(II(i),nx);
        ek = 0*ek;
        ek(II(i)) = dd1;            
        Y(:,nx) = Y(:,nx) + alpha(k)*Pp(ek);
    elseif (j==N+3)
        Y(II(i),nx-1) = Y(II(i),nx-1) + alpha(k)*K(II(i),nx-1);
        ek = 0*ek;
        ek(II(i)) = dd2;            
        Y(:,nx) = Y(:,nx) + alpha(k)*Pp(ek);
    else
        Y(II(i),JJ(j)) = Y(II(i),JJ(j)) + alpha(k)*K(II(i),JJ(j));
    end
        
end
X=Linv(C-Y);


end



function [ i,j ] = k2ij( k, N )
%K2IJ
% convert a single index k to two indeces
% reading the matrix X by rows; left to right and top to down

    j=rem(k,N+4)+(rem(k,N+4)==0)*(N+4);
    i=(k-j)/(N+4)+1;
    

end