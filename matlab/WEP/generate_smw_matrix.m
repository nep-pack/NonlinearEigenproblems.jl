function [ M ] = generate_smw_matrix( n, N, Linv, dd1, dd2, Pm, Pp, K, lufactorize )
%GENERATE_SMW_MATRIX generate the matrix M associate to the SMW solver
%   n is a parameter desciding the size of the problem. OBS: nz = n, and nx = n+4
%   N is the partition parameter. How many blocks in z-direction. OBS: N+4 blocks in x-direction
%   Linv is the inverse-operator of the linear matrix equation (Sylvester part)

L = n/N;             % Number of points in one dimanesion of the regions
LL = L*L;            % Number of points in the "square interior regions"
mm = (N^2 + 4*N);    % Number of elements in SMW-matrix



nz = n;
nx = n+4;

% block extract index
II=@(i) (i-1)*L+1:i*L;              % z-direction, always equal
JJ=@(j) ((j-3)*L+1:(j-2)*L) + 2;    % x-direction, different if boundary or interior
%JJ_2=@(j) (j==1)*1 + (j==2)*2 + (j==N+3)*(n+3) + (j==N+4)*(n+4);   %OBS: REPLACED BY A FUNCTION BELOW FOR JULIA COMPATIBILITY

% compute matrix M
M = zeros(mm, mm);

EEk = zeros(nz, nx);
ek = zeros(nz, 1);

for k=1:mm
    
    [i,j] = k2ij(k,N);
    
    % Ek tilde
    EEk = 0*EEk;
    if (j==1)
        EEk(II(i),JJ_2(j, n, N)) = K(II(i),JJ_2(j, n, N));
        ek = 0*ek;
        ek(II(i)) = dd1;            
        EEk(:,1) = EEk(:,1) + Pm(ek);
    elseif (j==2)
        EEk(II(i),JJ_2(j, n, N)) = K(II(i),JJ_2(j, n, N));
        ek = 0*ek;
        ek(II(i)) = dd2;            
        EEk(:,1) = EEk(:,1) + Pm(ek);
    elseif (j==N+4)
        EEk(II(i),JJ_2(j, n, N)) = K(II(i),JJ_2(j, n, N));
        ek = 0*ek;
        ek(II(i)) = dd1;            
        EEk(:,nx) = EEk(:,nx) + Pp(ek);
    elseif (j==N+3)
        EEk(II(i),JJ_2(j, n, N)) = K(II(i),JJ_2(j, n, N));
        ek = 0*ek;
        ek(II(i)) = dd2;            
        EEk(:,nx) = EEk(:,nx) + Pp(ek);
    else
        EEk(II(i),JJ(j)) = K(II(i),JJ(j));
    end
    
    % Sylvester solve of E tilde
    Fk = Linv(EEk);
    
    % Build this matrix element
    
    for kk=1:mm
        [i,j] = k2ij(kk,N);
        % evaluate the linear functional
        if((j==1)||(j==2)||(j==N+3)||(j==N+4))
            M(kk,k) = sum(sum(Fk(II(i),JJ_2(j, n, N))))/double(L);
        else
            M(kk,k) = sum(sum(Fk(II(i),JJ(j))))/double(LL);
        end
    end   
end

M = M + eye(mm);

% % Convert to LU-factorization
if(exist('lufactorize', 'var') && lufactorize == true)
  [L, U] = lu(M);
  M = [];
  M.L = L;
  M.U = U;
end


end



function [ i,j ] = k2ij( k, N )
%K2IJ 
% convert a single index k to two indeces
% reading the matrix X by rows; left to right and top to down

    j=rem(k,N+4)+(rem(k,N+4)==0)*(N+4);
    i=(k-j)/(N+4)+1;

end


function idx = JJ_2(j, n, N)
    idx = Inf;
    if(j==1)
        idx = 1;
    elseif(j==2)
        idx = 2;
    elseif(j==N+3)
        idx = (n+3);
    elseif(j==N+4)
        idx = (n+4);
    else
        idx = Inf;
    end
end
