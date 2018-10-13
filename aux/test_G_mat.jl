m=10;
# very compact way to generate the matrix G
G=zeros(m+1,m+1);
for i=1:m+1 G[i,1]=1/i end
for j=1:m
    for i=1:m+1
        G[i,j+1]=(G[i,j]*j)/(i+j);
    end
end
G
