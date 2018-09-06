using Random

a=rand(10000,50000);

tic()
b=a[10:20,1:10]
b[1,1]=0
toc()

tic()
b=view(a,10:20,1:10)
b[1,1]=0
toc()
