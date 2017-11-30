close all
clear all
clc

n=10;
a=rand(n,1)+1;
norm(a-cheb2mon(mon2cheb(a)))/norm(a)
norm(a-mon2cheb(cheb2mon(a)))/norm(a)


semilogy(mon2cheb(a))
figure
semilogy(abs(cheb2mon(a)))



function a=cheb2mon(b)
% a is a column vector

n=length(b);
b=flipud(b);
b(end+2)=0;
bb=zeros(n+2,1);
a=zeros(n,1);


for j=1:n 
   

   
   for k=n-j+1:-1:2
       bb(k)=2*b(k+1)-bb(k+2);
   end
   
   bb(1)=b(2)-bb(3)/2;
   a(j)=b(1)-bb(2)/2;
   
   
   b=bb;
   bb=0*bb;
end

a=flipud(a);

end


function bb=mon2cheb(a)
% a is a column vector

n=length(a);
a=flipud(a);
bb=zeros(n+2,1);
b=0*bb;


for j=n:-1:1 
    b(1)=a(j)+bb(2)/2;
    b(2)=bb(1)+bb(3)/2;
    for k=3:n-j+2
        b(k)=(bb(k-1)+bb(k+1))/2;
    end
    k=n-j+3;
	b(k)=b(k-1)/2;
    
    bb=b;   
end

bb=bb(1:n);
bb=flipud(bb);
end