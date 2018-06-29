%n=3;
%%M=randn(n); C=randn(n); K=randn(n);
%A0=randn(n); A1=randn(n);
%
%A0=A0; A1=A1;

n=4;
A0=[0.3000   -0.6000         0    0.4000
   -0.3000    0.4000   -0.8000    1.9000
    0.1000   -1.6000   -1.3000         0
   -1.4000   -0.9000    0.2000    0.9000];

A1=[0.8000    0.2000   -1.3000   -0.3000
   -1.1000    0.9000    1.2000    0.5000
    0.5000    0.2000   -1.6000   -1.3000
    0.7000    0.4000   -0.4000         0];
tau=1;



%%%  Compute a very exact solution 
a=-1;
b=0; k=2/(b-a) ; c=(a+b)/(a-b);
cheb_vect=@(t,Z) cos((0:(size(Z,2)-1))*acos(t))'; 
cheb2_vect_m1=@(Z)  (0:(size(Z,2)-1))';
Mterm=@(t,X) k*(X*((0:(size(X,2)-1))'.*cheb2_vect_m1(X)));

% QDEP: 
y0comp=@(X,Y) (A0+A1)\(Mterm(c,X)-A0*Y*cheb_vect(c,Y)-A1*Y*cheb_vect(-k*tau+c,Y));

[evps00,V2,H,V]=aof_abs(y0comp, a,b, n, 100);
hist2=[];
M=@(l) -l^2*eye(n)+A0+A1*exp(-tau*l);
Mp=@(l) -2*l*eye(n)-tau*A1*exp(-tau*l);
[evps0,kv]=correct_eigs(M,Mp,evps00)
evps0=evps0(~isnan(evps0))
[Y,I]=sortrows(abs(evps0)); 
evps0=evps0(I);

% evps0 contains a very exact solution.



av=-logspace(-2,2,300)
%av=[-1,logspace(0,2,50),logspace(-2,0,50)]
%av=[-1,-logspace(0,2,100),NaN,-logspace(0,-2,100)]

mineigs=[];

errv=[];
hist=[];
for a=av
    
    b=0; k=2/(b-a) ; c=(a+b)/(a-b);
    

    
    cheb_vect=@(t,Z) cos((0:(size(Z,2)-1))*acos(t))'; 
    cheb2_vect_m1=@(Z)  (0:(size(Z,2)-1))';
    Mterm=@(t,X) k*(X*((0:(size(X,2)-1))'.*cheb2_vect_m1(X)));

    % QDEP: 
    y0comp=@(X,Y) (A0+A1)\(Mterm(c,X)-A0*Y*cheb_vect(c,Y)-A1*Y*cheb_vect(-k*tau+c,Y));
    
    % DEP
    %y0comp=@(X,Y) (A0+A1)\(X*cheb_vect(c,X)-A0*Y*cheb_vect(c,Y)-A1*Y*cheb_vect(-k*tau+c,Y));
    
    
    %y0comp=@(X,Y) (A0+A1)\(sum(X,2)-A0*sum(Y,2)+A1*Y*((-1).^(1:(size(Y,2))))');
    
    
    if (length(mineigs)==0)
        %    [evps0,V2,H,V]=aof(y0comp, a,b, n, 150);
        mineigs=evps0;
        %        hist=mineigs.';
        
    end
    if (isnan(a))
        evps=evps0;
        %    evps=NaN*ones(30,1);
        hist=[hist;hist(1,:)];
    else
        [evps,V2,H,V]=aof(y0comp, a, b, n, 20);
        l=min(evps);
        [Y,I]=sortrows(abs(evps)); evps=evps(I);
        hist=addvectalignby(hist,evps.',evps0.');
        %        hist=addhistvect(hist,evps(1:20).');
    end
    
    
    % QDEP
    %shouldbezero=min(abs(eig(-l^2*eye(n)+A0+A1*exp(-l*tau))))
    % DEP
    %shouldbezero=min(abs(eig(-l*eye(n)+A0+A1*exp(-l*tau))))
end

f=figure(1);
clf;
%[Y,I]=sortrows(abs(av)'); av=av(I);
%hist=hist(1+I,:)
%for k=1:5
%    %loglog(-av(2:end),abs(errv(2:end)))
%    loglog(av,abs(hist(:,k)-hist(I(1)+1,k)),'k')
%    hold on;
%end

for k=1:size(hist,2)
    loglog(-av,abs(hist(:,k)-evps0(k)),'k')
    hold on;
end
grid on
xlabel('-a');
ylabel('|\lambda-\lambda_*|');
ax=axis();
ax(3)=eps/2; ax(4)=1;
axis(ax);


set(f,'Position',[0 0 270 240])



%% Produce eigenvalue plots


a=-1;
b=0; k=2/(b-a) ; c=(a+b)/(a-b);
cheb_vect=@(t,Z) cos((0:(size(Z,2)-1))*acos(t))'; 
cheb2_vect_m1=@(Z)  (0:(size(Z,2)-1))';
Mterm=@(t,X) k*(X*((0:(size(X,2)-1))'.*cheb2_vect_m1(X)));

% QDEP: 
y0comp=@(X,Y) (A0+A1)\(Mterm(c,X)-A0*Y*cheb_vect(c,Y)-A1*Y*cheb_vect(-k*tau+c,Y));

[evps,V2,H,V]=aof_abs(y0comp, a,b, n, 20);


f=figure(2);
clf;
p1=plot(evps0,'k*'); hold on;
p2=plot(evps,'ko'); 
axis([-8,4,-13 13])
xlabel('Real');
ylabel('Imag');
grid on;
set(f,'Position',[0 0 270 240])
legend([p1,p2],'\lambda_*','\lambda, k=20');

% Good run 
a=-1;

b=0; k=2/(b-a) ; c=(a+b)/(a-b);

cheb_vect=@(t,Z) cos((0:(size(Z,2)-1))*acos(t))'; 
cheb2_vect_m1=@(Z)  (0:(size(Z,2)-1))';
Mterm=@(t,X) k*(X*((0:(size(X,2)-1))'.*cheb2_vect_m1(X)));

% QDEP: 
y0comp=@(X,Y) (A0+A1)\(Mterm(c,X)-A0*Y*cheb_vect(c,Y)-A1*Y*cheb_vect(-k*tau+c,Y));

%[evps00,V2,H,V]=aof_abs(y0comp, a,b, n, 100);
%%evps00=evps00(abs(evps00)>0.01);
%hist2=[];
%
%M=@(l) -l^2*eye(n)+A0+A1*exp(-tau*l);
%Mp=@(l) -2*l*eye(n)-tau*A1*exp(-tau*l);
%[evps0,kv]=correct_eigs(M,Mp,evps00)
%
%evps0=evps0(~isnan(evps0))

Nv=2:100;
for N=Nv
    [evps,V2,H,V]=aof(y0comp, a, b, n, N);
    [Y,I]=sortrows(abs(evps)); evps=evps(I);
    %    hist2=addhistvect(hist2,evps(1:min(end,6)).')
    evps=evps(abs(evps)>0.01);
    %    hist2=addhistvect(hist2,evps.')    
    hist2=addvectalignby(hist2,evps.',evps0.')
end


% Bad run

a=-5;

b=0; k=2/(b-a) ; c=(a+b)/(a-b);

cheb_vect=@(t,Z) cos((0:(size(Z,2)-1))*acos(t))'; 
cheb2_vect_m1=@(Z)  (0:(size(Z,2)-1))';
Mterm=@(t,X) k*(X*((0:(size(X,2)-1))'.*cheb2_vect_m1(X)));

% QDEP: 
y0comp=@(X,Y) (A0+A1)\(Mterm(c,X)-A0*Y*cheb_vect(c,Y)-A1*Y*cheb_vect(-k*tau+c,Y));

hist3=[];

Nv=2:100;
for N=Nv
    [evps,V2,H,V]=aof(y0comp, a, b, n, N);
    [Y,I]=sortrows(abs(evps)); evps=evps(I);
    %    hist2=addhistvect(hist2,evps(1:min(end,6)).')
    evps=evps(abs(evps)>0.01);
    %    hist2=addhistvect(hist2,evps.')    
    hist3=addvectalignby(hist3,evps.',evps0.')
end



f=figure(4);
clf;

errmat3=hist3-ones(size(hist3,1),1)*evps0.';
badrun=semilogy(Nv,abs(errmat3).','r'); hold on;


errmat=hist2-ones(size(hist2,1),1)*evps0.';
%errmat=errmat(abs(errmat)<0.5)
%semilogy(Nv,abs(errmat).','k')

for k=1:size(hist2,2)
    errmat
    %    semilogy(Nv,abs(evps0(k)-hist2),'k')

    errv=errmat(:,k);
    I=find(abs(errv)>0.9);
    errv(I)=NaN*I;
    goodrun=semilogy(Nv,abs(errv),'k')    
    hold on;
end
grid on;
ax=axis();
ax(3)=eps/2; ax(4)=1;
axis(ax);
xlabel('k'); ylabel('|\lambda-\lambda_*|');


%set(badrun,'color','gray');
set(badrun,'color',[0.6,0.6,0.6]);
%set(badrun,'color','k');
%set(badrun,'LineWidth',2.5);

%for k=1:length(badrun)
%    set(badrun(k),'color',[0.6,0.6,0.6]);
%    set(badrun(k),'LineWidth',2.5);
%end

set(f,'Position',[0 0 520 240])

lg=legend([goodrun(1), badrun(1)],'a=-1','a=-5');
%set(lg,'Location','Southwest');
set(lg,'Location','Northeast');


I=find(abs(errmat3(:,1))<1e-10); bad_k=Nv(I(1))
I=find(abs(errmat(:,1))<1e-10);  good_k=Nv(I(1))


good_found=length(find(abs(errmat(79,:))<1e-10))
bad_found=length(find(abs(errmat3(79,:))<1e-10))

