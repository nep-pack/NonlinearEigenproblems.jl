
function gallery_dep_distributed()
    # E. Jarlebring and W. Michiels and K. Meerbergen},
    # The infinite  {Arnoldi} method and an application to time-delay systems with distributed delays,
    # Delay Systems - Methods, Applications and New Trends
    # 2012
    #
    # Some correct eigenvalues:
    # -0.400236388049641 + 0.970633098237807i
    # -0.400236388049641 - 0.970633098237807i
    #  2.726146249832675 + 0.000000000000000i
    # -1.955643591177653 + 3.364550574688863i
    # -1.955643591177653 - 3.364550574688863i
    #  4.493937056300693 + 0.000000000000000i
    # -1.631513006819252 + 4.555484848248613i
    # -1.631513006819252 - 4.555484848248613i
    # -1.677320660400946 + 7.496870451838560i
    # -1.677320660400946 - 7.496870451838560i
    #
    # 

    
    N=1000
    
    A0=-eye(3);
    A1=[2.5    2.8   -0.5
        1.8    0.3    0.3
        -2.3   -1.4    3.5];
    A2=[1.7    0.7   -0.3
        -2.4   -2.1   -0.2
        2.0    0.7    0.4];
    A3=[1.4   -1.3    0.4
        1.4    0.7    1.0
        0.6    1.6    1.7];
    idop= S -> S
    oneop= S -> eye(size(S,1),size(S,2))
    f1= S -> expm(-full(S))
    f2= S -> distributed_kernel(full(S),N)
    return SPMF_NEP([A0,A1,A2,A3],[idop,oneop,f1,f2])
end    


function distributed_kernel(S,N0)
    function fS(x); return expm(x*S)*(exp((x+0.5)^2)-exp(1/4)); end

    F=zeros(eltype(S),size(S,1),size(S,2));

    h=1/N0;
    for i=1:(N0+1) # Trapezoidal rule 
        x=float((i-1)*h-1)
        if (i==1 || i==N0+1)
            F += 0.5*fS(x)*h;
        else
            F += fS(x)*h;
        end
    end
    return F
end
