m=10; kk=4;
λv=kk*rand(m,1) .+ kk*rand(m,1)*1im



σ=sum(λv)/length(λv);
r=maximum(abs.(σ.-λv));
θ=range(0,stop=2π,length=1000);
Σ=σ.+r*cos.(θ) + 1im*r*sin.(θ)
