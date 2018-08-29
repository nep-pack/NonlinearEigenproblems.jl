workspace()

function DD0_mat_fun(T,f,S)
	# evaluete the divided differences matrix function
	# f[S,0] by using the equality
	#
	# f(S I) = (f(S) 		f[S,0]  )
	#  (0 0)   (0			f(0)	)
	#
	# Notice that f[S,0] is defined also for S singular.
	# If S is not singular it holds f[S,0]=S^(-1)-(f(S)-f(0))
	# Example:
	# n=10; S=rand(n,n); T=ComplexF64; f=x->expm(x)+x^2
	# Y1=DD0_mat_fun(T,f,S); Y2=inv(S)*(f(S)-f(zeros(S)));
	# norm(Y1-Y2)

	n=size(S,1);
	A=zeros(T,2*n,2*n);
	A[1:n,1:n]=S;
	A[1:n,n+1:2*n]=eye(S);
	return f(A)[1:n,n+1:end];
end

n=10; S=rand(n,n); T=ComplexF64; f=x->expm(x)+x^2
Y1=DD0_mat_fun(T,f,S); Y2=inv(S)*(f(S)-f(zeros(S)));
norm(Y1-Y2)
