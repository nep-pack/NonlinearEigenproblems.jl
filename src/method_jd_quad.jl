



export jd_quad

jd_quad(nep::NEP;params...) = jd_quad(Complex128,nep;params...)
function jd_quad{T}(::Type{T},
                    nep::ProjectableNEP;
                    errmeasure::Function =
                    default_errmeasure(nep::NEP),
                    tolerance=eps(real(T))*100,
                    maxit=100,
                    λ=zero(T),
                    v0=randn(size(nep,1)),
                    displaylevel=0,
                    eigsolvertype::DataType=DefaultEigSolver)


    λ::T = T(λ)
    v0 = Array{T,1}(v0)
    v::Array{T,1} = v0/norm(v0)
    n = size(nep,1)

    proj_nep = create_proj_NEP(nep)

    V::Array{T,2} = zeros(T,size(nep,1),1)
    V[:,1] = v
    u::Array{T,1} = v
    theta::T = λ

   #loop...

   for kk=1:maxit

    	err=errmeasure(theta,u)

    if (displaylevel>0)
      println("Iteration:",k," errmeasure:",err)
   	end
  	if (err< tolerance)
      	return (theta,u)
    end

	# Projected matrices
	Ap0 = 	V'*(nep.A[1]*V)
	Ap1 = 	V'*(nep.A[2]*V)
	Ap2 = 	V'*(nep.A[3]*V)


	# Create the projected polynomial problem and
	# find the eigenvalue with smallest absolute value

	Ap = [Ap0,Ap1,Ap2]
	pep_temp = PEP(Ap)
	Dc,Vc = polyeig(T,pep_temp,DefaultEigSolver)
	c = sortperm(abs.(Dc))


	theta = Dc[c[1]]
	s = Vc[:,c[1]]
	s = s/norm(s)

	u = V*s

    Mdu::Array{T,1} = compute_Mlincomb(nep,theta,u,[1],1)
    P1 = eye(T,n) - Mdu*u'/(u'*Mdu)
    P2 = eye(T,n) - u*u'

    r = compute_Mlincomb(nep,theta,u)

    MP2 = zeros(T,size(nep,1),size(P2,2))
    for ii = 1:size(P2,2)
        MP2[:,ii] = compute_Mlincomb(nep,theta,P2[:,ii],[1],0)
    end
    X = P1*MP2

	# Least squares, -pseudo_inv(X)*r
	Q,R = qr(X)
	t = -R\(Q'*r)

	#Modified Gram-Schmidt
	for ii=1:kk
	    temp = dot(V[:,ii],t)
	    t += - temp*V[:,ii]
	end
	# reorthogonalization
	for ii=1:kk
	      temp = dot(V[:,ii],t)
	      t += - temp*V[:,ii]
	end
	v = t/norm(t);


	# Update the search space V
    V = [V v]


	println("Iteration: ",kk," norm of residual:", compute_resnorm(nep,theta,u))

  end

    err=errmeasure(theta,u)
    msg="Number of iterations exceeded. maxit=$(maxit)."
    throw(NoConvergenceException(λ,v,err,msg))
end
