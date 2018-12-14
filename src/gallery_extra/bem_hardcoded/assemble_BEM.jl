include("deHoop.jl");

# Precompute some guass quadrature coefficieints
function precompute_quad!(mesh,gauss_order)
    # define Gauss-Legendre quadrature points and weights
    local pt,wg;
    if (gauss_order==3)
	# points
	pt = [2/3 1/6 1/6;
	      1/6 2/3 1/6;
	      1/6 1/6 2/3];
	# weights
        wg = [1/3;
              1/3;
              1/3];
    else
        error("The Gauss quadrature order you specified is not implemented")
    end
    # pre-compute guassP and guassW for all triangles
    for tri = mesh
	VK = [tri.P1'; tri.P2'; tri.P3'];
	tri.aux.gaussP = transpose(pt * VK);
	tri.aux.gaussW = tri.area * wg;
    end
end

function assemble_BEM(lambda, mesh, gauss_order,der=0)

    # Number of triangles = size of matrix
    n = size(mesh,1);

    # Storage for the computed matrix
    T = zeros(ComplexF64,n, n);

    # prepare auxiliary vectors
    rowind = vec(repeat(1:gauss_order,1,gauss_order)'); # [1 1 1 2 2 2 3 3 3 ]
    colind = vec(repeat(1:gauss_order,1,gauss_order))  # [1 2 3 1 2 3 1 2 3]


    for row = 1:n
	@inbounds for col = row:n # Skip because of symmetry

	    # compute (row, col)-entry of T
            trir=mesh[row];
            tric=mesh[col]


            # Compute a distance vector where zero elements
            # are set to one and marked in the idx vector
            S1=trir.aux.gaussP[:, rowind];
            S2=tric.aux.gaussP[:, colind];
            s0=(S1 - S2).^2
	    dist = sqrt.(sum(s0, dims=1));
	    idx = (dist .== 0);
	    dist[idx] .= 1;


            # E = derivative of greens function
            local E
            if (der==0) # Zero'th derivative
	        E = exp.((1im * lambda) .* dist) .- 1;
	        E[vec(idx)] .= 1im * lambda; # Modify in dist-zero points
            elseif (der==1)
	        E = (1im*dist).*(exp.((1im * lambda) .* dist));
	        E[vec(idx)] .= 1im;
            else
	        E = ((1im*dist).^der).*(exp.((1im * lambda) .* dist));
	        E[vec(idx)] .= 0;
            end

            # Compute the matrix element
            mgrW=trir.aux.gaussW[rowind]
            mgcW=tric.aux.gaussW[colind]
            mgrW2=trir.aux.gaussW
            local aa=mgrW .* (mgcW ./ dist');
            T[row,col]=(E * aa)[1] / (4 * pi);
            if (der==0) # Note: this constant term disappears when we differentiate
                sing_kern=deHoop(trir.aux.gaussP, tric); # Compute singular integral
	        T[row, col] += sing_kern' * mgrW2 / (4 * pi);
            end

	    # Set element on other side of diagonal
	    T[col, row] = T[row, col];
	end
    end

    return T

end
