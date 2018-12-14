# adds a triangle to mesh where a b c d are the relative index coordinates
macro addtri(a,b,c,d)
    return esc(quote
               P2 = center;
               P1 = zeros(3); P1[fixdim[l]]=fixdim_val[l]; P3 = copy(P1);
               P1[freedim[1]]=grid[ii+$a]; P1[freedim[2]]=grid[jj+$b];
               P3[freedim[1]]=grid[ii+$c]; P3[freedim[2]]=grid[jj+$d];
               push!(mesh,Triangle(P1,P2,P3,area=area))
    end)
end
# generates a mesh of unit cube withe a fichera corner. N determines the size.
function gen_ficheramesh(N)

    if (mod(N, 2) != 0)   # Enforce N to be even
	N = N + 1;
    end

    nn = round(Int,N / 2); # Half the size

    # Prepare the mesh: Vector of Triangle
    mesh = Vector{Triangle}(undef,0);

    # compute the area of each triangle
    area = 0.25/N/N;

    # Uniform grid on [0,1]
    grid = Vector(0:N) / N;


    # Metadata for the discretization
    fixdim=[1;2;3;1;1;1;1;2;2;2;2;3;3;3;3]
    fixdim_val=[0;0;0;1;1;1;0.5;1;1;1;0.5;1;1;1;0.5]
    freedims=[2 3; 3 1; 1 2; 2 3; 2 3; 2 3 ; 2 3;3 1;3 1; 3 1; 3 1; 1 2; 1 2; 1 2; 1 2];
    Nvals=[1 N 1 N;
           1 N 1 N;
           1 N 1 N;
           1 nn 1 nn;
           nn+1 N 1 nn;
           1 nn nn+1 N;
           nn+1 N nn+1 N
           1 nn 1 nn
           nn+1 N 1 nn
           1 nn nn+1 N
           nn+1 N nn+1 N
           1 nn 1 nn
           nn+1 N 1 nn
           1 nn nn+1 N
           nn+1 N nn+1 N]
    # Ordering:
    # left face
    # front
    # bottom
    # right
    # (1/2, 1) x (0, 1/2)
    # (0, 1/2) x (1/2, 1)
    # (1/2, 1) x (1/2, 1)
    # back face: (0, 1/2) x (0, 1/2)
    # (1/2, 1) x (0, 1/2)
    # (0, 1/2) x (1/2, 1)
    # (1/2, 1) x (1/2, 1)
    # top face: (0, 1/2) x (0, 1/2)
    # (1/2, 1) x (0, 1/2)
    # (0, 1/2) x (1/2, 1)
    # (1/2, 1) x (1/2, 1)

    for l=1:size(fixdim,1);
        for ii = Nvals[l,1]:Nvals[l,2]
	    for jj = Nvals[l,3]:Nvals[l,4]

                center = zeros(3);
                center[fixdim[l]]=fixdim_val[l];         # Fix dimension
                freedim=freedims[l,:];
                # freedim now contains freedim[1] and freedim[2] which
                # correspond to the non-fixed dimensions
                center[freedim[1]]=(grid[ii]+grid[ii+1])/2;
                center[freedim[2]]=(grid[jj]+grid[jj+1])/2;

                # Now create the four triangles associated with this "center"
                if (l<4)
                    @addtri(0,0,1,0);
                    @addtri(1,0,1,1);
                    @addtri(1,1,0,1);
                    @addtri(0,1,0,0);
                else # For constency with Effenberger
                    @addtri(0,0,0,1);
                    @addtri(0,1,1,1);
                    @addtri(1,1,1,0);
                    @addtri(1,0,0,0);
	        end
            end
        end
    end

    # Auxiliary mesh data
    for tri = mesh
        # determine midpoint of triangle
        tri.aux.midpoint = (tri.P1 + tri.P2 + tri.P3) / 3;

        # normalized tangent vectors
        tri.aux.tau1 = normalize(tri.P2 - tri.P1);
        tri.aux.tau2 = normalize(tri.P3 - tri.P2);
        tri.aux.tau3 = normalize(tri.P1 - tri.P3);

        # exterior (normalized) normal of triangle
        tri.aux.normal = normalize(cross( tri.aux.tau1, tri.aux.tau2 ));

        # exterior (normalized) normals of the edges
        tri.aux.nu1 = normalize(cross(tri.aux.tau1, tri.aux.normal));
        tri.aux.nu2 = normalize(cross(tri.aux.tau2, tri.aux.normal));
        tri.aux.nu3 = normalize(cross(tri.aux.tau3, tri.aux.normal));
    end
    return mesh
end
