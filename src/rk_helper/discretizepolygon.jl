export discretizepolygon

"""
Discretize the polygon given by the (complex) entries of Z.

# Arguments
- `z`: Polygon vertices given as complex values. Pass an empty array for the unit disk.
  Pass a single value for a disk centered at that value. Pass two values for Chebyshev points.
- `include_interior_points`: Whether to return interior points.
- `npts`: Number of boundary points to return.
- `nptsint`: Minimum number of interior points (more than this may be returned).

# Return values
- `boundary_points`: Boundary points (`npts` of these), followed by the input
  polygon vertices, with an additional final vertex equal to the first vertex.
- `interior_points`: Interior points, if specified to be returned, otherwise an
  empty vector. If specified, there will be at least `nptsint` of these.
"""
function discretizepolygon(
    z::Vector{CT} = [],
    include_interior_points::Bool = false,
    npts::Int = 10000,
    nptsint::Int = 5) where CT<:Complex{<:Real}

    if isempty(z) # return unit disk
        z = [CT(0)]
    end

    if length(z) == 1 # return disk centered at z
        zz = z[1] .+ cis.(2 * pi * (1:npts)/npts)
    elseif length(z) == 2 # return Chebyshev points
        zz = (z[2]-z[1])/2 * (cos.(pi*(npts-1:-1:0)/(npts-1)).+1) .+ z[1]
    else
        z = [z[:]; z[1]] # close polygon
        L = sum(abs.(diff(z))) # length of polygon

        # 0<=alph<=1 is position of current point between z[ind] and z[ind+1]
        ind = 1
        alph = 0.0
        zz = [z[ind]]
        remL = L/npts
        while length(zz) < npts
            d = abs(z[ind+1] - z[ind])
            if (1-alph)*d < remL # can go to end of edge
                ind += 1
                remL -= (1-alph) * d
                alph = 0.0
            else # next point zz is between z[ind] and z[ind+1]
                alph += remL / d
                remL = L / npts
                push!(zz, z[ind] + alph * (z[ind+1] - z[ind]))
            end
        end
    end

    append!(zz, z)

    Z = Vector{CT}()
    if include_interior_points
        if length(z) == 2 # interval case
            xnr = 2*nptsint
            if (xnr&1) == 0 xnr += 1 end
            xpts = range(z[1], stop = z[2], length = xnr)
            Z = collect(xpts[2:2:end])
            return zz, Z
        end

        points = length(z) == 1 ? zz : z
        realz = real(points)
        imagz = imag(points)
        real_min, real_max = extrema(realz)
        imag_min, imag_max = extrema(imagz)

        iter = 0
        spacing = (real_max - real_min) / 2.0001 / sqrt(nptsint)
        while length(Z) < nptsint
            iter += 1
            if iter > 10
                error("Failed to find interior polygon points. Polygon too narrow? (Note that intervals should be given by their two endpoints only.)")
            end

            xnr = trunc(Int, (real_max - real_min) / (2*spacing))
            ynr = trunc(Int, (imag_max - imag_min) / (2*spacing))
            spacing /= sqrt(sqrt(2))
            if xnr <= 1 || ynr <= 1
                continue
            end

            xpts = range(real_min, stop = real_max, length = xnr)
            xpts = xpts[2:2:end]

            ypts = range(imag_min-eps(), stop = imag_max+eps(), length = ynr)
            ypts = ypts[2:2:end]

            Z = [x + y*im for x in xpts for y in ypts]::Vector{CT}
            Z = filter(p -> inpolygon(real(p), imag(p), realz, imagz), Z)
        end
    end

    return zz, Z
end
