include("inpolygon.jl")

# ZZ = DISCRETIZEPOLYGON(Z)
# Discretize the polygon given by the (complex) entries of Z.
# npts = nr of boundary points
# nsigma = nr of interior points
function discretizepolygon(
    z::Vector{T} = [],
    include_interior_points::Bool = false,
    npts::Int = 10000,
    nsigma::Int = 5) where T<:Number

    if isempty(z) # return unit disk
        z = [T(0)]
    end

    if length(z) == 1 # return disk centered at z
        zz = z[1] + cis.(2 * pi * (1:npts)/npts)
    elseif length(z) == 2 # return Chebyshev points
        zz = (z[2]-z[1])/2 * (cos.(pi*(npts-1:-1:0)/(npts-1))+1) + z[1]
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

    Z = Vector{T}(0)
    if include_interior_points
        if length(z) == 2 # interval case
            xnr = 2*nsigma
            if (xnr&1) == 0 xnr += 1 end
            xpts = linspace(z[1], z[2], xnr)
            Z = collect(xpts[2:2:end])
            return zz, Z
        end

        realz = real.(z)
        imagz = imag.(z)
        iter = 0
        spacing = (maximum(real(z)) - minimum(real(z))) / 2.0001 / sqrt(nsigma)
        while length(Z) < nsigma
            iter += 1
            if iter > 10
                error("Failed to find interior polygon points. Sigma too narrow? (Note that intervals should be given by their two endpoints only.)")
            end

            xnr = trunc(Int, (maximum(real(z)) - minimum(real(z))) / (2*spacing))
            ynr = trunc(Int, (maximum(imag(z)) - minimum(imag(z))) / (2*spacing))
            spacing /= sqrt(sqrt(2))
            if xnr <= 1 || ynr <= 1
                continue
            end

            xpts = linspace(minimum(real(z)), maximum(real(z)), xnr)
            xpts = xpts[2:2:end]

            ypts = linspace(minimum(imag(z))-eps(), maximum(imag(z))+eps(), ynr)
            ypts = ypts[2:2:end]

            Z = [x + y*im for x in xpts for y in ypts]::Vector{T}
            Z = filter(p -> inpolygon(real(p), imag(p), realz, imagz), Z)
        end
    end

    return zz, Z
end
