export inpolygon

"""
Return whether the given point (px, py) is inside or on the edge of the
polygon defined by vertices (polyx, polyy). Uses Hormann-Agathos
"Point in Polygon" algorithm to support arbitrary polygons.
Implementation inspired by 'isinside' in Luxor
(see https://github.com/JuliaGraphics/Luxor.jl/blob/master/src/polygons.jl)
"""
function inpolygon(px, py, polyx, polyy)
    if !isfinite(px) || !isfinite(py)
        return false
    end

    c = false
    @inbounds for idx in 1:length(polyx)
        q1x = polyx[idx]
        q1y = polyy[idx]

        nextidx = 1 + mod(idx, length(polyx))
        q2x = polyx[nextidx]
        q2y = polyy[nextidx]

        if q1x == px && q1y == py
            return true # on vertex
        end
        if q2y == py
            if q2x == px
                return true # on vertex
            elseif q1y == py && (q2x > px) == (q1x < px)
                return true # on edge
            end
        end
        if (q1y < py) != (q2y < py) # crossing
            if q1x >= px
                if q2x > px
                    c = !c
                else
                    det = det3p(q1x, q1y, q2x, q2y, px, py)
                    if isapprox(0, det)
                        return true # on edge
                    elseif (det > 0) == (q2y > q1y)
                        c = !c
                    end
                end
            elseif q2x > px
                det = det3p(q1x, q1y, q2x, q2y, px, py)
                if isapprox(0, det)
                    return true # on edge
                elseif (det > 0) == (q2y > q1y)
                    c = !c
                end
            end
        end
    end
    return c
end

function det3p(q1x, q1y, q2x, q2y, px, py)
    (q1x - px) * (q2y - py) - (q2x - px) * (q1y - py)
end
