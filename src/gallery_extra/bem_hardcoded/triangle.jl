# Extra variables of a triangle, that can be computed from the Triangle
mutable struct TriangleAux
    # Middle of the triangle
    midpoint::Vector{Float64}
    # Normal
    normal::Vector{Float64}
    # Tangent vectors to the edges
    tau1::Vector{Float64}
    tau2::Vector{Float64}
    tau3::Vector{Float64}
    # Normal vectors for the triangles
    nu1::Vector{Float64}
    nu2::Vector{Float64}
    nu3::Vector{Float64}

    # For quadrature
    gaussP::Matrix{Float64} # Gauss points
    gaussW::Vector{Float64} # Gauss weights
end

struct Triangle #  Note: unmutable, whereas aux is mutable
    # The vertex coordinates
    P1::Vector{Float64}
    P2::Vector{Float64}
    P3::Vector{Float64}
    # Area of the element
    area::Float64
    # Auxilary info that can be computed from P1,P2,P3
    aux::TriangleAux
end

# Constructor
function Triangle(P1,P2,P3;area=0)
    aux=TriangleAux([],[],[],[],[],[],[],[],zeros(3,3),[]);
    return Triangle(P1,P2,P3,area,aux);
end
