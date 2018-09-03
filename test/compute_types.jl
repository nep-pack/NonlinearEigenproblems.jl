
# Run tests for the compute type tests

push!(LOAD_PATH, string(@__DIR__, "/../src"))
using NEPCore, NEPTypes, NEPSolver, Gallery
using LinearAlgebra
using Random
using Test


@testset "compute_types" begin

    # Standard PEP (Float64)
    pep=nep_gallery("pep0",200);

    # A bigfloat PEP
    A0=Matrix(I*one(BigFloat),3,3)*3
    A1=ones(BigFloat,3,3)
    pep_bigfloat=PEP([A0,A1]);

    # a complex PEP
    A0c=Matrix(I*one(ComplexF32),3,3)*3
    A1c=ones(ComplexF32,3,3)
    pep_complex=PEP([A0c,A1c]);

    # A standard DEP (Float64)
    dep0=nep_gallery("dep0",200);

    # A bigfloat DEP
    dep_bigfloat=DEP([A0,A1]);


    # Put the NEPs to test in a list and specify which "element type"
    neplist=Array{NEP,1}([pep,pep_bigfloat,pep_complex,dep0,dep_bigfloat])
    neptype_list=Vector{Any}([Float64,BigFloat,ComplexF32,Float64,BigFloat]);
    nepreal_list=Array{Bool,1}([true,true,false,true,true]);


    # All the vector types that should be tested
    typelist=[Float64,Float32,Float16,BigFloat,ComplexF16,ComplexF32,ComplexF64,Complex{BigFloat}];


    for i in 1:size(neplist,1)
        nep=neplist[i]; stype="$(typeof(nep))"
        @testset "Testing compute nep$i: $stype" begin
            nepreal=nepreal_list[i];
            n=size(nep,1);
            eltype_nep=neptype_list[i];

            displaylevel=1

            @testset "compute_Mder. NEP eltype, typeof(λ)" begin
                @testset "$eltype_nep, $Tλ" for Tλ in typelist
                    local λ::Tλ=one(Tλ);

                    M=compute_Mder(nep,λ)
                    # predicted type: greatest of typeof(λ) and eltype_nep
                    predict_type=promote_type(typeof(λ),eltype_nep)

                    @ifd(println("typeof(λ)=",typeof(λ)))
                    @ifd(println(predict_type, " =? ",eltype(M)))

                    if (nepreal)
                        @test predict_type==eltype(M)
                    else
                        @test complex(predict_type)==eltype(M)
                    end
                end
            end


            @testset "compute_Mlincomb. NEP eltype, eltype(V), typeof(λ)" begin
                @testset "$eltype_nep, $TV, $Tλ" for Tλ in typelist, TV in typelist
                    local V::Matrix{TV}=ones(TV,n,3);
                    local λ::Tλ=one(Tλ);
                    y=compute_Mlincomb(nep,λ,V)

                    # the predicted type is the "greatest" of
                    # NEP type, typeof(λ) and eltype(V)
                    predict_type=promote_type(promote_type(typeof(λ),eltype(V)),eltype_nep)

                    @ifd(println("typeof(λ)=",typeof(λ), " eltype(V)=",eltype(V)))
                    @ifd(println(predict_type, " =? ",eltype(y)))
                    if (nepreal)
                        @test predict_type==eltype(y)
                    else
                        @test complex(predict_type)==eltype(y)
                    end
                end

            end
        end




    end


end
