
# Run tests for the compute type tests

push!(LOAD_PATH, string(@__DIR__, "/../src"))
using NEPCore, NEPTypes, NEPSolver, Gallery
using LinearAlgebra
using Random
using Test


struct compute_types_metadata
    nep::NEP
    nep_numbertype::DataType
    nep_real::Bool
    skip_types_Mder::Vector{DataType}
    skip_types_Mlincomb::Vector{DataType}
    skip_types_MM::Vector{DataType}
    check_type_stability::Bool
end


function test_one_nep(metadata::compute_types_metadata;fulltest=true)
    # All the vector types that should be tested

    typelist=[Float64,Float32,Float16,BigFloat,ComplexF16,ComplexF32,ComplexF64,Complex{BigFloat}];

    if (!fulltest)
        typelist=[typelist[1],typelist[5]];
    end

    nep=metadata.nep; stype="$(typeof(nep))"
    @testset "Testing compute nep: $stype" begin
        nepreal=metadata.nep_real;
        n=size(nep,1);
        eltype_nep=metadata.nep_numbertype;
        skip_types_Mder=metadata.skip_types_Mder;
        skip_types_Mlincomb=metadata.skip_types_Mlincomb;
        skip_types_MM=metadata.skip_types_MM;


        check_type_stability=metadata.check_type_stability;


        println("Testing compute functions for NEP:$(stype)");
        displaylevel=0

        ## Test compute_Mder
        @testset "compute_Mder. NEP:$eltype_nep, typeof(λ)" begin
            @testset "$Tλ" for Tλ in typelist
                if (!in(Tλ,skip_types_Mder)  && !in(eltype_nep,skip_types_MM))
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

                    # Check type stability as well
                    if check_type_stability
                        @inferred compute_Mder(nep,λ)
                    end

                end
            end
        end


        ## Test compute_Mlincomb
        @testset "compute_Mlincomb. NEP:$eltype_nep. eltype(V), typeof(λ)" begin
            @testset "$TV, $Tλ" for Tλ in typelist, TV in typelist
                if (!in(Tλ,skip_types_Mlincomb)  && !in(eltype_nep,skip_types_Mlincomb) &&
                    !in(TV,skip_types_Mlincomb))
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

                    if check_type_stability
                        @inferred compute_Mlincomb(nep,λ,V) # Type stability
                    end

                    # Test also when v is just a vector
                    local v::Vector{TV}=ones(TV,n);
                    y2=compute_Mlincomb(nep,λ,v);
                    if (nepreal)
                        @test predict_type==eltype(y2)
                    else
                        @test complex(predict_type)==eltype(y2)
                    end

                    if check_type_stability
                        @inferred compute_Mlincomb(nep,λ,v) # Type stability
                    end
                end
            end
        end

        ## Test compute_MM
        @testset "compute_MM. NEP:$eltype_nep. eltype(S), eltype(V)" begin
            @testset "$TS, $TV" for TS in typelist, TV in typelist
                # Exclude some tests since expm(::BigFloat) is not implemented
                if (!in(TS,skip_types_MM)  && !in(TV,skip_types_MM) && !in(eltype_nep,skip_types_MM))
                    local S::Matrix{TS}=Matrix(one(TS)*I,2,2);
                    local V::Matrix{TV}=ones(TV,n,2);
                    Y=compute_MM(nep,S,V)

                    # the predicted type is the "greatest" of
                    # NEP type, eltype(S) and eltype(V)
                    predict_type=promote_type(promote_type(eltype(S),eltype(V)),eltype_nep)

                    if (nepreal)
                        @test predict_type==eltype(Y)
                    else
                        @test complex(predict_type)==eltype(Y)
                    end
                    if check_type_stability
                        @inferred compute_MM(nep,S,V) # Type stability
                    end

                end
            end
        end

    end


end



@testset "compute_types" begin

    testlist=Vector{compute_types_metadata}()

    # Standard PEP (Float64)
    pep=nep_gallery("pep0",5);

    push!(testlist,compute_types_metadata(pep,Float64,true,[],[],[],false));

    # A bigfloat PEP
    A0=Matrix(I*one(BigFloat),3,3)*3
    A1=ones(BigFloat,3,3)
    pep_bigfloat=PEP([A0,A1]);

    push!(testlist,compute_types_metadata(pep_bigfloat,BigFloat,true,[],[],[],false));

    # a complex PEP
    A0c=Matrix(I*one(ComplexF32),3,3)*3
    A1c=ones(ComplexF32,3,3)
    pep_complex=PEP([A0c,A1c]);

    push!(testlist,compute_types_metadata(pep_complex,ComplexF32,false,[],[],[],false));

    # A standard DEP (Float64)
    dep0=nep_gallery("dep0",10);

    bigfloat_types=Vector{DataType}([BigFloat,Complex{BigFloat}])
    push!(testlist,compute_types_metadata(dep0,Float64,true,[],[],bigfloat_types,true));

    # A bigfloat DEP
    dep_bigfloat=DEP([A0,A1]);

    push!(testlist,compute_types_metadata(dep_bigfloat,BigFloat,true,[],[],bigfloat_types,true));

    A0bc=Matrix{Complex{BigFloat}}(A0);
    A1bc=Matrix{Complex{BigFloat}}(A1);
    dep_cbigfloat=DEP([A0bc,A1bc]);


    push!(testlist,compute_types_metadata(dep_cbigfloat,Complex{BigFloat},false,[],[],bigfloat_types,false));


    # SPMF 1
    B0=randn(3,3); B1=randn(3,3);
    oneop= S-> S;
    sqrop= S-> Matrix(S)^2;
    spmf_nep=SPMF_NEP([B0,B1],[oneop,sqrop])
    push!(testlist,compute_types_metadata(spmf_nep,Float64,true,[],[],[],false));

    # SPMF 2
    B0=randn(3,3); B1=randn(3,3);
    oneop= S-> 1im*S;
    sqrop= S-> Matrix(S)^2;
    spmf_nep=SPMF_NEP([B0,B1],[oneop,sqrop])
    push!(testlist,compute_types_metadata(spmf_nep,ComplexF64,false,[],[],[],false));

    # SPMF 3
    expmop= S -> exp(Matrix(S))
    B0b=Matrix{BigFloat}(B0);
    B1b=Matrix{BigFloat}(B1);
    spmf_nep2=SPMF_NEP([B0b,B1b],[oneop,expmop])

    # Skip these since expm not supported
    bigfloats_and_float16=[BigFloat,Complex{BigFloat},Float16,Complex{Float16}];

    push!(testlist,compute_types_metadata(spmf_nep2,BigFloat,true,
                                          bigfloats_and_float16,bigfloats_and_float16,
                                          bigfloats_and_float16,false));



    testlist=testlist[1:end]

    for i in 1:size(testlist,1)
        test_one_nep(testlist[i],fulltest=true)
    end


end
