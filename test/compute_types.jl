# Run tests for the compute type tests

using NonlinearEigenproblems
using Test,Random
using LinearAlgebra
using SpecialFunctions
using SparseArrays

import Base.exp

struct compute_types_metadata
    nep::NEP
    nep_numbertype::DataType
    nep_real::Bool
    skip_types_Mder::Vector{DataType}
    skip_types_Mlincomb::Vector{DataType}
    skip_types_MM::Vector{DataType}
    check_type_stability::Bool
    name::AbstractString
end

compute_types_metadata(nep, nep_numbertype, nep_real, skip_types_Mder, skip_types_Mlincomb, skip_types_MM, check_type_stability) =
    compute_types_metadata(nep, nep_numbertype, nep_real, skip_types_Mder, skip_types_Mlincomb, skip_types_MM, check_type_stability, "")

exp(A::Matrix{Float16})=Matrix{Float16}(exp(Matrix{Float32}(A))) # Hack which makes exp(::Matrix{Float16}) available
exp(A::Matrix{ComplexF16})=Matrix{ComplexF16}(exp(Matrix{ComplexF32}(A))) # Hack which makes exp(::Matrix{Float32}) available



function test_one_nep(metadata::compute_types_metadata,typelist::Vector{DataType})
    # All the vector types that should be tested



    nep=metadata.nep;
    eltype_nep=metadata.nep_numbertype;
    name = string(typeof(nep).name) * (isempty(metadata.name) ? "" : " " * metadata.name)
    @bench @testset "$name ($eltype_nep)" begin
        nepreal=metadata.nep_real;
        n=size(nep,1);
        skip_types_Mder=metadata.skip_types_Mder;
        skip_types_Mlincomb=metadata.skip_types_Mlincomb;
        skip_types_MM=metadata.skip_types_MM;


        check_type_stability=metadata.check_type_stability;


        @info "Testing compute functions for NEP:$name ($eltype_nep)"
        displaylevel=0

        ## Test compute_Mder
        @testset "compute_Mder. NEP:$eltype_nep, typeof(λ)" begin

            @testset "$Tλ" for Tλ in typelist
                if (!in(Tλ,skip_types_Mder)  && !in(eltype_nep,skip_types_MM))
                    local λ::Tλ=one(Tλ);

                    if (!(Tλ<:Real))
                        λ += 1im*one(Tλ)
                    end

                    M=compute_Mder(nep,λ)
                    # predicted type: greatest of typeof(λ) and eltype_nep
                    predict_type=promote_type(typeof(λ),eltype_nep)

                    @debug "typeof(λ) = $(typeof(λ))"
                    @debug "$predict_type =? $(eltype(M))"

                    if (nepreal)
                        @test(predict_type==eltype(M));
                    else
                        @test(complex(predict_type)==eltype(M))
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

                    # the predicted type is the "greatest" of
                    # NEP type, typeof(λ) and eltype(V)
                    predict_type=promote_type(promote_type(typeof(λ),eltype(V)),eltype_nep)
                    y=compute_Mlincomb(nep,λ,V)

                    @debug "typeof(λ) = $(typeof(λ)), eltype(V) = $(eltype(V))"
                    @debug "$predict_type =? $(eltype(y))"

                    if (nepreal)
                        @test(predict_type==eltype(y))
                    else
                        @test(complex(predict_type)==eltype(y))
                    end

                    y=compute_Mlincomb(nep,λ,V)
                    if (nepreal)
                        @test(predict_type==eltype(y))
                    else
                        @test(complex(predict_type)==eltype(y))
                    end

                    if check_type_stability
                        @inferred compute_Mlincomb(nep,λ,V) # Type stability
                    end

                    # Test also when v is just a vector
                    local v::Vector{TV}=ones(TV,n);
                    y2=compute_Mlincomb(nep,λ,v);
                    if (nepreal)
                        @test(predict_type==eltype(y2))
                    else
                        @test(complex(predict_type)==eltype(y2))
                    end

                    if check_type_stability # Type stability
                        @inferred compute_Mlincomb(nep,λ,v)
                        @inferred compute_Mlincomb!(nep,λ,v)
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
                        @test(predict_type==eltype(Y))
                    else
                        @test(complex(predict_type)==eltype(Y))
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





    full_typelist=[Float64,Float32,Float16,BigFloat,ComplexF16,ComplexF32,ComplexF64,Complex{BigFloat}];
    reduced_typelist=[Float64,ComplexF16];



    full_test=false;

    local typelist
    if (full_test)
        typelist=full_typelist
    else
        typelist=reduced_typelist
    end



    bigfloat_types=Vector{DataType}([BigFloat,Complex{BigFloat}])

    testlist=Vector{compute_types_metadata}()


    if (full_test)
        for i=1:size(typelist,1)
            T=typelist[i];
            A0=Matrix(I*one(T),3,3)*3
            A1=ones(T,3,3)
            A2=ones(T,3,3)*3
            if (!(T<:Real))
                A0=A0+1im*ones(T,3,3)
                A1=A1+1im*ones(T,3,3)
                A2=A2+1im*Matrix(I*one(T),3,3)
            end
            A0_sparse=sparse(A0);
            A1_sparse=sparse(A1);
            A2_sparse=sparse(A1);

            dep_i=DEP([A0,A1]);
            push!(testlist,compute_types_metadata(dep_i,T,T<:Real,[],[],bigfloat_types,true));

            dep_i=DEP([A0_sparse,A1_sparse]);
            push!(testlist,compute_types_metadata(dep_i,T,T<:Real,[],[],bigfloat_types,true));

            pep_i=PEP([A0,A1,A2]);
            push!(testlist,compute_types_metadata(pep_i,T,T<:Real,[],[],[],false));

            pep_i=PEP([A0_sparse,A1_sparse,A2_sparse]);
            push!(testlist,compute_types_metadata(pep_i,T,T<:Real,[],[],[],false));


            oneop= S-> S;
            sqrop= S-> S^2;
            spmf_nep=SPMF_NEP([A0,A2],[oneop,sqrop])
            # Disable test of compute_MM + compute_Mlincomb:  since they currently fails.
            #
            push!(testlist,compute_types_metadata(spmf_nep,T,T<:Real,
                                                  [],full_typelist,full_typelist,false));


            spmf_nep_sparse=SPMF_NEP([A0_sparse,A2_sparse],[oneop,sqrop])
            # Disable test of compute_MM + compute_Mlincomb:  since they currently fails.
            #
            push!(testlist,compute_types_metadata(spmf_nep_sparse,T,T<:Real,
                                                  [],full_typelist,full_typelist,false));

        end

    else

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

        push!(testlist,compute_types_metadata(dep0,Float64,true,[],[],bigfloat_types,true));


        # A standard DEP (Float64 + sparse)
        A0_sparse=Matrix(I*one(Float64),3,3)*3
        A1_sparse=ones(Float64,3,3);
        dep0_sparse=DEP([A0_sparse,A1_sparse])
        push!(testlist,compute_types_metadata(dep0_sparse,Float64,true,[],[],bigfloat_types,true, "sparse"));


        # A bigfloat DEP
        dep_bigfloat=DEP([A0,A1]);

        push!(testlist,compute_types_metadata(dep_bigfloat,BigFloat,true,[],[],bigfloat_types,true));

        A0bc=Matrix{Complex{BigFloat}}(A0);
        A1bc=Matrix{Complex{BigFloat}}(A1);
        dep_cbigfloat=DEP([A0bc,A1bc]);


        push!(testlist,compute_types_metadata(dep_cbigfloat,Complex{BigFloat},false,[],[],bigfloat_types,false));


        # SPMF 1
        Random.seed!(0);
        B0=randn(3,3); B1=randn(3,3);
        oneop= S-> S;
        sqrop= S-> S^2;
        spmf_nep=SPMF_NEP([B0,B1],[oneop,sqrop],Ftype=Float64)
        push!(testlist,compute_types_metadata(spmf_nep,Float64,true,[],[],[],false, "#1"));

        # SPMF 2
        B0=randn(3,3); B1=randn(3,3);
        oneop= S-> 1im*S;
        sqrop= S-> S^2;
        spmf_nep2=SPMF_NEP([B0,B1],[oneop,sqrop],Ftype=ComplexF64)
        push!(testlist,compute_types_metadata(spmf_nep2,ComplexF64,false,[],[],[],false, "#2"));

        # SPMF 3
        expmop= S -> exp(S)
        B0b=Matrix{BigFloat}(B0);
        B1b=Matrix{BigFloat}(B1);
        spmf_nep3=SPMF_NEP([B0b,B1b],[oneop,expmop],Ftype=Float64)


        # Skip these since expm not supported
        bigfloats_and_float16=[BigFloat,Complex{BigFloat}];

        push!(testlist,compute_types_metadata(spmf_nep3,BigFloat,true,
                                              bigfloats_and_float16,bigfloats_and_float16,
                                              bigfloats_and_float16,false, "#3"));

        # SPMF nep 4
        oneop= S-> 1im*S;
        sqrop= S-> S^2;
        spmf_nep4=SPMF_NEP([A0_sparse,A1_sparse],[oneop,sqrop],Ftype=ComplexF64)

        push!(testlist,compute_types_metadata(spmf_nep4,ComplexF64,false,[],[],[],false, "#4"));


        # A complex sumnep
        sumnep=SumNEP(pep_complex,spmf_nep2);
        push!(testlist,compute_types_metadata(sumnep,ComplexF64,false,
                                              [],[],[],false));

    end


    testlist=testlist[1:end]



    for i in 1:size(testlist,1)
        test_one_nep(testlist[i],typelist)
    end


end
