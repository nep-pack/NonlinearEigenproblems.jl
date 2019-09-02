
export estimate_error
export Errmeasure
export BackwardErrmeasure
export ResidualErrmeasure
export DefaultErrmeasure
export EigvalReferenceErrmeasure
export ErrmeasureType


"""
    abstract type Errmeasure; end

Concrete subtypes of `Errmeasure` represent specific ways of
measuring the error of an eigenpair. NEP-solvers take such
an object as input. As a NEP-solver user,
you use the type as follows
```julia
julia> quasinewton(nep,errmeasure=ResidualErrmeasure(nep))
```
User-specified ways of measuring error can be given by
creating a new subtype of `Errmeasure` and using it as a
`errmeasure` keyword. You need to specify the way to
measure the error in the method `estimate_error`.

Note that in practice a `Function` can essentially
be used instead of a `Errmeasure`-object, which is a simple
way to have user-defined error measures.
See [`ErrmeasureType`](@ref).


See also: [`ErrmeasureType`](@ref), [`DefaultErrmeasure`](@ref), [`ResidualErrmeasure`](@ref), [`BackwardErrmeasure`](@ref), [`estimate_error`](@ref), [`EigvalReferenceErrmeasure`](@ref).

"""
abstract type Errmeasure; end


"""
    ErrmeasureType = Union{Type{<:Errmeasure}, Function}

The `ErrmeasureType` represents (essentially) what you can insert in the
`errmeasure` keyword argument for most NEP-solvers.
If it is a `Function` this function will be used to obtain
error estimate.

# Example
This shows how to compute a reference solution and
then use this as a reference solution. The error
in the second run will be effectively the
eigenvector error (appropriately normalized).
```julia
julia> using LinearAlgebra
julia> nep=nep_gallery("qdep0");
julia> (λref,vref)=quasinewton(nep,λ=-1,v=ones(size(nep,1)));
julia> myerrmeasure=(λ,v) -> norm(vref/vref[1]-v/v[1])
julia> (λ,v)=quasinewton(nep,errmeasure=myerrmeasure,λ=-1.0 ,logger=1,tol=5e-13,v=ones(size(nep,1)))
Precomputing linsolver
iter 1 err:0.0024669885857651064 λ=-1.0 + 0.0im
iter 2 err:0.2961339774298044 λ=-0.7063330111559607 + 0.0im
iter 3 err:0.11050908031267426 λ=-0.8919579082730908 + 0.0im
iter 4 err:0.007291415670313883 λ=-1.009758404256079 + 0.0im
iter 5 err:8.460128136422718e-5 λ=-1.0023823873044009 + 0.0im
iter 6 err:9.01533362851481e-7 λ=-1.0024660870524023 + 0.0im
iter 7 err:8.006004341698514e-7 λ=-1.0024677891861993 + 0.0im
iter 8 err:3.889644784038637e-8 λ=-1.0024669496893173 + 0.0im
iter 9 err:3.2391431759037914e-9 λ=-1.0024669918249083 + 0.0im
iter 10 err:2.418489852828998e-10 λ=-1.0024669883439161 + 0.0im
iter 11 err:2.0229151687090052e-11 λ=-1.0024669886059943 + 0.0im
iter 12 err:0.0 λ=-1.002466988585765 + 0.0im
```
The eigenvalue error can be measured with the
[`EigvalReferenceErrmeasure`](@ref).
"""
ErrmeasureType = Union{Errmeasure, Function}

"""
    struct DefaultErrmeasure <: Errmeasure
    function DefaultErrmeasure(nep::NEP)

When you specify this `Errmeasure`, NEP-PACK tries to determine
a suitable `Errmeasure` based on the type of the `NEP`.
Note that this behavior may change in future versions.

See also: [`Errmeasure`](@ref)

"""
struct DefaultErrmeasure{X<:Errmeasure}<:Errmeasure;
    errm::X
    function DefaultErrmeasure(nep::NEP)
        if (nep isa AbstractSPMF)
            errm=BackwardErrmeasure(nep);
        else
            errm=ResidualErrmeasure(nep);
        end
        return new{typeof(errm)}(errm);
    end
end


"""
    struct ResidualErrmeasure <: Errmeasure
    function ResidualErrmeasure(nep::NEP)

This `Errmeasure` species that the residual norm should be
used to measure the error.

See also: [`Errmeasure`](@ref)

"""
struct ResidualErrmeasure <: Errmeasure;
    nep::NEP
end


"""
    function estimate_error(E::ErrmeasureType,λ,v)

Returns the error estimate for the eigenpair `(λ,v)`. The way to measure
the error is specified in `E`, which can be an `Errmeasure` or a `Function`.

See also: [`Errmeasure`](@ref), [`ErrmeasureType`](@ref)

"""
function estimate_error(errm::ResidualErrmeasure, λ,v)
   return norm(compute_Mlincomb(errm.nep,λ,v))/norm(v);
end


function estimate_error(errm::DefaultErrmeasure, λ,v)
   return estimate_error(errm.errm,λ,v);
end




function estimate_error(errm::Function, λ,v)
   return errm(λ,v);
end


"""
    struct BackwardErrmeasure <: Errmeasure
    function BackwardErrmeasure(nep::NEP)

This `Errmeasure` provides a way to compute the backward error.
The backward error estimate are only given for NEPs which are
subtypes of  `AbstractSPMF`. We use the Frobenius norm as the
matrix norm, since it is much cheaper to compute than the spectral
norm.

# Example
```julia
julia> nep=nep_gallery("qdep0");
julia> (λ,v)=quasinewton(nep,λ=-1,errmeasure=BackwardErrmeasure(nep),tol=1e-10);
```
See also: [`Errmeasure`](@ref)

"""
struct BackwardErrmeasure{X<:Real} <: Errmeasure
    nep::NEP
    coeffs::Vector{X};
    function BackwardErrmeasure(nep::NEP)
        Av=get_Av(nep);
        # Note: norm(A) is a the frobenius norm in Julia
        coeffs=map(i->norm(Av[i]),1:size(Av,1));
        return new{eltype(coeffs)}(nep,coeffs);
    end
end


function estimate_error(errm::BackwardErrmeasure, λ,v)
    Av=get_Av(errm.nep);
    fv=get_fv(errm.nep);
    denom=mapreduce(i->errm.coeffs[i]*abs(fv[i](λ)), +, 1:size(Av,1));
    return norm(compute_Mlincomb(errm.nep,λ,v))/(norm(v)*denom);
end



"""
    struct EigvalReferenceErrmeasure{X<:Number} <: Errmeasure
    function EigvalReferenceErrmeasure(λref,nep)

Use the difference between a precomputed λ-value (reference solution)
and the eigenvalue estimate
as the error measure.


# Example
```julia
julia> using LinearAlgebra
julia> nep=nep_gallery("qdep0");
julia> (λref,vref)=quasinewton(nep,λ=-1,v=ones(size(nep,1)));
julia> (λ,v)=quasinewton(nep,errmeasure=EigvalReferenceErrmeasure(nep,λref),λ=-1.0 ,logger=1,tol=5e-13,v=ones(size(nep,1)))
Precomputing linsolver
iter 1 err:0.0024669885857651064 λ=-1.0 + 0.0im
iter 2 err:0.2961339774298044 λ=-0.7063330111559607 + 0.0im
iter 3 err:0.11050908031267426 λ=-0.8919579082730908 + 0.0im
iter 4 err:0.007291415670313883 λ=-1.009758404256079 + 0.0im
iter 5 err:8.460128136422718e-5 λ=-1.0023823873044009 + 0.0im
iter 6 err:9.01533362851481e-7 λ=-1.0024660870524023 + 0.0im
iter 7 err:8.006004341698514e-7 λ=-1.0024677891861993 + 0.0im
iter 8 err:3.889644784038637e-8 λ=-1.0024669496893173 + 0.0im
iter 9 err:3.2391431759037914e-9 λ=-1.0024669918249083 + 0.0im
iter 10 err:2.418489852828998e-10 λ=-1.0024669883439161 + 0.0im
iter 11 err:2.0229151687090052e-11 λ=-1.0024669886059943 + 0.0im
iter 12 err:0.0 λ=-1.002466988585765 + 0.0im
```

See also: [`Errmeasure`](@ref)

"""
struct EigvalReferenceErrmeasure{X<:Number} <: Errmeasure
    nep::NEP
    λref::X
    function EigvalReferenceErrmeasure(nep,λref)
        return new{typeof(λref)}(nep,λref);
    end
end


function estimate_error(errm::EigvalReferenceErrmeasure, λ,v)
    return abs(errm.λref-λ)
end
