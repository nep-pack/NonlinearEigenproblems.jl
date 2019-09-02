export init_errmeasure
export estimate_error
export Errmeasure
export BackwardErrmeasure
export ResidualErrmeasure
export DefaultErrmeasure
export ErrmeasureType


"""
    abstract type Errmeasure; end

Concrete subtypes of `Errmeasure` represent specific ways of
measuring the error of an eigenpair. NEP-solvers take a type as
input and then instantiate a `Errmeasure`. As a NEP-solver user,
you use the type as follows
```julia
julia> quasinewton(nep,errmeasure=ResidualErrmeasure)
```
User-specified ways of measuring error can be given by
creating a new subtype of `Errmeasure` and using it as a
`errmeasure` keyword. You need to specify the way to
measure the error in the method `estimate_error` and (optionally)
`init_errmeasure`.


See also: [`DefaultErrmeasure`](@ref), [`ResidualErrmeasure`](@ref), [`BackwardErrmeasure`](@ref), [`estimate_error`](@ref), [`init_errmeasure`](@ref).

"""
abstract type Errmeasure; end

"""
    ErrmeasureType = Union{Type{<:Errmeasure}, Function}


# Example
This shows how to compute a reference solution and
then use this as a reference solution. The error
in the second run will be effectively the
eigenvalue error.
```julia
julia> nep=nep_gallery("qdep0");
julia> (λref,vref)=quasinewton(nep,λ=-1,v=ones(size(nep,1)));
julia> myerrmeasure=(λ,v) -> abs(λref-λ)
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
"""
ErrmeasureType = Union{Type{<:Errmeasure}, Function}

"""
    function init_errmeasure(E::Errmeasure,nep)

This function is called in a precomputation phase for an error measure.
For user defined `Errmeasure`, you do not need to overload this
if your `Errmeasure` contains only one field which is a `NEP`.

See also: [`Errmeasure`](@ref)
"""
function init_errmeasure(E::Type{<:Errmeasure},nep::NEP)
    return E(nep);
end

"""
    function estimate_error(E::Errmeasure,λ,v)

Returns the error estimate for the eigenpair `(λ,v)`. The way to measure
the error is specified in `Errmeasure`.

See also: [`Errmeasure`](@ref)


"""
function estimate_error(errdata::Errmeasure,λ,v)
    error("Not implemented:",errdata);
end


"""
    struct ResidualErrmeasure <: Errmeasure

This `Errmeasure` species that the residual norm should be
used to measure the error.

See also: [`Errmeasure`](@ref)

"""
struct ResidualErrmeasure <: Errmeasure;
    nep::NEP
end

function estimate_error(errdata::ResidualErrmeasure, λ,v)
   return norm(compute_Mlincomb(errdata.nep,λ,v))/norm(v);
end


"""
    struct BackwardErrmeasure <: Errmeasure

This `Errmeasure` provides a way to compute the backward error.
The backward error estimate are only given for NEPs which are
subtypes of  `AbstractSPMF`.

# Example
```julia
julia> nep=nep_gallery("qdep0");
julia> (λ,v)=quasinewton(nep,λ=-1,errmeasure=BackwardErrmeasure,tol=1e-10);
```
See also: [`Errmeasure`](@ref)

"""
struct BackwardErrmeasure{X<:Real} <: Errmeasure
    nep::NEP
    coeffs::Vector{X};
end

function init_errmeasure(E::Type{BackwardErrmeasure},nep::AbstractSPMF)
    Av=get_Av(nep);
    # Note: norm(A) is a the frobenius norm in Julia
    coeffs=map(i->norm(Av[i]),1:size(Av,1));
    return BackwardErrmeasure(nep,coeffs);
end

function estimate_error(errdata::BackwardErrmeasure, λ,v)
    Av=get_Av(errdata.nep);
    fv=get_fv(errdata.nep);
    denom=mapreduce(i->errdata.coeffs[i]*abs(fv[i](λ)), +, 1:size(Av,1));
    return norm(compute_Mlincomb(errdata.nep,λ,v))/(norm(v)*denom);
end


"""
    struct DefaultErrmeasure <: Errmeasure

When you specify this `Errmeasure`, NEP-PACK tries to determine
a suitable `Errmeasure` based on the type of the `NEP`.
Note that this behavior may change in future versions.

See also: [`Errmeasure`](@ref)

"""
abstract type DefaultErrmeasure <: Errmeasure; end # Default behavior: If AbstractSPMF -> do backward error. Otherwise residual norm.

init_errmeasure(E::Type{DefaultErrmeasure},nep::AbstractSPMF)=init_errmeasure(BackwardErrmeasure,nep)
init_errmeasure(E::Type{DefaultErrmeasure},nep::NEP)=init_errmeasure(ResidualErrmeasure,nep)



struct UserDefinedErrmeasure <: Errmeasure
   nep::NEP
   errmeasure_fun::Function
end

function init_errmeasure(f::Function, nep::NEP)
    return UserDefinedErrmeasure(nep,f);
end
function estimate_error(e::UserDefinedErrmeasure, λ,v)
    return e.errmeasure_fun(λ,v)
end
