
export estimate_error
export Errmeasure
export StandardSPMFErrmeasure
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


See also: [`ErrmeasureType`](@ref), [`DefaultErrmeasure`](@ref), [`ResidualErrmeasure`](@ref), [`StandardSPMFErrmeasure`](@ref), [`estimate_error`](@ref), [`EigvalReferenceErrmeasure`](@ref).

"""
abstract type Errmeasure; end


"""
    ErrmeasureType = Union{Type{<:Errmeasure}, Function}

The `ErrmeasureType` represents (essentially) what you can insert in the
`errmeasure` keyword argument for most NEP-solvers. It can be
a function or an  [`Errmeasure`](@ref) object.
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
julia> myerrmeasure=(λ,v) -> norm(vref/vref[1]-v/v[1]);
julia> (λ,v)=quasinewton(nep,errmeasure=myerrmeasure,λ=-1.0 ,logger=1,tol=5e-13,v=ones(size(nep,1)));
Precomputing linsolver
iter 1 err:46.40296482739195 λ=-1.0 + 0.0im
iter 2 err:2.1592671533657883 λ=-0.7063330111559607 + 0.0im
iter 3 err:0.17079231439255405 λ=-0.8919579082730457 + 0.0im
iter 4 err:0.1633846991066227 λ=-1.0097584042560848 + 0.0im
iter 5 err:0.003434042059262583 λ=-1.0023823873044 + 0.0im
iter 6 err:0.0003182517281689052 λ=-1.0024660870524031 + 0.0im
iter 7 err:2.0105257231740345e-5 λ=-1.0024677891861997 + 0.0im
iter 8 err:1.618661190265619e-6 λ=-1.0024669496893164 + 0.0im
iter 9 err:1.233489068442819e-7 λ=-1.0024669918249076 + 0.0im
iter 10 err:9.44707811957546e-9 λ=-1.0024669883439166 + 0.0im
iter 11 err:7.867601351698812e-10 λ=-1.0024669886059947 + 0.0im
iter 12 err:0.0 λ=-1.002466988585764 + 0.0im

```
The eigenvalue error can be measured with the
[`EigvalReferenceErrmeasure`](@ref).

See also: [`Errmeasure`](@ref)
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
            errm=StandardSPMFErrmeasure(nep);
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
    struct StandardSPMFErrmeasure <: Errmeasure
    function StandardSPMFErrmeasure(nep::AbstractSPMF)

This `Errmeasure` provides a way to compute the backward error.
The backward error estimate are only given for NEPs which are
subtypes of  `AbstractSPMF`. We use the Frobenius norm as the
matrix norm, since it is much cheaper to compute than the spectral
norm.

# Example
```julia
julia> nep=nep_gallery("qdep0");
julia> (λ,v)=quasinewton(nep,λ=-1,v=ones(size(nep,1)),errmeasure=StandardSPMFErrmeasure(nep),tol=1e-10,logger=1);
Precomputing linsolver
iter 1 err:0.022010375110869937 λ=-1.0 + 0.0im
iter 2 err:0.002515422247048546 λ=-0.7063330111559607 + 0.0im
iter 3 err:0.000892354247568813 λ=-0.8919579082730457 + 0.0im
iter 4 err:5.445678793151584e-5 λ=-1.0097584042560848 + 0.0im
iter 5 err:6.649967517409105e-7 λ=-1.0023823873044 + 0.0im
iter 6 err:1.0557281809769784e-8 λ=-1.0024660870524031 + 0.0im
iter 7 err:6.420125566431444e-9 λ=-1.0024677891861997 + 0.0im
iter 8 err:3.181093707909799e-10 λ=-1.0024669496893164 + 0.0im
iter 9 err:2.6368050026394416e-11 λ=-1.0024669918249076 + 0.0im

```
See also: [`Errmeasure`](@ref)

"""
struct StandardSPMFErrmeasure{X<:Real} <: Errmeasure
    nep::AbstractSPMF
    coeffs::Vector{X};
    function StandardSPMFErrmeasure(nep::AbstractSPMF)
        Av=get_Av(nep);
        # Note: norm(A) is a the frobenius norm in Julia
        coeffs=map(i->norm(Av[i]),1:size(Av,1));
        return new{eltype(coeffs)}(nep,coeffs);
    end
end


function estimate_error(errm::StandardSPMFErrmeasure, λ,v)
    Av=get_Av(errm.nep);
    fv=get_fv(errm.nep);
    denom=mapreduce(i->errm.coeffs[i]*abs(fv[i](λ)), +, 1:size(Av,1));
    return norm(compute_Mlincomb(errm.nep,λ,v))/(norm(v)*denom);
end



"""
    struct EigvalReferenceErrmeasure{X<:Number} <: Errmeasure
    function EigvalReferenceErrmeasure(nep,λref)

Use the difference between a precomputed λ-value (reference solution)
and the eigenvalue estimate
as the error measure.


# Example
```julia
julia> using LinearAlgebra
julia> nep=nep_gallery("qdep0");
julia> (λref,vref)=quasinewton(nep,λ=-1,v=ones(size(nep,1)));
julia> (λ,v)=quasinewton(nep,errmeasure=EigvalReferenceErrmeasure(nep,λref),λ=-1.0 ,logger=1,tol=5e-13,v=ones(size(nep,1)));
Precomputing linsolver
iter 1 err:0.002466988585763996 λ=-1.0 + 0.0im
iter 2 err:0.2961339774298033 λ=-0.7063330111559607 + 0.0im
iter 3 err:0.11050908031271833 λ=-0.8919579082730457 + 0.0im
iter 4 err:0.007291415670320767 λ=-1.0097584042560848 + 0.0im
iter 5 err:8.460128136400513e-5 λ=-1.0023823873044 + 0.0im
iter 6 err:9.015333608530796e-7 λ=-1.0024660870524031 + 0.0im
iter 7 err:8.006004357241636e-7 λ=-1.0024677891861997 + 0.0im
iter 8 err:3.889644761834177e-8 λ=-1.0024669496893164 + 0.0im
iter 9 err:3.2391436199930013e-9 λ=-1.0024669918249076 + 0.0im
iter 10 err:2.418474309706653e-10 λ=-1.0024669883439166 + 0.0im
iter 11 err:2.0230705999324528e-11 λ=-1.0024669886059947 + 0.0im
iter 12 err:0.0 λ=-1.002466988585764 + 0.0im
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
