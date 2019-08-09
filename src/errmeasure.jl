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

# Example
This shows how to compute a reference solution and
then use this as a reference solution. The error
in the second run will be effectively the
eigenvalue error. 
```julia
julia> nep=nep_gallery("qdep0");
julia> (λref,vref)=quasinewton(nep,λ=-1);
julia> struct EigvalError <: Errmeasure; nep::NEP; end
julia> function NonlinearEigenproblems.estimate_error(E::EigvalError,λ,v)
return abs(λref-λ);
end
julia> (λ,v)=quasinewton(nep,errmeasure=EigvalError,λ=-1.0 ,logger=1,tol=5e-13)
Precomputing linsolver
Iteration:  1 errmeasure:2.466988587467300320e-03, λ=-1.0 + 0.0im
Iteration:  2 errmeasure:4.625160667012763183e-01, λ=-0.539950921886191 + 0.0im
Iteration:  3 errmeasure:2.799848726755167494e-01, λ=-1.282451861262984 + 0.0im
Iteration:  4 errmeasure:3.422925625951256379e-02, λ=-1.0366962448469799 + 0.0im
Iteration:  5 errmeasure:5.530128437585268841e-04, λ=-1.0019139757437088 + 0.0im
Iteration:  6 errmeasure:1.159388768512403800e-04, λ=-1.002351049710616 + 0.0im
Iteration:  7 errmeasure:2.658434455016234210e-06, λ=-1.0024643301530123 + 0.0im
Iteration:  8 errmeasure:1.726871190488310503e-07, λ=-1.0024668159003483 + 0.0im
Iteration:  9 errmeasure:4.819693533164581822e-09, λ=-1.0024669934071608 + 0.0im
Iteration: 10 errmeasure:5.234268574128009277e-10, λ=-1.0024669880640404 + 0.0im
Iteration: 11 errmeasure:3.762568034915148019e-11, λ=-1.002466988625093 + 0.0im
Iteration: 12 errmeasure:3.205657961302676995e-12, λ=-1.0024669885842616 + 0.0im
Iteration: 13 errmeasure:3.352873534367972752e-14, λ=-1.0024669885874338 + 0.0im
```
Note that this can also be achieved by providing a function handle:
```julia
julia> myerrmeasure= (λ,v) -> abs(λref-λ);
julia> (λ,v)=quasinewton(nep,errmeasure=myerrmeasure,λ=-1.0 ,logger=1,tol=5e-13)
...
Iteration: 12 errmeasure:3.205657961302676995e-12, λ=-1.0024669885842616 + 0.0im
Iteration: 13 errmeasure:3.352873534367972752e-14, λ=-1.0024669885874338 + 0.0im
```

See also: [`DefaultErrmeasure`](@ref), [`ResidualErrmeasure`](@ref), [`BackwardErrmeasure`](@ref), [`estimate_error`](@ref), [`init_errmeasure`](@ref).

"""
abstract type Errmeasure; end

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
