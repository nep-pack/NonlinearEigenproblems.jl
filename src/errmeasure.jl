export init_errmeasure
export estimate_error
export Errmeasure
export BackwardErrmeasure
export ResidualErrmeasure
export DefaultErrmeasure



abstract type Errmeasure; end


function init_errmeasure(E::Type{<:Errmeasure},nep::NEP)
    return E(nep);
end

function estimate_error(errdata::Errmeasure,λ,v)
    error("Not implemented:",errdata);
end

struct ResidualErrmeasure <: Errmeasure;
    nep::NEP
end

function estimate_error(errdata::ResidualErrmeasure, λ,v)
   return norm(compute_Mlincomb(errdata.nep,λ,v))/norm(v);
end


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

# Default behavior: If AbstractSPMF -> do backward error. Otherwise residual norm.
abstract type DefaultErrmeasure <: Errmeasure; end
init_errmeasure(E::Type{DefaultErrmeasure},nep::AbstractSPMF)=init_errmeasure(BackwardErrmeasure,nep)
init_errmeasure(E::Type{DefaultErrmeasure},nep::NEP)=init_errmeasure(ResidualErrmeasure,nep)
