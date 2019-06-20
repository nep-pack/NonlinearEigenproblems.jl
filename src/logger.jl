# Logging functionality

export PrintLogger, Logger, ErrorLogger
export push_info!, push_iteration_info!
import Base.push!
using Printf

abstract type Logger ; end

function push_info!(logger::Logger,v::String; continues::Bool=false)
    push_info!(logger,1,v;continues=continues);
end
function push_iteration_info!(logger::Logger,iter; kwargs...)
    push_iteration_info!(logger,1,iter;kwargs...);
end



# PrintLogger: printouts only
struct PrintLogger <: Logger ;
    displaylevel::Int
end


function push_info!(logger::PrintLogger,level,v::String;continues::Bool=false)
    if (logger.displaylevel>=level)
        print(v);
        if (!continues)
            println();
        end
    end
end

function push_iteration_info!(logger::PrintLogger,level,iter;
                              err=Inf,λ=NaN,v=NaN,
                              continues::Bool=false)
    if (logger.displaylevel>=level)
        print("iter ",iter, " err:",err, " λ=",λ);
        if (!continues)
            println()
        end
    end
end



# PrintLogger: printouts only
struct ErrorLogger{T} <: Logger  where {T};
    errs::Matrix{T}
    printlogger::Logger
end


function ErrorLogger(nof_eigvals=100,nof_iters=100,displaylevel=1)
    s=Matrix{Float64}(undef,nof_iters,nof_eigvals)
    s[:] .= NaN;
    printlogger=PrintLogger(displaylevel)
    return ErrorLogger{eltype(s)}(s,printlogger);
end



function push_iteration_info!(logger::ErrorLogger,level,iter;
                              err=Inf,λ=NaN,v=NaN,
                              continues::Bool=false)
    if (iter<=size(logger.errs,1))
        if (size(err,1)<=size(logger.errs,2))
            logger.errs[iter,1:size(err,1)]=err;
        else
            if (logger.printlogger.displaylevel>1)
                println("Warning: Insufficient space in logger matrix");
            end
        end

    end

    push_iteration_info!(logger.printlogger,level,iter;
                              err=err,λ=λ,v=v,
                              continues=continues)

end
function push_info!(logger::ErrorLogger,level::Int,
                    v::String;continues::Bool=false)
    push_info!(logger.printlogger,level,v,continues=continues)
end
