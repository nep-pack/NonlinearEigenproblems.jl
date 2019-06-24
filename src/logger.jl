# Logging functionality

export Logger, PrintLogger, ErrorLogger
export push_info!, push_iteration_info!
import Base.push!
using Printf

    """
    abstract type Logger ; end

The type represents a way to log information throughout the algorithms,
as well as a storage. The most common `Logger`s are `PrintLogger` and
`ErrorLogger`.

As a method developer you want to use `push_info!` and
`push_iteration_info!`.

 See also: [`PrintLogger`](@ref) and [`ErrorLogger`](@ref).
"""
abstract type Logger ; end

   """
    function push_info!(logger, [level,] v; continues::Bool=false)

Pushes a string v to the logger. If continues=true, the next
`push_info!` (or `push_iteration_info!`) is connected to this,
e.g. line-feed will be omitted.
"""
function push_info!(logger::Logger,v::String; continues::Bool=false)
    push_info!(logger,1,v;continues=continues);
end
   """
    function push_iteration_info!(logger, [level,] iter; kwargs)

Pushes information about a specific iteration `iter`. The supported
keyword arguments include `λ`, `err` and `v`.

"""
function push_iteration_info!(logger::Logger,iter; kwargs...)
    push_iteration_info!(logger,1,iter;kwargs...);
end


   """
   struct PrintLogger <: Logger ;

When you use this logger, you will obtain printouts in stdout,
no other logging. The displaylevel parameter specified,
how much should be printed.
"""
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


    """
    struct ErrorLogger <: Logger

When you use this logger, the error of `push_iteration_info!`-calls
will be stored in `logger.errs`. It can also print to stdout,
if `displaylevel` is set to a value greater than zero.

"""
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
            err_vec=err;
            # Make sure err_vec is a vector
            if (err_vec isa Number)
                err_vec=[err_vec]
            end
            logger.errs[iter,1:size(err,1)]=err_vec;
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
