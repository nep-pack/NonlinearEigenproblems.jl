# Logging functionality

export PrintLogger, Logger
import Base.push!
using Printf

abstract type Logger ; end

function push_info!(logger::Logger,v::String,continues::Bool=false)
    push_info!(logger,1,v,continues);
end


struct PrintLogger <: Logger ;
    displaylevel::Int
end


function push_info!(logger::PrintLogger,level,v::String,continues::Bool=false)
    if (logger.displaylevel>=level)
        print(v);
        if (!continues)
            println();
        end
    end
end

function push_iteration_info!(logger::PrintLogger,iter,continues::Bool=false;
                              err=Inf,λ=NaN,v=NaN)
    if (logger.displaylevel>=1)
        print("iter ",iter, " err:",err, " λ=",λ);
        if (!continues)
            println()
        end
    end
end
