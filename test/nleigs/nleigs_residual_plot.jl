using NonlinearEigenproblems.RKHelper
using Plots

function nleigs_residual_plot(title, solution_info, Σ; plot_attributes...)
    Lam = solution_info.Lam
    Res = solution_info.Res

    n = size(Lam, 1)
    LamC = zero(Lam)
    ResC = zero(Res)

    # first iteration
    LamC[1,1] = Lam[1,1]
    ResC[1,1] = Res[1,1]

    # next iterations
    for i = 2:n
        sLam = Lam[1:i,i]
        sRes = Res[1:i,i]
        for j = 1:i-1
            _,ii = findmin(abs.(sLam .- LamC[j,i-1]))
            LamC[j,i] = sLam[ii]
            ResC[j,i] = sRes[ii]
            deleteat!(sLam, ii)
            deleteat!(sRes, ii)
        end
        LamC[i,i] = sLam[1]
        ResC[i,i] = sRes[1]
    end

    z = map(p -> inpolygon(real(p), imag(p), real(Σ), imag(Σ)), LamC[:,end])

    ResIn = ResC[z,:]
    ResOut = ResC[.!z,:]
    LamIn = LamC[z,:]
    LamOut = LamC[.!z,:]

    p = plot(title=title, xlabel="iteration", ylabel="residual", leg=false; plot_attributes...)

    for i = 1:n
        color, style = isempty(Σ) || z[i] ? (:red, :solid) : (:black, :dot)
        p = plot!(p, i:n, ResC[i,i:n], yscale=:log10, linecolor=color, linestyle=style)
    end

    display(p)
end
