function prettypercent(v, pad=false)
    if isnan(v)
        str = pad ? "--- " : "---"
    else
        i = round(Int, v)
        str = (i == 0 ? "0" : i > 999 ? ">999" : i < -999 ? "<-999" : @sprintf("%+d", i)) * "%"
    end
    return pad ? lpad(str, 5, " ") : str
end

const TIME_UNITS = [(1e3, 1, "ns"), (1e3, 1e3, "μs"), (1e3, 1e6, "ms"), (Inf, 1e9, "s")]

prettytime(v, pad=false) = pretty(TIME_UNITS, v, pad)

const MEMORY_UNITS = [(1e3, 1, "B"), (1e3, 1024, "KB"), (1e3, 1024^2, "MB"), (Inf, 1024^3, "GB")]

prettymemory(v, pad=false) = pretty(MEMORY_UNITS, v, pad)

function pretty(units, v, pad)
    _, scale, unit = units[findfirst(x -> round(Int, v/x[2]) < x[1], units)]
    str = three_significant_digits(v / scale) * " $unit"
    return pad ? lpad(str, 7, " ") : str
end

three_significant_digits(v) =
    round(v) >= 100 ? @sprintf("%.0f", v) :
    round(v * 10) >= 100 ? @sprintf("%.1f", v) : @sprintf("%.2f", v)
