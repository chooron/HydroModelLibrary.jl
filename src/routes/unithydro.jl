@variables t

function uh_1_half(; input, output, lag, name=:uh_1_half)
    return HydroModels.@unithydro :maxbas_uh begin
        uh_vars = input => output
        uh_func = begin
            lag => (t / lag)^2.5
        end
    end
end

function uh_2_full(; input, output, lag, name=:uh_1_full)
    return HydroModels.@unithydro name begin
        uh_vars = input => output
        uh_func = begin
            2 * (lag) => 1 - 0.5 * (2 - t / lag)^(2.5)
            (lag) => 0.5 * (t / lag)^(2.5)
        end
    end
end

function uh_3_half(; input, output, lag, name=:uh_3_half)
    return HydroModels.@unithydro name begin
        uh_vars = input => output
        uh_func = begin
            lag => 1 / (0.5 * lag^2) * (0.5 * min(t, lag)^2 - 0.5 * (t - 1)^2)
        end
    end
end

function uh_4_full(; input, output, lag, name=:uh_4_full)
    return HydroModels.@unithydro name begin
        uh_vars = input => output
        uh_func = begin
            lag => 1 / (0.5 * lag^2) * (0.5 * min(t, lag)^2 - 0.5 * (t - 1)^2)
        end
    end
end

function uh_5_half(; input, output, lag, name=:uh_5_half)
    return HydroModels.@unithydro name begin
        uh_vars = input => output
        uh_func = begin
            lag => -exp(-t)
        end
    end
end

function uh_7_uniform(; input, output, lag, name=:uh_7_uniform)
    uh_func(t, pas) = begin
        lag = pas.params.lag
        if t == ceil(lag)
            return (lag % (t - 1)) / lag
        else
            return 1 / lag
        end
    end

    max_lag_func(pas) = begin
        return ceil(pas.params.lag)
    end

    return HydroModels.UnitHydrograph(
        [input], [output], [lag], uh_func, max_lag_func, name=name
    )
end

function uh_8_delay(; input, output, lag, name=:uh_8_delay)
    uh_func(t, pas) = begin
        lag = pas.params.lag
        if t == floor(lag) + 1
            return 1 - lag + floor(lag)
        elseif t == floor(lag) + 2
            return lag - floor(lag)
        else
            return 0
        end
    end

    max_lag_func(pas) = begin
        return pas.params.lag
    end

    return HydroModels.UnitHydrograph(
        [input], [output], [lag], uh_func, max_lag_func, name=name
    )
end

function gamma_uh(; input, output, alpha, beta, name=:gamma_uh)
    return HydroModels.@unithydro name begin
        uh_func = begin
            10 => (1 / t) * (beta * t)^alpha / gamma(alpha) * exp(-beta * t)
        end
        uh_vars = input => output
    end
end