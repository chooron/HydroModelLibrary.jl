module ihacres
using ..HydroModels
using ..HydroModels: step_func

# Model variables
@variables CMD P Ea U Ep Qt Uq Us Xs Xq q
@variables t
model_variables = [CMD, P, Ea, U, Ep, Qt, q ,Uq, Us, Xs, Xq]

# Model parameters
@parameters lp [description = "Wilting point", bounds = (1, 2000), unit = "mm"]
@parameters d [description = "Deficit threshold for flow from rain", bounds = (1, 2000), unit = "mm"]
@parameters p [description = "Deficit non-linearity", bounds = (0, 10), unit = "-"]
@parameters alpha [description = "Fraction flow to quick routing", bounds = (0, 1), unit = "-"]
@parameters Tq [description = " Unit Hydrograph time base", bounds = (1, 700), unit = "d"]
@parameters Ts [description = "Unit Hydrograph time base", bounds = (1, 700), unit = "d"]
@parameters Td [description = "Unit Hydrograph delay", bounds = (0, 119), unit = "d"]

@parameters tau_q [description = "Fast flow routing delay", bounds = (1, 10), unit = "d"]
@parameters tau_s [description = "Slow flow routing delay", bounds = (1, 10), unit = "d"]
@parameters tau_d [description = "Delay", bounds = (1, 10), unit = "d"]

model_parameters = [lp, d, p, alpha, Tq, Ts, Td, tau_q, tau_s, tau_d]
bucket1 = @hydrobucket :bucket1 begin
    fluxes = begin
        @hydroflux Ea ~ Ep * min(1, exp(2 * (1 - CMD / lp)))
        @hydroflux U ~ P * (1 - clamp(CMD / d, 0, 1)^p)
        @hydroflux Uq ~ alpha * U
        @hydroflux Us ~ (1 - alpha) * U
    end
    dfluxes = begin
        @stateflux CMD ~ -P + Ea + U
    end
end

uh_1 = @unithydro :uh_1 begin
    uh_vars = Uq => Xq
    uh_func = begin
        tau_q => -exp(-t)
    end
end

uh_2 = @unithydro :uh_2 begin
    uh_vars = Us => Xs
    uh_func = begin
        tau_s => -exp(-t)
    end
end

q_flux = @hydroflux q ~ Xq + Xs

uh_3 = HydroModels.UnitHydrograph(
    [q], [Qt], [tau_d],
    (t, p) -> begin
        lag = p.params.tau_d
        if t == floor(lag)
            return 1 - lag + floor(lag)
        elseif t == floor(lag) + 1
            return lag - floor(lag)
        else
            return 0
        end
    end, (p) -> p.params.tau_d, name=:uh_3
)

model = @hydromodel :ihacres begin
    bucket1
    uh_1
    uh_2
    q_flux
    uh_3
end

end