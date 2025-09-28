module mopex2
using ..HydroModels
using ..HydroModels: step_func

# Model variables
@variables P [description = "precipitation", unit = "mm/d"]
@variables T [description = "temperature", unit = "oC"]
@variables Ep [description = "potential evapotransporation", unit = "mm/d"]
@variables Ps Pr Qn Et1 Q1f Qw Et2 Q2u Qf Qs Qt Sn S1 S2 Sc1 Sc2 Qu Qf Qt Q2f
model_variables = [P, T, Ep, Ps, Pr, Qn, Et1, Q1f, Qw, Et2, Q2u, Qf, Qs, Qt, Sn, S1, S2, Sc1, Sc2, Qu, Qf, Qt, Q2f]
# Model parameters
@parameters tcrit [description = "Snowfall & snowmelt temperature", bounds = (-3, 3), unit = "oC"]
@parameters ddf [description = "Degree-day factor for snowmelt", bounds = (0, 20), unit = "mm/oC/d"]
@parameters Sb1 [description = "Maximum soil moisture storage", bounds = (1, 2000), unit = "mm"]
@parameters tw [description = "Groundwater leakage time", bounds = (0, 1), unit = "d-1"]
@parameters tu [description = "Slow flow routing response time", bounds = (0, 1), unit = "d-1"]
@parameters Se [description = "Root zone storage capacity", bounds = (1, 2000), unit = "mm"]
@parameters tc [description = "Mean residence time", bounds = (0, 1), unit = "d-1"]
model_parameters = [tcrit, ddf, Sb1, tw, tu, Se, tc]

# Soil water component
bucket_1 = @hydrobucket :bucket_1 begin
    fluxes = begin
        @hydroflux Ps ~ step_func(tcrit - T) * P
        @hydroflux Qn ~ min(Sn, ddf * max(0.0, T - tcrit))
    end

    dfluxes = begin
        @stateflux Sn ~ Ps - Qn
    end
end

bucket_2 = @hydrobucket :bucket_2 begin
    fluxes = begin
        @hydroflux Pr ~ step_func(T - tcrit) * P
        @hydroflux Et1 ~ clamp(S1 / Sb1, 0.0, 1.0) * Ep
        @hydroflux Et2 ~ clamp(S2 / Se, 0.0, 1.0) * Ep
        @hydroflux Q1f ~ step_func(S1 - Sb1) * (Pr + Qn)
        @hydroflux Q2u ~ tu * S2
        @hydroflux Qw ~ tw * S1
    end

    dfluxes = begin
        @stateflux S1 ~ Pr + Qn - Et1 - Q1f - Qw
        @stateflux S2 ~ Qw - Et2 - Q2u
    end
end

bucket_3 = @hydrobucket :bucket_4 begin
    fluxes = begin
        @hydroflux Qf ~ tc * Sc1
        @hydroflux Qu ~ tc * Sc2
        @hydroflux Qt ~ Qf + Qu
    end

    dfluxes = begin
        @stateflux Sc1 ~ Q1f - Qf
        @stateflux Sc2 ~ Q2u - Qu
    end
end

model = @hydromodel :mopex1 begin
    bucket_1
    bucket_2
    bucket_3
end

end