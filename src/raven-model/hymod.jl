module RavenHYMOD

using ..HydroModels
using ..HydroModelLibrary: SOILEVAP_TOPMODEL, INF_PDM, BASE_LINEAR

export build_raven_hymod, model, model_parameters, model_variables

@variables P [description = "Precipitation", unit = "mm/d"]
@variables Ep [description = "Potential evapotranspiration", unit = "mm/d"]

@variables Sm [description = "Soil moisture storage", unit = "mm"]
@variables S [description = "Slow reservoir storage", unit = "mm"]
@variables Sf1 [description = "Fast reservoir 1 storage", unit = "mm"]
@variables Sf2 [description = "Fast reservoir 2 storage", unit = "mm"]
@variables Sf3 [description = "Fast reservoir 3 storage", unit = "mm"]

@variables Ea [description = "Actual evaporation", unit = "mm/d"]
@variables infil [description = "Soil storage increment", unit = "mm/d"]
@variables Pe [description = "Effective precipitation", unit = "mm/d"]
@variables Pf [description = "Fast flow input", unit = "mm/d"]
@variables Ps [description = "Slow flow input", unit = "mm/d"]
@variables Qf1 [description = "Fast reservoir 1 outflow", unit = "mm/d"]
@variables Qf2 [description = "Fast reservoir 2 outflow", unit = "mm/d"]
@variables Qf3 [description = "Fast reservoir 3 outflow", unit = "mm/d"]
@variables Qs [description = "Slow reservoir outflow", unit = "mm/d"]
@variables Qt [description = "Total runoff", unit = "mm/d"]

model_variables = [P, Ep, Sm, S, Sf1, Sf2, Sf3, Ea, infil, Pe, Pf, Ps, Qf1, Qf2, Qf3, Qs, Qt]

@parameters Smax [description = "Maximum soil moisture capacity", bounds = (10.0, 2000.0), unit = "mm"]
@parameters b [description = "PDM shape parameter", bounds = (0.0, 10.0), unit = "-"]
@parameters a [description = "Fraction of effective precipitation to fast reservoirs", bounds = (0.0, 1.0), unit = "-"]
@parameters kf [description = "Fast reservoir coefficient", bounds = (0.0, 1.0), unit = "d-1"]
@parameters ks [description = "Slow reservoir coefficient", bounds = (0.0, 1.0), unit = "d-1"]

model_parameters = [Smax, b, a, kf, ks]

soil_bucket = @hydrobucket :raven_hymod_soil begin
    fluxes = begin
        SOILEVAP_TOPMODEL(;
            soil_evaporation=Ea,
            waterstorage=Sm,
            potential_evaporation=Ep,
            storage_tension=Smax,
            pet_corr=1.0,
            flux_name=:raven_hymod_ea,
        )
        INF_PDM(;
            infiltration=infil,
            prcpeff=P,
            waterstorage=Sm,
            max_waterstorage=Smax,
            pdm_b=b,
            flux_name=:raven_hymod_infil,
        )
        @hydroflux Pe ~ max(P - infil, 0.0)
        @hydroflux Pf ~ a * Pe
        @hydroflux Ps ~ (1 - a) * Pe
    end
    dfluxes = begin
        @stateflux Sm ~ infil - Ea
    end
end

fast_bucket = @hydrobucket :raven_hymod_fast begin
    fluxes = begin
        BASE_LINEAR(; baseflow=Qf1, waterstorage=Sf1, baseflow_coeff=kf, flux_name=:raven_hymod_qf1)
        BASE_LINEAR(; baseflow=Qf2, waterstorage=Sf2, baseflow_coeff=kf, flux_name=:raven_hymod_qf2)
        BASE_LINEAR(; baseflow=Qf3, waterstorage=Sf3, baseflow_coeff=kf, flux_name=:raven_hymod_qf3)
    end
    dfluxes = begin
        @stateflux Sf1 ~ Pf - Qf1
        @stateflux Sf2 ~ Qf1 - Qf2
        @stateflux Sf3 ~ Qf2 - Qf3
    end
end

slow_bucket = @hydrobucket :raven_hymod_slow begin
    fluxes = begin
        BASE_LINEAR(; baseflow=Qs, waterstorage=S, baseflow_coeff=ks, flux_name=:raven_hymod_qs)
        @hydroflux Qt ~ Qs + Qf3
    end
    dfluxes = begin
        @stateflux S ~ Ps - Qs
    end
end

model = @hydromodel :raven_hymod begin
    soil_bucket
    fast_bucket
    slow_bucket
end

build_raven_hymod() = model

end


