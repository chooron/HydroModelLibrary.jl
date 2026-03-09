module RavenHBVLight

using ..HydroModels
using ..HydroModelLibrary: RAINSNOW_HBV,
                           SNOBAL_HBV,
                           SOILEVAP_HBV,
                           CAPRISE_HBV,
                           PERC_CONSTANT,
                           BASE_LINEAR,
                           BASE_POWER_LAW

export build_raven_hbv_light, model, model_parameters, model_variables

@variables P [description = "Precipitation", unit = "mm/d"]
@variables T [description = "Air temperature", unit = "degC"]
@variables Ep [description = "Potential evapotranspiration", unit = "mm/d"]

@variables SP [description = "Snow water equivalent", unit = "mm"]
@variables WC [description = "Liquid water stored in snowpack", unit = "mm"]
@variables SM [description = "Soil moisture storage", unit = "mm"]
@variables UZ [description = "Upper groundwater storage", unit = "mm"]
@variables LZ [description = "Lower groundwater storage", unit = "mm"]

@variables sf [description = "Snowfall", unit = "mm/d"]
@variables rf [description = "Rainfall", unit = "mm/d"]
@variables ma [description = "Snowmelt factor", unit = "mm/(d*degC)"]
@variables potmelt [description = "Potential melt", unit = "mm/d"]
@variables refreeze [description = "Refreezing", unit = "mm/d"]
@variables snow_release [description = "Liquid release from snowpack", unit = "mm/d"]
@variables melt [description = "Snowmelt", unit = "mm/d"]
@variables cmelt [description = "Cumulative melt index", unit = "mm"]
@variables evap [description = "Actual soil evaporation", unit = "mm/d"]
@variables capillary [description = "Capillary rise from upper zone", unit = "mm/d"]
@variables recharge [description = "Recharge from soil to upper zone", unit = "mm/d"]
@variables perc [description = "Percolation to lower zone", unit = "mm/d"]
@variables q0 [description = "Upper zone outflow", unit = "mm/d"]
@variables q1 [description = "Lower zone outflow", unit = "mm/d"]
@variables Qt [description = "Total runoff", unit = "mm/d"]

model_variables = [
    P, T, Ep,
    SP, WC, SM, UZ, LZ,
    sf, rf, ma, potmelt, refreeze, snow_release, melt, cmelt,
    evap, capillary, recharge, perc, q0, q1, Qt,
]

@parameters TT [description = "Rain-snow threshold temperature", bounds = (-3.0, 5.0), unit = "degC"]
@parameters TTI [description = "Rain-snow transition interval", bounds = (0.0, 17.0), unit = "degC"]
@parameters TTM [description = "Snowmelt threshold temperature", bounds = (-3.0, 3.0), unit = "degC"]
@parameters CFR [description = "Refreezing coefficient multiplier", bounds = (0.0, 1.0), unit = "-"]
@parameters CFMAX [description = "Degree-day factor", bounds = (0.0, 20.0), unit = "mm/(d*degC)"]
@parameters WHC [description = "Snow liquid water holding fraction", bounds = (0.0, 1.0), unit = "-"]
@parameters CFLUX [description = "Capillary rise coefficient", bounds = (0.0, 4.0), unit = "mm/d"]
@parameters FC [description = "Maximum soil moisture storage", bounds = (1.0, 2000.0), unit = "mm"]
@parameters LP [description = "Soil evaporation threshold as fraction of FC", bounds = (0.05, 0.95), unit = "-"]
@parameters BETA [description = "Soil recharge exponent", bounds = (0.0, 10.0), unit = "-"]
@parameters K0 [description = "Upper-zone outflow coefficient", bounds = (0.0, 1.0), unit = "d-1"]
@parameters ALPHA [description = "Upper-zone outflow nonlinearity", bounds = (0.0, 4.0), unit = "-"]
@parameters PPERC [description = "Maximum percolation rate", bounds = (0.0, 20.0), unit = "mm/d"]
@parameters K1 [description = "Lower-zone outflow coefficient", bounds = (0.0, 1.0), unit = "d-1"]

model_parameters = [TT, TTI, TTM, CFR, CFMAX, WHC, CFLUX, FC, LP, BETA, K0, ALPHA, PPERC, K1]

snowfall_flux, rainfall_flux = RAINSNOW_HBV(;
    snowfall=sf,
    rainfall=rf,
    precipitation=P,
    temp=T,
    tt=TT,
    tti=TTI,
    snow_flux_name=:raven_hbv_light_snowfall,
    rain_flux_name=:raven_hbv_light_rainfall,
)

snow_bucket = SNOBAL_HBV(;
    snowstorage=SP,
    liquidstorage=WC,
    temp=T,
    potential_melt=potmelt,
    cumulmelt=cmelt,
    refreeze=refreeze,
    overflow=snow_release,
    snowmelt=melt,
    snowfall=sf,
    rainfall=rf,
    ma=ma,
    ka=CFR * CFMAX,
    mamax=CFMAX,
    mamin=CFMAX,
    malpha=0.0,
    tbm=TTM,
    tbf=TTM,
    const_swi=WHC,
    name=:raven_hbv_light_snow,
)

soil_bucket = @hydrobucket :raven_hbv_light_soil begin
    fluxes = begin
        SOILEVAP_HBV(;
            soil_evaporation=evap,
            snow_depth=SP,
            waterstorage=SM,
            potential_evaporation=Ep,
            storage_tension=FC * LP,
            pet_corr=1.0,
            flux_name=:raven_hbv_light_evap,
        )
        CAPRISE_HBV(;
            capillary_rise=capillary,
            lower_storage=UZ,
            soil_storage=SM,
            max_soil_storage=FC,
            capillary_coeff=CFLUX,
            flux_name=:raven_hbv_light_capillary,
        )
        @hydroflux recharge ~ snow_release * clamp(SM / FC, 0.0, 1.0)^BETA
        BASE_POWER_LAW(;
            baseflow=q0,
            waterstorage=UZ,
            baseflow_coeff=K0,
            baseflow_n=1 + ALPHA,
            flux_name=:raven_hbv_light_q0,
        )
        PERC_CONSTANT(;
            percolation=perc,
            waterstorage=UZ,
            max_perc_rate=PPERC,
            flux_name=:raven_hbv_light_perc,
        )
        BASE_LINEAR(;
            baseflow=q1,
            waterstorage=LZ,
            baseflow_coeff=K1,
            flux_name=:raven_hbv_light_q1,
        )
        @hydroflux Qt ~ q0 + q1
    end
    dfluxes = begin
        @stateflux SM ~ snow_release + capillary - evap - recharge
        @stateflux UZ ~ recharge - capillary - q0 - perc
        @stateflux LZ ~ perc - q1
    end
end

model = @hydromodel :raven_hbv_light begin
    snowfall_flux
    rainfall_flux
    snow_bucket
    soil_bucket
end

build_raven_hbv_light() = model

end



