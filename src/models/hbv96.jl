module hbv96
using ..HydroModels
using ..HydroModels: step_func

@variables P [description = "Precipitation", unit = "mm/d"]
@variables T [description = "Temperature", unit = "°C"]
@variables Ep [description = "Potential evapotranspiration", unit = "mm/d"]

@variables SP [description = "Snow storage", unit = "mm"]
@variables WC [description = "Liquid water in snowpack", unit = "mm"]
@variables SM [description = "Soil moisture", unit = "mm"]
@variables UZ [description = "Upper zone storage", unit = "mm"]
@variables LZ [description = "Lower zone storage", unit = "mm"]

@variables sf rf refr melt infi excess evap cf r perc q0 q1 q Qt t
model_variables = [P, T, Ep, SP, WC, SM, UZ, LZ, sf, rf, refr, melt, infi, excess, evap, cf, r, perc, q0, q1, q, Qt, t]

@parameters TT [description = "Threshold temperature for snowfall", bounds = (-3.0, 5.0), unit = "°C"]
@parameters TTI [description = "Interval length of rain-snow spectrum", bounds = (0.0, 17.0), unit = "°C"]
@parameters TTM [description = "Threshold temperature for snowmelt", bounds = (-3.0, 3.0), unit = "°C"]
@parameters CFR [description = "Coefficient of refreezing of melted snow", bounds = (0.0, 1.0), unit = "-"]
@parameters CFMAX [description = "Degree-day factor of snowmelt and refreezing", bounds = (0.0, 20.0), unit = "mm/°C/d"]
@parameters WHC [description = "Maximum water holding content of snow pack", bounds = (0.0, 1.0), unit = "-"]
@parameters CFLUX [description = "Maximum rate of capillary rise", bounds = (0.0, 4.0), unit = "mm/d"]
@parameters FC [description = "Maximum soil moisture storage", bounds = (1.0, 2000.0), unit = "mm"]
@parameters LP [description = "Wilting point as fraction of FC", bounds = (0.05, 0.95), unit = "-"]
@parameters BETA [description = "Non-linearity coefficient of upper zone recharge", bounds = (0.0, 10.0), unit = "-"]
@parameters K0 [description = "Runoff coefficient from upper zone", bounds = (0.0, 1.0), unit = "d⁻¹"]
@parameters ALPHA [description = "Non-linearity coefficient of runoff from upper zone", bounds = (0.0, 4.0), unit = "-"]
@parameters PERC [description = "Maximum rate of percolation to lower zone", bounds = (0.0, 20.0), unit = "mm/d"]
@parameters K1 [description = "Runoff coefficient from lower zone", bounds = (0.0, 1.0), unit = "d⁻¹"]
model_parameters = [TT, TTI, TTM, CFR, CFMAX, WHC, CFLUX, FC, LP, BETA, K0, ALPHA, PERC, K1]

snow_bucket = @hydrobucket :hbv_snow begin
    fluxes = begin
        @hydroflux sf ~ P * clamp((TT + TTI / 2) - T, 0, 1)
        @hydroflux rf ~ P * clamp(T - (TT + TTI / 2), 0, 1)
        @hydroflux refr ~ clamp(CFR * CFMAX * (TTM - T), 0, WC)
        @hydroflux melt ~ clamp(CFMAX * (T - TTM), 0, SP)
        @hydroflux infi ~ step_func(WC - WHC * SP) * (rf + melt)
        @hydroflux excess ~ max(0.0, WC - WHC * SP)
    end
    dfluxes = begin
        @stateflux SP ~ sf + refr - melt
        @stateflux WC ~ rf + melt - refr - infi - excess
    end
end

soil_bucket = @hydrobucket :hbv_soil begin
    fluxes = begin
        @hydroflux cf ~ min(UZ, CFLUX * (1 - SM / FC))
        @hydroflux evap ~ min(Ep * min(1.0, SM / (FC * LP)), SM)
        @hydroflux r ~ (infi + excess) * clamp(SM / FC, 0, 1)^BETA
        @hydroflux q0 ~ min(K0 * (max(0.0, UZ)^(1 + ALPHA)), UZ)
        @hydroflux perc ~ min(PERC, UZ)
        @hydroflux q1 ~ K1 * LZ
        @hydroflux q ~ q0 + q1
    end
    dfluxes = begin
        @stateflux SM ~ infi + excess + cf - evap - r
        @stateflux UZ ~ r - cf - q0 - perc
        @stateflux LZ ~ perc - q1
    end
end

model = @hydromodel :hbv96 begin
    snow_bucket
    soil_bucket
end

end