module hmets

using ..HydroModels
using ..GammaUnitHydro
#--------------------------------------------------------------------------------
# Parameters
#--------------------------------------------------------------------------------

# Snow parameters
@parameters ddfmin [description = "Minimum degree-day-factor in mm/°C/day", bounds = (0.0, 20.0)]
@parameters ddfplus [description = "Maximum degree-day-factor in mm/°C/day (ddfmin + ddfplus = ddfmax)", bounds = (0.0, 20.0)]
@parameters Tbm [description = "Base melting temperature in °C", bounds = (-2.0, 3.0)]
@parameters Kcum [description = "Empirical parameter for the calculation of the degree-day-factor in mm⁻¹", bounds = (0.01, 0.2)]
@parameters fcmin [description = "Minimum fraction for the snowpack water retention capacity", bounds = (0.0, 0.1)]
@parameters fcplus [description = "Maximum fraction of the snowpack water retention capacity (fcmin + fcplus = fcmax)", bounds = (0.01, 0.25)]
@parameters Ccum [description = "Parameter for the calculation of water retention capacity in mm⁻¹", bounds = (0.005, 0.05)]
@parameters Tbf [description = "Base refreezing temperature in °C", bounds = (-5.0, 2.0)]
@parameters Kf [description = "Degree-day factor for refreezing in mm/°C/day", bounds = (0.0, 5.0)]
@parameters Fe [description = "Empirical exponent for the freezing equation", bounds = (0.0, 1.0)]

# Real evapotranspiration parameter
@parameters ETeff [description = "Fraction of the potential evapotranspiration", bounds = (0.0, 3.0)]

# Subsurface parameters
@parameters cr [description = "Fraction of the water for surface and delayed runoff", bounds = (0.0, 1.0)]
@parameters cvp [description = "Fraction of the water for groundwater recharge", bounds = (0.00001, 0.02)]
@parameters cv [description = "Fraction of the water for hypodermic flow", bounds = (0.0, 0.1)]
@parameters cp [description = "Fraction of the water for groundwater flow", bounds = (0.00001, 0.01)]
@parameters LVmax [description = "Maximum level of the vadose zone in mm", bounds = (0.001, 500.0)]
@parameters LPmax [description = "Maximum level of the phreatic zone in mm", bounds = (0.001, 500.0)]

# Unit hydrograph parameters
@parameters α1 [description = "Shape parameter for the gamma distribution used on the surface unit hydrograph", bounds = (0.0, 1.0)]
@parameters β1 [description = "Rate parameter for the gamma distribution used on the surface unit hydrograph", bounds = (0.0, 1.0)]
@parameters α2 [description = "Shape parameter for the gamma distribution used on the delayed unit hydrograph", bounds = (0.0, 1.0)]
@parameters β2 [description = "Rate parameter for the gamma distribution used on the delayed unit hydrograph", bounds = (0.0, 1.0)]
model_parameters = [ddfmin, ddfplus, Tbm, Kcum, fcmin, fcplus, Ccum, Tbf, Kf, Fe, ETeff, cr, cvp, cv, cp, LVmax, LPmax, α1, β1, α2, β2]
#--------------------------------------------------------------------------------
# Variables
#--------------------------------------------------------------------------------
# Forcings
@variables T [description = "Mean daily temperature in °C"]
@variables P [description = "Precipitation in mm/day"]
@variables SNF [description = "Snowfall in mm/day"]
@variables RAIN [description = "Rainfall in mm/day"]
@variables Ep [description = "Potential Evapotranspiration in mm/day"]

# State variables
@variables SNW [description = "Snow water equivalent in mm"]
@variables LV [description = "Vadose zone level in mm"]
@variables LP [description = "Phreatic zone level in mm"]

# Fluxes and other internal variables
@variables TD [description = "Meandiurnal temperature in °C"]
@variables POR [description = "Potential amount of overnight refreezing in mm"]
@variables DDF [description = "Actual degree day factor in mm/°C/day"]
@variables PSM [description = "Potential snowmelt in mm"]
@variables WRF [description = "Water retention fraction"]
@variables WAR [description = "Water available for runoff and infiltration in mm"]
@variables RET [description = "Real Evapotranspiration in mm"]
@variables INF [description = "Infiltration in mm"]
@variables GR [description = "Groundwater recharge in mm"]
@variables H1 [description = "Surface runoff component from rainfall in mm"]
@variables H2 [description = "Surface runoff component from snowmelt in mm"]
@variables H2_routed [description = "Surface runoff component from snowmelt routed through unit hydrograph in mm"]
@variables H3 [description = "Hypodermic flow in mm"]
@variables H4 [description = "Groundwater flow in mm"]
@variables H3_routed [description = "Hypodermic flow routed through unit hydrograph in mm"]
@variables Qt [description = "Total runoff in mm"]
model_variables = [T, P, SNF, RAIN, Ep, SNW, LV, LP, TD, POR, DDF, PSM, WRF, WAR, RET, INF, GR, H1, H2, H3, H4, Qt, H2_routed, H3_routed]
#--------------------------------------------------------------------------------
# Model definition
#--------------------------------------------------------------------------------

# Parameter-derived values
ddfmax = ddfmin + ddfplus
fcmax = fcmin + fcplus

snow_bucket = @hydrobucket begin
    fluxes = begin
        @hydroflux begin
            RAIN ~ step_func(T) * P
            SNF ~ step_func(-T) * P
        end
        @hydroflux POR ~ max(0, Kf * max(0,Tbf - (T / 2))^Fe)
        @hydroflux DDF ~ ddfmin * (1 + Kcum * SNW)
        @hydroflux PSM ~ max(0, DDF * (T - Tbm))
        @hydroflux WRF ~ max(fcmin, fcmax * (1 - Ccum * SNW))
        @hydroflux WAR ~ max(0, SNW + RAIN - WRF)
    end

    dfluxes = begin
        @stateflux SNW ~ POR + SNF - PSM
    end
end

soil_bucket = @hydrobucket begin
    fluxes = begin
        @hydroflux RET ~ ETeff * Ep
        @hydroflux H1 ~ min(LV, cr * (LV / LVmax) * WAR)
        @hydroflux INF ~ WAR - H1 - RET
        @hydroflux GR ~ cvp * LV
        @hydroflux H2 ~ min(LV, cr * INF * (LV / LVmax))
        @hydroflux H3 ~ cv * LV
        @hydroflux H4 ~ cp * LP
    end

    dfluxes = begin
        @stateflux LV ~ INF - RET - H1 - H2 - H3 - GR
        @stateflux LP ~ GR - H4
    end
end

gamma_uh1 = GammaHydrograph([H2], [H2_routed], α1, β1)
gamma_uh2 = GammaHydrograph([H3], [H3_routed], α2, β2)
q_flux = @hydroflux Qt ~ H2_routed + H3_routed + H4 + H1

model = @hydromodel begin
    snow_bucket
    soil_bucket
    gamma_uh1
    gamma_uh2
    q_flux
end

# todo unit hydrograph
end