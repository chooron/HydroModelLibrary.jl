module cemaneigegr4j
using ..HydroModels

# Model variables
@variables P [description = "Precipitation input", unit = "mm/d"]
@variables Ep [description = "Potential evapotranspiration input", unit = "mm/d"]
@variables T [description = "Mean temperature input", unit = "°C"]

# Cemaneige snow module variables
@variables frac_solid [description = "Fraction of solid precipitation", unit = "-"]
@variables snow [description = "Solid precipitation", unit = "mm/d"]
@variables rain [description = "Liquid precipitation", unit = "mm/d"]
@variables pot_melt [description = "Potential snowmelt", unit = "mm/d"]
@variables G_ratio [description = "Snow cover ratio", unit = "-"]
@variables melt [description = "Actual snowmelt", unit = "mm/d"]
@variables liquid_water [description = "Liquid water output from snow module", unit = "mm/d"]
@variables G [description = "Snow pack storage", unit = "mm"]
@variables eTG [description = "Thermal state of snow pack", unit = "°C"]

# GR4J variables
@variables ps [description = "Net rainfall that enters the production store", unit = "mm/d"]
@variables pn [description = "Net rainfall (precipitation minus evapotranspiration when positive)", unit = "mm/d"]
@variables es [description = "Actual evaporation from the production store", unit = "mm/d"]
@variables en [description = "Net evaporation (evapotranspiration minus precipitation when positive)", unit = "mm/d"]
@variables perc [description = "Percolation from production store to routing", unit = "mm/d"]
@variables Q9 [description = "Slow flow component", unit = "mm/d"]
@variables Q1 [description = "Fast flow component", unit = "mm/d"]
@variables Q9_routed [description = "Slow flow component routed through unit hydrograph 1", unit = "mm/d"]
@variables Q1_routed [description = "Fast flow component routed through unit hydrograph 2", unit = "mm/d"]
@variables Qroute [description = "Outflow from routing store", unit = "mm/d"]
@variables Qt [description = "Total runoff", unit = "mm/d"]
@variables exch [description = "Water exchange between groundwater and surface water", unit = "mm/d"]
@variables S [description = "Production store level", unit = "mm"]
@variables R [description = "Routing store level", unit = "mm"]

@variables t
model_variables = [P, Ep, T, frac_solid, snow, rain, pot_melt, G_ratio, melt, liquid_water,
                   G, eTG, ps, pn, es, en, perc, Q9, Q1, Q9_routed, Q1_routed, Qroute, Qt, exch, S, R]

# Model parameters
# Cemaneige parameters
@parameters CTG [description = "Thermal state weighting coefficient", bounds = (0, 1), unit = "-"]
@parameters Kf [description = "Degree-day melt factor", bounds = (0, 10), unit = "mm/°C/d"]
@parameters G_thresh [description = "Snow pack threshold for snow cover calculation", bounds = (1, 1000), unit = "mm"]

# GR4J parameters
@parameters x1 [description = "Maximum soil moisture storage", bounds = (1, 2000), unit = "mm"]
@parameters x2 [description = "Subsurface water exchange", bounds = (-20, 20), unit = "mm/d"]
@parameters x3 [description = "Routing store depth", bounds = (1, 300), unit = "mm"]
@parameters x4 [description = "Unit Hydrograph time base", bounds = (0, 10), unit = "d"]

model_parameters = [CTG, Kf, G_thresh, x1, x2, x3, x4]

# Cemaneige snow module
snow_module = @hydrobucket :snow_module begin
    fluxes = begin
        # Calculate fraction of solid precipitation based on temperature
        # Linear transition between 0°C and 2°C: clamp to [0, 1]
        @hydroflux frac_solid ~ clamp((2.0 - T) / 2.0, 0.0, 1.0)

        # Split precipitation into solid and liquid
        @hydroflux snow ~ P * frac_solid
        @hydroflux rain ~ P * (1.0 - frac_solid)

        # Calculate potential melt (only when temperature is positive)
        # Limited by available snow pack
        @hydroflux pot_melt ~ max(0.0, min(Kf * max(0.0, T), G))

        # Calculate snow cover ratio
        @hydroflux G_ratio ~ min(G / G_thresh, 1.0)

        # Calculate actual melt with snow cover effect
        @hydroflux melt ~ (0.9 * G_ratio + 0.1) * pot_melt

        # Liquid water output: rain + melt
        @hydroflux liquid_water ~ rain + melt
    end
    dfluxes = begin
        # Snow pack accumulation and melt
        @stateflux G ~ snow - melt

        # Thermal state update (exponential smoothing with temperature)
        # eTG converges to min(T, 0) to ensure eTG <= 0
        # This implements: eTG_new = min(CTG * eTG_old + (1-CTG) * T, 0)
        @stateflux eTG ~ (1.0 - CTG) * (min(T, 0.0) - eTG)
    end
end

# GR4J production store (using liquid_water instead of P)
bucket1 = @hydrobucket :bucket1 begin
    fluxes = begin
        @hydroflux pn ~ max(0.0, liquid_water - Ep)
        @hydroflux en ~ max(0.0, Ep - liquid_water)
        @hydroflux ps ~ x1 * (1 - (S / x1)^2) * tanh(pn / x1) / (1 + S / x1 * tanh(pn / x1))
        @hydroflux es ~ S * (2 - S / x1) * tanh(en / x1) / (1 + (1 - S / x1) * tanh(en / x1))
        @hydroflux perc ~ S * (1 - (1 + ((4 / 9) * (S / x1))^4)^(-0.25))
    end
    dfluxes = begin
        @stateflux S ~ ps - es - perc
    end
end

# Split flux
split_flux = @hydroflux begin
    Q9 ~ (perc + pn - ps) * 0.9
    Q1 ~ (perc + pn - ps) * 0.1
end

# Unit hydrographs
uh_1 = @unithydro begin
    uh_func = begin
        x4 => (t / x4)^2.5
    end
    uh_vars = Q9 => Q9_routed
end

uh_2 = @unithydro begin
    uh_func = begin
        2x4 => (1 - 0.5 * (2 - t / x4)^2.5)
        x4 => (0.5 * (t / x4)^2.5)
    end
    uh_vars = Q1 => Q1_routed
end

# GR4J routing store
bucket2 = @hydrobucket :bucket2 begin
    fluxes = begin
        @hydroflux exch ~ x2 * clamp(R / x3, 0, 1)^3.5
        @hydroflux Qroute ~ R * clamp(1 - (1 + clamp(R / x3, 0, Inf)^4)^(-0.25), 0, 1)
        @hydroflux Qt ~ Qroute + max(Q1_routed + exch, 0.0)
    end
    dfluxes = begin
        @stateflux R ~ Q9_routed + exch - Qroute
    end
end

# Complete model: Cemaneige + GR4J
model = @hydromodel :cemaneigegr4j begin
    snow_module
    bucket1
    split_flux
    uh_1
    uh_2
    bucket2
end

export model
end
