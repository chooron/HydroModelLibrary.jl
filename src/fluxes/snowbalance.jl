module SnowBalance

using ..HydroModels
using ..HydroModels: @parameters, @variables, Number

export SNOBAL_SIMPLE_MELT,
       SNOBAL_HBV,
       SNOBAL_HMETS,
       SNOBAL_CEMA_NIEGE,
       SNOBAL_CEMA_NEIGE,
       SNOBAL_COLD_CONTENT,
       SNOBAL_TWO_LAYER

"""Simple snow balance using a cumulative-melt melt factor."""
function SNOBAL_SIMPLE_MELT(;
    snowstorage::Number=first(@variables snowstorage),
    snowmelt::Number=first(@variables snowmelt),
    cumulmelt::Number=first(@variables cumulmelt),
    overflow::Number=first(@variables overflow),
    potential_melt::Number=first(@variables potential_melt),
    snowfall::Number=first(@variables snowfall),
    rainfall::Number=first(@variables rainfall),
    temp::Number=first(@variables temp),
    ma::Number=first(@variables ma),
    mamax::Number=first(@parameters mamax [description = "Maximum melt factor", bounds = (0, 20), unit = "mm/(d*degC)"]),
    mamin::Number=first(@parameters mamin [description = "Minimum melt factor", bounds = (0, 20), unit = "mm/(d*degC)"]),
    malpha::Number=first(@parameters malpha [description = "Cumulative melt factor coefficient", bounds = (0, 1), unit = "-"]),
    tbm::Number=first(@parameters tbm [description = "Melt temperature threshold", bounds = (-10, 10), unit = "degC"]),
    name::Union{Nothing, Symbol}=nothing,
)
    return HydroModels.@hydrobucket name begin
        fluxes = begin
            @hydroflux ma ~ min(mamax, mamin * (1 + malpha * max(0.0, cumulmelt)))
            @hydroflux potential_melt ~ ma * max(0.0, temp - tbm)
            @hydroflux snowmelt ~ min(max(0.0, snowstorage), potential_melt)
            @hydroflux overflow ~ snowmelt + rainfall
        end
        dfluxes = begin
            @stateflux snowstorage ~ snowfall - snowmelt
            @stateflux cumulmelt ~ ifelse(snowstorage > 0, snowmelt, -cumulmelt)
        end
    end
end

"""HBV snow balance."""
function SNOBAL_HBV(;
    snowstorage::Number=first(@variables snowstorage),
    liquidstorage::Number=first(@variables liquidstorage),
    temp::Number=first(@variables temp [description = "Air temperature"]),
    potential_melt::Number=first(@variables potential_melt),
    cumulmelt::Number=first(@variables cumulmelt),
    refreeze::Number=first(@variables refreeze),
    overflow::Number=first(@variables overflow),
    snowmelt::Number=first(@variables snowmelt),
    snowfall::Number=first(@variables snowfall),
    rainfall::Number=first(@variables rainfall),
    ma::Number=first(@variables ma),
    ka::Number=first(@parameters ka [description = "Refreezing coefficient", bounds = (0, 10), unit = "mm/(d*degC)"]),
    mamax::Number=first(@parameters mamax [description = "Maximum melt factor", bounds = (0, 20), unit = "mm/(d*degC)"]),
    mamin::Number=first(@parameters mamin [description = "Minimum melt factor", bounds = (0, 20), unit = "mm/(d*degC)"]),
    malpha::Number=first(@parameters malpha [description = "Cumulative melt factor coefficient", bounds = (0, 1), unit = "-"]),
    tbm::Number=first(@parameters tbm [description = "Melt temperature threshold", bounds = (-10, 10), unit = "degC"]),
    tbf::Number=first(@parameters tbf [description = "Refreezing temperature threshold", bounds = (-10, 10), unit = "degC"]),
    const_swi::Number=first(@parameters const_swi [description = "Maximum snow liquid water fraction", bounds = (0, 1), unit = "-"]),
    name::Union{Nothing, Symbol}=nothing,
)
    return HydroModels.@hydrobucket name begin
        fluxes = begin
            @hydroflux ma ~ min(mamax, mamin * (1 + malpha * max(0.0, cumulmelt)))
            @hydroflux potential_melt ~ ma * max(0.0, temp - tbm)
            @hydroflux snowmelt ~ min(max(0.0, snowstorage), potential_melt)
            @hydroflux refreeze ~ min(max(0.0, liquidstorage), ka * max(tbf - temp, 0.0))
            @hydroflux overflow ~ max(0.0, liquidstorage + rainfall + snowmelt - const_swi * max(0.0, snowstorage))
        end
        dfluxes = begin
            @stateflux snowstorage ~ snowfall - snowmelt + refreeze
            @stateflux liquidstorage ~ snowmelt + rainfall - refreeze - overflow
            @stateflux cumulmelt ~ ifelse(snowstorage > 0, snowmelt, -cumulmelt)
        end
    end
end

"""HMETS snow balance."""
function SNOBAL_HMETS(;
    snowstorage::Number=first(@variables snowstorage),
    liquidstorage::Number=first(@variables liquidstorage),
    cumulmelt::Number=first(@variables cumulmelt),
    tmean::Number=first(@variables tmean),
    potential_melt::Number=first(@variables potential_melt),
    swi::Number=first(@variables swi),
    refreeze::Number=first(@variables refreeze),
    snowmelt::Number=first(@variables snowmelt),
    snowfall::Number=first(@variables snowfall),
    rainfall::Number=first(@variables rainfall),
    overflow::Number=first(@variables overflow),
    ma::Number=first(@variables ma),
    alpha::Number=first(@parameters alpha [description = "Snow water retention decline coefficient", bounds = (0, 1), unit = "-"]),
    swimax::Number=first(@parameters swimax [description = "Maximum snow liquid water fraction", bounds = (0, 1), unit = "-"]),
    swimin::Number=first(@parameters swimin [description = "Minimum snow liquid water fraction", bounds = (0, 1), unit = "-"]),
    kf::Number=first(@parameters kf [description = "Refreezing coefficient", bounds = (0, 10), unit = "mm/(d*degC^f)"]),
    tbm::Number=first(@parameters tbm [description = "Melt temperature threshold", bounds = (-10, 10), unit = "degC"]),
    tbf::Number=first(@parameters tbf [description = "Refreezing temperature threshold", bounds = (-10, 10), unit = "degC"]),
    f::Number=first(@parameters f [description = "HMETS refreezing exponent", bounds = (0, 5), unit = "-"]),
    mamax::Number=first(@parameters mamax [description = "Maximum melt factor", bounds = (0, 20), unit = "mm/(d*degC)"]),
    mamin::Number=first(@parameters mamin [description = "Minimum melt factor", bounds = (0, 20), unit = "mm/(d*degC)"]),
    malpha::Number=first(@parameters malpha [description = "Cumulative melt factor coefficient", bounds = (0, 1), unit = "-"]),
    name::Union{Nothing, Symbol}=nothing,
)
    return HydroModels.@hydrobucket name begin
        fluxes = begin
            @hydroflux ma ~ min(mamax, mamin * (1 + malpha * max(0.0, cumulmelt)))
            @hydroflux potential_melt ~ ma * max(0.0, tmean - tbm)
            @hydroflux snowmelt ~ min(max(0.0, snowstorage), potential_melt)
            @hydroflux refreeze ~ min(max(0.0, liquidstorage), kf * max(0.0, tbf - tmean)^f)
            @hydroflux swi ~ max(swimin, swimax * (1 - alpha * max(0.0, cumulmelt)))
            @hydroflux overflow ~ max(0.0, liquidstorage + rainfall + snowmelt - swi * max(0.0, snowstorage))
        end
        dfluxes = begin
            @stateflux snowstorage ~ snowfall - snowmelt + refreeze
            @stateflux liquidstorage ~ snowmelt + rainfall - refreeze - overflow
            @stateflux cumulmelt ~ ifelse(snowstorage > 0, snowmelt, -cumulmelt)
        end
    end
end

"""CemaNeige-style snow balance with a thermal-state index and snow-cover scaling."""
function SNOBAL_CEMA_NIEGE(;
    snowstorage::Number=first(@variables snowstorage),
    thermal_state::Number=first(@variables thermal_state),
    snow_cover::Number=first(@variables snow_cover),
    potential_melt::Number=first(@variables potential_melt),
    snowmelt::Number=first(@variables snowmelt),
    overflow::Number=first(@variables overflow),
    snowfall::Number=first(@variables snowfall),
    rainfall::Number=first(@variables rainfall),
    temp::Number=first(@variables temp),
    ctg::Number=first(@parameters ctg [description = "Thermal state weighting coefficient", bounds = (0, 1), unit = "-"]),
    melt_factor::Number=first(@parameters melt_factor [description = "Degree-day melt factor", bounds = (0, 20), unit = "mm/(d*degC)"]),
    g_thresh::Number=first(@parameters g_thresh [description = "Snow-cover threshold storage", bounds = (1.0e-6, 5000), unit = "mm"]),
    name::Union{Nothing, Symbol}=nothing,
)
    return HydroModels.@hydrobucket name begin
        fluxes = begin
            @hydroflux snow_cover ~ clamp(max(0.0, snowstorage) / max(g_thresh, 1.0e-12), 0.0, 1.0)
            @hydroflux potential_melt ~ max(0.0, min(melt_factor * max(0.0, temp), max(0.0, snowstorage)))
            @hydroflux snowmelt ~ (0.9 * snow_cover + 0.1) * potential_melt
            @hydroflux overflow ~ max(0.0, rainfall) + snowmelt
        end
        dfluxes = begin
            @stateflux snowstorage ~ snowfall - snowmelt
            @stateflux thermal_state ~ (1 - ctg) * (min(temp, 0.0) - thermal_state)
        end
    end
end

"""Alias matching the Raven Chapter 3 process name."""
function SNOBAL_CEMA_NEIGE(; kwargs...)
    return SNOBAL_CEMA_NIEGE(; kwargs...)
end

"""Cold-content snow balance with explicit energy deficit before melt can occur."""
function SNOBAL_COLD_CONTENT(;
    snowstorage::Number=first(@variables snowstorage),
    liquidstorage::Number=first(@variables liquidstorage),
    cold_content::Number=first(@variables cold_content),
    temp::Number=first(@variables temp),
    snowfall::Number=first(@variables snowfall),
    rainfall::Number=first(@variables rainfall),
    potential_melt::Number=first(@variables potential_melt),
    warming::Number=first(@variables warming),
    snowmelt::Number=first(@variables snowmelt),
    refreeze::Number=first(@variables refreeze),
    overflow::Number=first(@variables overflow),
    tf::Number=first(@parameters tf [description = "Freeze point", bounds = (-10, 10), unit = "degC"]),
    kf::Number=first(@parameters kf [description = "Refreezing coefficient", bounds = (0, 10), unit = "mm/(d*degC)"]),
    cold_content_coeff::Number=first(@parameters cold_content_coeff [description = "Cold-content accumulation coefficient", bounds = (0, 100), unit = "mm/(d*degC)"]),
    swi::Number=first(@parameters swi [description = "Maximum snow liquid water fraction", bounds = (0, 1), unit = "-"]),
    name::Union{Nothing, Symbol}=nothing,
)
    return HydroModels.@hydrobucket name begin
        fluxes = begin
            @hydroflux warming ~ min(max(0.0, potential_melt), max(0.0, cold_content))
            @hydroflux snowmelt ~ min(max(0.0, snowstorage), max(0.0, potential_melt - warming))
            @hydroflux refreeze ~ min(max(0.0, liquidstorage), kf * max(tf - temp, 0.0))
            @hydroflux overflow ~ max(0.0, liquidstorage + rainfall + snowmelt - swi * max(0.0, snowstorage))
        end
        dfluxes = begin
            @stateflux snowstorage ~ snowfall - snowmelt + refreeze
            @stateflux liquidstorage ~ rainfall + snowmelt - refreeze - overflow
            @stateflux cold_content ~ cold_content_coeff * max(tf - temp, 0.0) - warming
        end
    end
end

"""Two-layer snowpack with upper-to-lower liquid transfer and final outflow from the lower layer."""
function SNOBAL_TWO_LAYER(;
    upper_snowstorage::Number=first(@variables upper_snowstorage),
    lower_snowstorage::Number=first(@variables lower_snowstorage),
    upper_liquidstorage::Number=first(@variables upper_liquidstorage),
    lower_liquidstorage::Number=first(@variables lower_liquidstorage),
    temp::Number=first(@variables temp),
    snowfall::Number=first(@variables snowfall),
    rainfall::Number=first(@variables rainfall),
    potential_melt::Number=first(@variables potential_melt),
    upper_melt::Number=first(@variables upper_melt),
    lower_melt::Number=first(@variables lower_melt),
    upper_refreeze::Number=first(@variables upper_refreeze),
    lower_refreeze::Number=first(@variables lower_refreeze),
    transfer_to_lower::Number=first(@variables transfer_to_lower),
    overflow::Number=first(@variables overflow),
    tf::Number=first(@parameters tf [description = "Freeze point", bounds = (-10, 10), unit = "degC"]),
    kf::Number=first(@parameters kf [description = "Refreezing coefficient", bounds = (0, 10), unit = "mm/(d*degC)"]),
    upper_swi::Number=first(@parameters upper_swi [description = "Upper-layer liquid water fraction", bounds = (0, 1), unit = "-"]),
    lower_swi::Number=first(@parameters lower_swi [description = "Lower-layer liquid water fraction", bounds = (0, 1), unit = "-"]),
    name::Union{Nothing, Symbol}=nothing,
)
    return HydroModels.@hydrobucket name begin
        fluxes = begin
            @hydroflux upper_melt ~ min(max(0.0, upper_snowstorage), max(0.0, potential_melt))
            @hydroflux lower_melt ~ min(max(0.0, lower_snowstorage), max(0.0, potential_melt - upper_melt))
            @hydroflux upper_refreeze ~ min(max(0.0, upper_liquidstorage), kf * max(tf - temp, 0.0))
            @hydroflux lower_refreeze ~ min(max(0.0, lower_liquidstorage), kf * max(tf - temp, 0.0))
            @hydroflux transfer_to_lower ~ max(0.0, upper_liquidstorage + rainfall + upper_melt - upper_refreeze - upper_swi * max(0.0, upper_snowstorage))
            @hydroflux overflow ~ max(0.0, lower_liquidstorage + transfer_to_lower + lower_melt - lower_refreeze - lower_swi * max(0.0, lower_snowstorage))
        end
        dfluxes = begin
            @stateflux upper_snowstorage ~ snowfall - upper_melt + upper_refreeze
            @stateflux lower_snowstorage ~ -lower_melt + lower_refreeze
            @stateflux upper_liquidstorage ~ rainfall + upper_melt - upper_refreeze - transfer_to_lower
            @stateflux lower_liquidstorage ~ transfer_to_lower + lower_melt - lower_refreeze - overflow
        end
    end
end

end

