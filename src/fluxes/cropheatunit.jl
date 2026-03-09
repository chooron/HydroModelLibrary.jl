module CropHeatUnit

using ..HydroModels
using ..HydroModels: @parameters, @variables, Number

export CHU_ONTARIO

function CHU_ONTARIO(;
    chu_increment::Number=first(@variables chu_increment),
    max_temp::Number=first(@variables max_temp),
    min_temp::Number=first(@variables min_temp),
    flux_name::Symbol=:chu_ontario,
)
    chu_day = 3.33 * (max_temp - 10) - 0.084 * (max_temp - 10)^2
    chu_night = 1.8 * (min_temp - 4.4)
    @hydroflux flux_name chu_increment ~ max(0.0, 0.5 * (chu_day + chu_night))
end

end
