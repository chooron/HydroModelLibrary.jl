module Seepage

using ..HydroModels
using ..HydroModels: @parameters, @variables, Number

export SEEPAGE_LINEAR,
       SEEP_LINEAR,
       SEEPAGE_THRESHOLD

"""Linear seepage from depression, wetland, or groundwater storage."""
function SEEPAGE_LINEAR(;
    seepage::Number=first(@variables seepage),
    waterstorage::Number=first(@variables waterstorage),
    datum_storage::Number=first(@parameters datum_storage [description = "Reference datum storage", bounds = (-5000, 5000), unit = "mm"]),
    seepage_coeff::Number=first(@parameters seepage_coeff [description = "Linear seepage coefficient", bounds = (0, 100), unit = "mm/d"]),
    flux_name::Symbol=:seepage_linear,
)
    @hydroflux flux_name seepage ~ max(0.0, seepage_coeff * (waterstorage - datum_storage))
end

"""Alias matching the Raven Chapter 3 process name."""
function SEEP_LINEAR(; kwargs...)
    params = (; kwargs...)
    return SEEPAGE_LINEAR(; params..., flux_name=get(params, :flux_name, :seep_linear))
end

"""Threshold seepage above a storage threshold."""
function SEEPAGE_THRESHOLD(;
    seepage::Number=first(@variables seepage),
    waterstorage::Number=first(@variables waterstorage),
    seepage_coeff::Number=first(@parameters seepage_coeff [description = "Seepage coefficient", bounds = (0, 100), unit = "mm/d"]),
    storage_threshold::Number=first(@parameters storage_threshold [description = "Seepage threshold", bounds = (-5000, 5000), unit = "mm"]),
    flux_name::Symbol=:seepage_threshold,
)
    @hydroflux flux_name seepage ~ max(0.0, seepage_coeff * max(waterstorage - storage_threshold, 0.0))
end

end
