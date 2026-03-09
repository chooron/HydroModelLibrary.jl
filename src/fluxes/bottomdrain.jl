module BottomDrain

using ..HydroModels
using ..HydroModels: @parameters, @variables, Number

export BOTTOMDRAIN_LINEAR, BOTTOMDRAIN_POWER, BOTTOMDRAIN_THRESH

"""Linear bottom drainage."""
function BOTTOMDRAIN_LINEAR(;
    bottom_drainage::Number=first(@variables bottom_drainage),
    waterstorage::Number=first(@variables waterstorage),
    drain_coeff::Number=first(@parameters drain_coeff [description = "Bottom drainage coefficient", bounds = (0, 10), unit = "d-1"]),
    flux_name::Symbol=:bottomdrain_linear,
)
    @hydroflux flux_name bottom_drainage ~ clamp(drain_coeff * max(0.0, waterstorage), 0.0, max(0.0, waterstorage))
end

"""Power-law bottom drainage."""
function BOTTOMDRAIN_POWER(;
    bottom_drainage::Number=first(@variables bottom_drainage),
    waterstorage::Number=first(@variables waterstorage),
    drain_coeff::Number=first(@parameters drain_coeff [description = "Bottom drainage coefficient", bounds = (0, 100), unit = "mm^(1-n)/d"]),
    drain_n::Number=first(@parameters drain_n [description = "Bottom drainage exponent", bounds = (0, 10), unit = "-"]),
    flux_name::Symbol=:bottomdrain_power,
)
    @hydroflux flux_name bottom_drainage ~ clamp(drain_coeff * max(0.0, waterstorage)^drain_n, 0.0, max(0.0, waterstorage))
end

"""Threshold bottom drainage."""
function BOTTOMDRAIN_THRESH(;
    bottom_drainage::Number=first(@variables bottom_drainage),
    waterstorage::Number=first(@variables waterstorage),
    drain_coeff::Number=first(@parameters drain_coeff [description = "Bottom drainage coefficient", bounds = (0, 10), unit = "d-1"]),
    storage_threshold::Number=first(@parameters storage_threshold [description = "Bottom drainage threshold", bounds = (0, 5000), unit = "mm"]),
    flux_name::Symbol=:bottomdrain_thresh,
)
    raw_drainage = drain_coeff * max(max(0.0, waterstorage) - storage_threshold, 0.0)
    @hydroflux flux_name bottom_drainage ~ clamp(raw_drainage, 0.0, max(0.0, waterstorage))
end

end