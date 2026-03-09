module Overflow

using ..HydroModels
using ..HydroModels: @parameters, @variables, Number

export OVERFLOW_ALL, OVERFLOW_THRESHOLD, OVERFLOW_LINEAR, OVERFLOW_NONLINEAR, OVERFLOW_GR4J

"""Overflow once storage exceeds the maximum."""
function OVERFLOW_ALL(;
    overflow::Number=first(@variables overflow),
    waterstorage::Number=first(@variables waterstorage),
    max_waterstorage::Number=first(@parameters max_waterstorage [description = "Maximum storage", bounds = (0, 5000), unit = "mm"]),
    flux_name::Symbol=:overflow_all,
)
    @hydroflux flux_name overflow ~ max(max(0.0, waterstorage) - max_waterstorage, 0.0)
end

"""Threshold overflow."""
function OVERFLOW_THRESHOLD(;
    overflow::Number=first(@variables overflow),
    waterstorage::Number=first(@variables waterstorage),
    storage_threshold::Number=first(@parameters storage_threshold [description = "Overflow threshold", bounds = (0, 5000), unit = "mm"]),
    flux_name::Symbol=:overflow_threshold,
)
    @hydroflux flux_name overflow ~ max(max(0.0, waterstorage) - storage_threshold, 0.0)
end

"""Linear overflow above a storage threshold."""
function OVERFLOW_LINEAR(;
    overflow::Number=first(@variables overflow),
    waterstorage::Number=first(@variables waterstorage),
    overflow_coeff::Number=first(@parameters overflow_coeff [description = "Overflow coefficient", bounds = (0, 10), unit = "d-1"]),
    storage_threshold::Number=first(@parameters storage_threshold [description = "Overflow threshold", bounds = (0, 5000), unit = "mm"]),
    flux_name::Symbol=:overflow_linear,
)
    raw_overflow = overflow_coeff * max(max(0.0, waterstorage) - storage_threshold, 0.0)
    @hydroflux flux_name overflow ~ clamp(raw_overflow, 0.0, max(0.0, waterstorage))
end

"""Nonlinear overflow above a storage threshold."""
function OVERFLOW_NONLINEAR(;
    overflow::Number=first(@variables overflow),
    waterstorage::Number=first(@variables waterstorage),
    overflow_coeff::Number=first(@parameters overflow_coeff [description = "Overflow coefficient", bounds = (0, 100), unit = "mm^(1-n)/d"]),
    storage_threshold::Number=first(@parameters storage_threshold [description = "Overflow threshold", bounds = (0, 5000), unit = "mm"]),
    overflow_n::Number=first(@parameters overflow_n [description = "Overflow exponent", bounds = (0, 10), unit = "-"]),
    flux_name::Symbol=:overflow_nonlinear,
)
    raw_overflow = overflow_coeff * max(max(0.0, waterstorage) - storage_threshold, 0.0)^overflow_n
    @hydroflux flux_name overflow ~ clamp(raw_overflow, 0.0, max(0.0, waterstorage))
end

"""GR4J routing-store overflow."""
function OVERFLOW_GR4J(;
    overflow::Number=first(@variables overflow),
    waterstorage::Number=first(@variables waterstorage),
    reference_waterstorage::Number=first(@parameters reference_waterstorage [description = "Reference routing storage", bounds = (0, 5000), unit = "mm"]),
    flux_name::Symbol=:overflow_gr4j,
)
    storage_ratio = max(0.0, waterstorage) / max(reference_waterstorage, 1.0e-12)
    raw_overflow = max(0.0, waterstorage) * (1 - (1 + ((4 / 9) * storage_ratio)^4)^(-1 / 4))
    @hydroflux flux_name overflow ~ clamp(raw_overflow, 0.0, max(0.0, waterstorage))
end

end