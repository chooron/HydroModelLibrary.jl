module Interflow

using ..HydroModels
using ..HydroModels: @parameters, @variables, Number

export INTERFLOW_PRMS

"""PRMS-style storage-dependent interflow release."""
function INTERFLOW_PRMS(;
    interflow::Number=first(@variables interflow),
    waterstorage::Number=first(@variables waterstorage),
    field_capacity::Number=first(@parameters field_capacity [description = "Field capacity threshold", bounds = (0, 5000), unit = "mm"]),
    interflow_coeff::Number=first(@parameters interflow_coeff [description = "Interflow coefficient", bounds = (0, 10), unit = "d-1"]),
    slope_factor::Number=first(@parameters slope_factor [description = "Slope or travel-time scaling factor", bounds = (0, 100), unit = "-"]),
    flux_name::Symbol=:interflow_prms,
)
    drainable_storage = max(max(0.0, waterstorage) - field_capacity, 0.0)
    raw_interflow = interflow_coeff * slope_factor * drainable_storage
    @hydroflux flux_name interflow ~ clamp(raw_interflow, 0.0, max(0.0, waterstorage))
end

end