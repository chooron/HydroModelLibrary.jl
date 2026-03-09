module LakeRelease

using ..HydroModels
using ..HydroModels: @parameters, @variables, Number

export LAKEREL_LINEAR

"""Linear lake release from lake storage."""
function LAKEREL_LINEAR(;
    lake_release::Number=first(@variables lake_release),
    lake_storage::Number=first(@variables lake_storage),
    release_coeff::Number=first(@parameters release_coeff [description = "Lake release coefficient", bounds = (0, 10), unit = "d-1"]),
    flux_name::Symbol=:lakerel_linear,
)
    @hydroflux flux_name lake_release ~ max(0.0, release_coeff * lake_storage)
end

end
