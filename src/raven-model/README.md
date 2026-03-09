# Raven Template Models

This directory assembles Raven Appendix F template models from reusable process
families in `src/fluxes`.

## Implemented

- `hbv_light.jl`: HBV-Light style snow-soil-groundwater template
- `gr4j.jl`: GR4J production and routing template
- `hymod.jl`: HYMOD probability-distributed soil moisture and reservoir chain

## TODO

The remaining Appendix F template families are registered in
`todo_models.jl`. Those entries expose placeholder builder functions so the
missing process gaps are explicit and discoverable from `HydroModelLibrary`.
