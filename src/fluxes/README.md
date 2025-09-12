# Hydrological Fluxes Implementation

This module implements hydrological flux calculations based on the Raven hydrological model, providing core algorithms for various hydrological processes.

## Features

- **Snowmelt Calculation**: Potential snowmelt computation using degree-day method
- **Evapotranspiration**: Potential evapotranspiration correction factor calculation
- **Temperature Control**: Freeze/melt temperature threshold handling
- **Parameterization**: Flexible parameter configuration with boundary constraints

## Module Structure

- `PotentialMelt`: Potential snowmelt calculation module
  - `POTMELT_DEGREE_DAY`: Degree-day method snowmelt calculation function