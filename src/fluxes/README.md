# Hydrological Fluxes

This directory contains reusable Raven-style process flux implementations and a Chapter 3 formula catalog.

## Modules

- `rainsnow.jl`: precipitation partitioning formulations
- `precipinterception.jl`: rainfall/snowfall interception partitioning using LAI-based and exponential-LAI canopy formulations
- `abstraction.jl`: abstraction, interception, and derived residual storage uptake
- `baseflow.jl`: baseflow algorithms such as constant, linear, GR4J, TOPMODEL, VIC and threshold formulations
- `infiltration.jl`: rational, SCS, Green-Ampt, VIC, VIC/ARNO, HBV, PRMS, UBC, GR4J, HMETS, AWBM and PDM infiltration formulations
- `capillaryrise.jl`: HBV-style capillary rise with `CRISE_HBV` alias matching the Raven Chapter 3 naming
- `soilbalance.jl`: generic soil-balance, PRMS-style two-layer, and Penman-style root-zone/deficit buckets 
- `canopyevap.jl`: canopy evaporation formulations including all-available, linear-storage, PRMS, maximum canopy evaporation, and Rutter variants
- `canopydrip.jl`: canopy drip and throughfall formulations including excess-storage, PRMS-style, linear-drainage, and slow-drainage variants
- `openwaterevap.jl`: open-water evaporation formulations including Hargreaves, Penman, basic PET-scaled, riparian, and UWFS variants
- `overflow.jl`: saturation-excess and threshold overflow formulations
- `depressionstorage.jl`: generic depression/wetland storage bucket plus threshold-power, linear, and weir overflow formulations
- `seepage.jl`: linear and threshold seepage formulations, including `SEEP_LINEAR` alias matching the Raven Chapter 3 naming
- `lakerelease.jl`: linear lake-release formulation
- `percolation.jl`: percolation and recharge formulations including GR4J exchange aliases
- `potentialmelt.jl`: prescribed melt, degree-day melt, signed freeze/melt, HBV rain-on-snow, HMETS and restricted melt formulations
- `interflow.jl`: interflow formulations
- `bottomdrain.jl`: bottom-drainage formulations
- `snowmelt.jl`: actual snowmelt helpers derived from potential melt
- `snowsublimation.jl`: Kuzmin, Central Sierra, and bulk-aerodynamic snow sublimation formulations
- `snowfreeze.jl`: degree-day and HMETS refreezing formulations
- `snowalbedo.jl`: UBC, CRHM-Essery, and Baker snow albedo evolution formulations
- `snowbalance.jl`: simple melt, HBV, HMETS, CemaNeige, cold-content, and two-layer snow-balance buckets
- `glaciermelt.jl`: degree-day, potential-melt-driven, and GSMSOCONT-style glacier melt formulations
- `glacierrelease.jl`: linear-storage, analytic linear-storage, and HBV-EC glacier-release formulations
- `glacierinfiltration.jl`: direct, linear, and GSMSOCONT-style glacier infiltration formulations
- `lakefreeze.jl`: basic lake-freezing formulation
- `cropheatunit.jl`: Ontario crop heat unit evolution
- `specialprocesses.jl`: special-process helpers including GR4J and gamma convolution kernels
- `soilevap.jl`: soil evaporation formulations including sequential, CHU, AWBM, HYMOD2, HBV-Oresund, PRMS, Sacramento, UBC, PDM, GR4J and HYPR variants
- `raven_chapter3_process_formula_catalog.md`: normalized Chapter 3 process and formula checklist


