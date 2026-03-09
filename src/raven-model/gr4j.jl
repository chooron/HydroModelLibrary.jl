module RavenGR4J

using ..HydroModels
using ..HydroModelLibrary: INF_GR4J, SOILEVAP_GR4J, BASE_GR4J

export build_raven_gr4j, model, model_parameters, model_variables

@variables P [description = "Precipitation", unit = "mm/d"]
@variables Ep [description = "Potential evapotranspiration", unit = "mm/d"]
@variables S [description = "Production store storage", unit = "mm"]
@variables R [description = "Routing store storage", unit = "mm"]

@variables pn [description = "Net rainfall", unit = "mm/d"]
@variables en [description = "Net evaporation", unit = "mm/d"]
@variables ps [description = "Production-store inflow", unit = "mm/d"]
@variables es [description = "Production-store evaporation", unit = "mm/d"]
@variables perc [description = "Percolation from production store", unit = "mm/d"]
@variables peff [description = "Effective rainfall leaving production store", unit = "mm/d"]
@variables Q9 [description = "Slow branch inflow", unit = "mm/d"]
@variables Q1 [description = "Fast branch inflow", unit = "mm/d"]
@variables Q9_routed [description = "Routed slow branch inflow", unit = "mm/d"]
@variables Q1_routed [description = "Routed fast branch inflow", unit = "mm/d"]
@variables exch [description = "Groundwater exchange", unit = "mm/d"]
@variables Qroute [description = "Routing store outflow", unit = "mm/d"]
@variables Qt [description = "Total runoff", unit = "mm/d"]
@variables t

model_variables = [P, Ep, S, R, pn, en, ps, es, perc, peff, Q9, Q1, Q9_routed, Q1_routed, exch, Qroute, Qt, t]

@parameters x1 [description = "Production store capacity", bounds = (1.0, 2000.0), unit = "mm"]
@parameters x2 [description = "Groundwater exchange coefficient", bounds = (-20.0, 20.0), unit = "mm/d"]
@parameters x3 [description = "Routing store capacity", bounds = (1.0, 300.0), unit = "mm"]
@parameters x4 [description = "Unit hydrograph time base", bounds = (0.0, 10.0), unit = "d"]

model_parameters = [x1, x2, x3, x4]

production_bucket = @hydrobucket :raven_gr4j_production begin
    fluxes = begin
        @hydroflux pn ~ max(0.0, P - Ep)
        @hydroflux en ~ max(0.0, Ep - P)
        INF_GR4J(;
            infiltration=ps,
            prcpeff=pn,
            waterstorage=S,
            max_waterstorage=x1,
            flux_name=:raven_gr4j_ps,
        )
        SOILEVAP_GR4J(;
            soil_evaporation=es,
            waterstorage=S,
            potential_evaporation=en,
            max_waterstorage=x1,
            pet_corr=1.0,
            flux_name=:raven_gr4j_es,
        )
        BASE_GR4J(;
            baseflow=perc,
            waterstorage=S,
            reference_waterstorage=x1,
            flux_name=:raven_gr4j_perc,
        )
        @hydroflux peff ~ perc + pn - ps
    end
    dfluxes = begin
        @stateflux S ~ ps - es - perc
    end
end

split_flux = @hydroflux begin
    Q9 ~ 0.9 * peff
    Q1 ~ 0.1 * peff
end

uh_1 = @unithydro begin
    uh_func = begin
        x4 => (t / x4)^2.5
    end
    uh_vars = Q9 => Q9_routed
end

uh_2 = @unithydro begin
    uh_func = begin
        2x4 => 1 - 0.5 * (2 - t / x4)^2.5
        x4 => 0.5 * (t / x4)^2.5
    end
    uh_vars = Q1 => Q1_routed
end

routing_bucket = @hydrobucket :raven_gr4j_routing begin
    fluxes = begin
        @hydroflux exch ~ x2 * clamp(R / x3, 0.0, Inf)^3.5
        BASE_GR4J(;
            baseflow=Qroute,
            waterstorage=R,
            reference_waterstorage=x3,
            flux_name=:raven_gr4j_qroute,
        )
        @hydroflux Qt ~ Qroute + max(Q1_routed + exch, 0.0)
    end
    dfluxes = begin
        @stateflux R ~ Q9_routed + exch - Qroute
    end
end

model = @hydromodel :raven_gr4j begin
    production_bucket
    split_flux
    uh_1
    uh_2
    routing_bucket
end

build_raven_gr4j() = model

end

