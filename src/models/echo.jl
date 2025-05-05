
module echo
# Model variables
@variables I [description = "current interception storage", unit = "mm"]
@variables T [description = "temperature", unit = "oC"]
@variables P [description = "precipitation", unit = "mm/d"]
@variables Ei [description = "evaporation", unit = "mm/d"]
@variables Pn [description = "net precipitation", unit = "mm/d"]
@variables Ep [description = "the potential rate", unit = "mm/d"]
@variables Hs [description = "current storage in the snow pack", unit = "mm"]
@variables Ps [description = "precipitation-as-snow", unit = "mm/d"]
@variables Fs [description = "refreezing of melted snow", unit = "mm/d"]
@variables Ms [description = "snowmelt", unit = "mm/d"]
@variables Gs [description = "ground-heat flux", unit = "mm/d"]
@variables Hw [description = "current storage of liquid water in the snow pack", unit = "mm"]
@variables Pr [description = "precipitation-as-rain", unit = "mm/d"]
@variables Mw [description = "outflow of melt water", unit = "mm/d"]
@variables S [description = "current storage in the soil moisture zone", unit = "mm"]
@variables Fi [description = "infiltration", unit = "mm/d"]
@variables RD [description = "Dunne-type runoff", unit = "mm/d"]
@variables Et_pot [description = "potential evapotranspiration", unit = "mm/d"]
@variables Et [description = "evapotranspiration", unit = "mm/d"]
@variables L [description = "leakage", unit = "mm/d"]
@variables Peq [description = "equivalent precipitation", unit = "mm/d"]
@variables RH [description = "Horton-type runoff", unit = "mm/d"]
@variables Sfast [description = "current storage in the fast runoff reservoir", unit = "mm"]
@variables Lf [description = "leakage-to-fast-flow", unit = "mm/d"]
@variables Rf [description = "fast runoff", unit = "mm/d"]
@variables Ls [description = "leakage-to-slow-flow", unit = "mm/d"]
@variables Lmax [description = "maximum leakage rate", unit = "mm/d"]
@variables Sslow [description = "current storage in the slow runoff reservoir", unit = "mm"]
@variables Ls [description = "leakage-to-slow-flow", unit = "mm/d"]
@variables Rs [description = "slow runoff", unit = "mm/d"]
@variables Q [description = "Total flow", unit = "mm/d"]


# Model parameters
@parameters rho [description = "Maximum interception storage", bounds = (0,5), unit = "mm"]
@parameters Ts [description = "Threshold temperature for snowfall", bounds = (-3,5), unit = "oC"]
@parameters Tm [description = "Threshold temperature for snowmelt", bounds = (-3,3), unit = "oC"]
@parameters as [description = "Degree-day factor for snowmelt", bounds = (0,20), unit = "mm/oC/d"]
@parameters af [description = "Degree-day factor reduction factor for refreezing", bounds = (0,1), unit = "-"]
@parameters Gmax [description = "Snow melt through ground heat flux rate", bounds = (0,2), unit = "mm/d"]
@parameters the [description = "Water holding capacity as fraction of current snow pack", bounds = (0,1), unit = "-"]
@parameters phi [description = "Maximum Horton type flow rate", bounds = (0,200), unit = "mm/d"]
@parameters Smax [description = "Maximum soil moisture storage", bounds = (1,2000), unit = "mm"]
@parameters sw [description = "Wilting point", bounds = (0.05,0.95), unit = "-"]
@parameters sm [description = "Plant stress point", bounds = (0.05,0.95), unit = "-"]
@parameters Ksat [description = "Runoff coefficient", bounds = (0,1), unit = "d-1"]
@parameters c [description = "Runoff non-linearity", bounds = (0,5), unit = "-"]
@parameters Lmax [description = "Maximum leakage rate", bounds = (0,20), unit = "mm/d"]
@parameters kf [description = "Runoff coefficient", bounds = (0,1), unit = "d-1"]
@parameters ks [description = "Runoff coefficient", bounds = (0,1), unit = "d-1"]

# Soil water component
soil_bucket_1 = @hydrobucket :interception begin
    fluxes = begin
        @hydroflux Ei ~ ifelse(I > 0, Ep, 0)
        @hydroflux Pn ~ ifelse(I == rho, P, 0)
    end

    dfluxes = begin
        @stateflux I ~ P - Ei - Pn
    end
end

soil_bucket_2 = @hydrobucket :snowstorage begin
    fluxes = begin
        @hydroflux Ps ~ ifelse(T <= Ts, Pn, 0)
        @hydroflux Ms ~ ifelse((T > Tm) && (Hs > 0), as * (T - Tm), 0)
        @hydroflux Fs ~ ifelse((T < Tm) && (Hw > 0), af * as * (Tm - T), 0)
        @hydroflux Gs ~ ifelse(Hs > 0, Gmax, 0)
    end

    dfluxes = begin
        @stateflux Hs ~ Ps + Fs - Ms - Gs
    end
end

soil_bucket_3 = @hydrobucket :waterlayer begin
    fluxes = begin
        @hydroflux Pr ~ ifelse(T > Ts, Pn, 0)
        @hydroflux Mw ~ ifelse(Hw == the * Hs, Pr + Ms, 0)
    end

    dfluxes = begin
        @stateflux Hw ~ Pr + Ms - Fs - Mw
    end
end

soil_bucket_4 = @hydrobucket :soilstorage begin
    fluxes = begin
        @hydroflux Peq ~ Mw + Gs
        @hydroflux RH ~ ifelse(S < Smax, max(Peq - phi, 0), 0)
        @hydroflux Fi ~ Peq - RH
        @hydroflux RD ~ ifelse(S == Smax, Peq, 0)
        @hydroflux Et_pot ~ Ep - Ei
        @hydroflux Et ~ min(max(0, Et_pot * (S - sw) / (sm - sw)), Et_pot)
        @hydroflux L ~ Ksat * S^c
    end

    dfluxes = begin
        @stateflux S ~ Fi - RD - Et - L
    end
end

soil_bucket_5 = @hydrobucket :fastflow begin
    fluxes = begin
        @hydroflux Ls ~ min(L, Lmax)
        @hydroflux Lf ~ L - Ls
        @hydroflux Rf ~ kf * Sfast
    end

    dfluxes = begin
        @stateflux Sfast ~ Lf - Rf
    end
end

soil_bucket_6 = @hydrobucket :flow begin
    fluxes = begin
        @hydroflux Rs ~ ks*Sslow
        @hydroflux Q ~ RH+RD+Rf+Rs
    end

    dfluxes = begin
        @stateflux Sslow ~ Ls-Rs
    end
end
model = @hydromodel :echo begin
    soil_bucket
end

end
