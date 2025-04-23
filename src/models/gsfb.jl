
module gsfb
# Model variables
@variables S [description = "current storage in the upper zone", unit = "mm"]
@variables P [description = "precipitation", unit = "mm/d"]
@variables Qdr [description = "deep groundwater", unit = "mm/d"]
@variables Ea [description = "evaporation", unit = "mm/d"]
@variables Qs [description = "surface runoff", unit = "mm/d"]
@variables F [description = "infiltration", unit = "mm/d"]
@variables Ep [description = "potential rate", unit = "mm/d"]
@variables DS [description = "current deep storage", unit = "mm"]
@variables SS [description = "current storage in the subsurface store", unit = "mm"]
@variables Qb [description = "baseflow", unit = "mm/d"]
@variables Dp [description = "deep percolation", unit = "mm/d"]
@variables Qt [description = "Total flow", unit = "mm/d"]

# Model parameters
@parameters C [description = "Recharge coefficient", bounds = (0, 1), unit = "d-1"]
@parameters NDC [description = "Recharge coefficient", bounds = (0.05, 0.95), unit = "-"]
@parameters Smax [description = "Maximum soil moisture storage", bounds = (1, 2000), unit = "mm"]
@parameters Emax [description = "Maximum evaporation rate", bounds = (0,20), unit = "mm/d"]
@parameters Frate [description = "Recharge rate", bounds = (0,200), unit = "mm/d"]
@parameters B [description = "Fraction subsurface flow to stream", bounds = (0,1), unit = "-"]
@parameters DPF [description = "Runoff coefficient", bounds = (0,1), unit = "d-1"]
@parameters SDRmax [description = "Threshold for subsurface flow generation", bounds = (1,300), unit = "mm"]

# Soil water component
soil_bucket_1 = @hydrobucket :soil begin
    fluxes = begin
        @hydroflux Qdr ~ C*DS*(1-min(S/(NDC*Smax),1.00))
        @hydroflux Ea ~ min(Ep,Emax*min(S/(NDC*Smax),1.00))
        @hydroflux Qs ~ S == Smax ? P : 0.0
        @hydroflux F ~ Frate*step_func(S-NDC*Smax)
        
    end

    dfluxes = begin
        @stateflux S ~ P+Qdr-Ea-Qs-F
    end
end

soil_bucket_2 = @hydrobucket :soil begin
    fluxes = begin
        @hydroflux Qb ~ B*DPF*step_func(SS-SDRmax)
        @hydroflux Dp ~ (1-B)*DPF*SS
        
    end

    dfluxes = begin
        @stateflux SS ~ F-Qb-Dp
    end
end

soil_bucket_3 = @hydrobucket :soil begin
    fluxes = begin
        @hydroflux Qt ~ Qs+Qb
        
    end

    dfluxes = begin
        @stateflux DS ~ Dp-Qdr
    end
end

model = @hydromodel :gsfb begin
    soil_bucket
end

end