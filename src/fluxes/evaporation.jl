@variables evaporation

EvaporationFlux(input::namedtuple(:S, :pet), params::namedtuple(), ::Val{1}; evaporation=evaporation) = begin
    HydroFlux(collect(input) => [evaporation], collect(params), exprs=[ifelse(input.S > 0, input.pet, 0)])
end

EvaporationFlux(input::namedtuple(:S,), params::namedtuple(:p1, :Smax), ::Val{2}; evaporation=evaporation) = begin
    HydroFlux(collect(input) => [evaporation], collect(params), exprs=[params.p1 * input.S / params.Smax])
end

EvaporationFlux(input::namedtuple(:S, :pet), params::namedtuple(:p1, :Smax), ::Val{3}; evaporation=evaporation) = begin
    HydroFlux(collect(input) => [evaporation], collect(params),
        exprs=[ifelse(input.S < params.p1 * params.Smax, input.pet * input.S / params.p1 / params.Smax, input.pet)])
end

EvaporationFlux(input::namedtuple(:S, :pet), params::namedtuple(:p1, :p2, :Smax), ::Val{4}; evaporation=evaporation) = begin
    HydroFlux(collect(input) => [evaporation], collect(params),
        exprs=[input.pet * max(0, params.p1 * (input.S - params.p2 * params.Smax) / (params.Smax - params.p2 * params.Smax))])
end

EvaporationFlux(input::namedtuple(:S, :pet), params::namedtuple(:p1, :Smax), ::Val{5}; evaporation=evaporation) = begin
    HydroFlux(collect(input) => [evaporation], collect(params),
        exprs=[(1 - params.p1) * input.S / params.Smax * input.pet])
end

EvaporationFlux(input::namedtuple(:S, :pet), params::namedtuple(:p1, :p2, :Smax), ::Val{6}; evaporation=evaporation) = begin
    HydroFlux(collect(input) => [evaporation], collect(params),
        exprs=[min(1.0, input.S / params.p2 * params.Smax) * params.p1 * input.pet])
end

EvaporationFlux(input::namedtuple(:S, :pet), params::namedtuple(:Smax), ::Val{7}; evaporation=evaporation) = begin
    HydroFlux(collect(input) => [evaporation], collect(params),
        exprs=[input.S / params.Smax * input.pet])
end

EvaporationFlux(input::namedtuple(:S1, :S2, :pet), params::namedtuple(:Smax, :p1, :p2), ::Val{8}; evaporation=evaporation) = begin
    HydroFlux(collect(input) => [evaporation], collect(params),
        exprs=[input.S1 / (input.S1 + input.S2) * params.p1 * input.pet * min(1, input.S1 / params.p2)])
end

EvaporationFlux(input::namedtuple(:S1, :S2, :pet), params::namedtuple(:Smax, :p1), ::Val{9}; evaporation=evaporation) = begin
    HydroFlux(collect(input) => [evaporation], collect(params),
        exprs=[input.S1 / (input.S1 + input.S2) * (1 - params.p1) * (input.S1 / (params.Smax - input.S2))])
end

EvaporationFlux(input::namedtuple(:S, :pet), params::namedtuple(:Smax, :p1), ::Val{10}; evaporation=evaporation) = begin
    HydroFlux(collect(input) => [evaporation], collect(params),
        exprs=[input.S / params.Smax * params.p1 * input.pet])
end

EvaporationFlux(input::namedtuple(:S, :pet), params::namedtuple(:Smax), ::Val{11}; evaporation=evaporation) = begin
    HydroFlux(collect(input) => [evaporation], collect(params),
        exprs=[(2 * input.S / params.Smax - (input.S / params.Smax)^2) * input.pet])
end

EvaporationFlux(input::namedtuple(:S, :pet), params::namedtuple(:p1), ::Val{12}; evaporation=evaporation) = begin
    HydroFlux(collect(input) => [evaporation], collect(params),
        exprs=[min(1, exp(2 * (1 - input.S / params.p1))) * input.pet])
end

EvaporationFlux(input::namedtuple(:pet), params::namedtuple(:p1, :p2), ::Val{13}; evaporation=evaporation) = begin
    HydroFlux(collect(input) => [evaporation], collect(params),
        exprs=[(params.p1^params.p2) * input.pet])
end

EvaporationFlux(input::namedtuple(:S, :pet), params::namedtuple(:p1, :p2, :Smin), ::Val{14}; evaporation=evaporation) = begin
    HydroFlux(collect(input) => [evaporation], collect(params),
        exprs=[ifelse(input.S > params.Smin, 0.0, (params.p1^params.p2) * input.pet)])
end

EvaporationFlux(input::namedtuple(:S1, :S2, :pet), params::namedtuple(:p1, :Smax), ::Val{15}; evaporation=evaporation) = begin
    HydroFlux(collect(input) => [evaporation], collect(params),
        exprs=[ifelse(input.S2 < params.p1, input.S1 / params.Smax * input.pet, 0.0)])
end

EvaporationFlux(input::namedtuple(:S, :pet), params::namedtuple(:p1, :p2), ::Val{16}; evaporation=evaporation) = begin
    HydroFlux(collect(input) => [evaporation], collect(params),
        exprs=[ifelse(input.S < params.p2, params.p1 * input.pet, 0.0)])
end

EvaporationFlux(input::namedtuple(:S, :pet), params::namedtuple(:p1), ::Val{17}; evaporation=evaporation) = begin
    HydroFlux(collect(input) => [evaporation], collect(params),
        exprs=[1 / (1 + exp(-params.p1 * input.S)) * input.pet])
end

EvaporationFlux(input::namedtuple(:S, :pet), params::namedtuple(:p1, :p2, :p3), ::Val{18}; evaporation=evaporation) = begin
    HydroFlux(collect(input) => [evaporation], collect(params),
        exprs=[params.p1 * exp(-params.p2 * input.S / params.p3) * input.pet])
end

EvaporationFlux(input::namedtuple(:S, :pet), params::namedtuple(:p1, :p2, :Smax), ::Val{19}; evaporation=evaporation) = begin
    HydroFlux(collect(input) => [evaporation], collect(params),
        exprs=[params.p1 * (input.S / params.Smax)^params.p2 * input.pet])
end

EvaporationFlux(input::namedtuple(:S, :pet), params::namedtuple(:p1, :p2, :Smax), ::Val{20}; evaporation=evaporation) = begin
    HydroFlux(collect(input) => [evaporation], collect(params),
        exprs=[ifelse(input.S < params.p2 * params.Smax, params.p1 * input.S / (params.Smax * params.p2), input.pet)])
end

EvaporationFlux(input::namedtuple(:S, :pet), params::namedtuple(:p1, :p2), ::Val{21}; evaporation=evaporation) = begin
    HydroFlux(collect(input) => [evaporation], collect(params),
        exprs=[ifelse(input.S < params.p1, input.pet,
            ifelse(input.S < params.p2 * params.p1, input.S / params.p1 * input.pet, params.p2 * input.pet))
        ])
end

EvaporationFlux(input::namedtuple(:S, :pet), params::namedtuple(:p1, :p2), ::Val{22}; evaporation=evaporation) = begin
    HydroFlux(collect(input) => [evaporation], collect(params),
        exprs=[ifelse(input.S < params.p1, input.pet,
            ifelse(input.S < params.p2 * params.p1, (input.S - params.p1) / (params.p1 - params.p2) * input.pet, 0))
        ])
end

EvaporationFlux(input::namedtuple(:S, :pet), params::namedtuple(:p1, :p2, :Smax), ::Val{23}; evaporation=evaporation) = begin
    HydroFlux(collect(input) => [evaporation], collect(params),
        exprs=[min(input.S / params.p2 * params.Smax, 1.0) * params.p1 * input.pet + (1 + params.p1) * input.S / params.Smax * input.pet])
end

export EvaporationFlux