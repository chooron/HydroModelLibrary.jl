@variables infiltration

InfiltrationFlux(input::namedtuple(:S), params::namedtuple(:p1, :p2, :Smax), ::Val{:1}; infiltration=infiltration) = begin
    HydroFlux(collect(input) => [infiltration], collect(params),
        exprs=[params.p1 * exp(-params.p2 * input.S / params.Smax)]
    )
end

InfiltrationFlux(input::namedtuple(:S, :used), params::namedtuple(:p1, :p2, :Smax), ::Val{:2}; infiltration=infiltration) = begin
    HydroFlux(collect(input) => [infiltration], collect(params),
        exprs=[params.p1 * exp(-params.p2 * input.S / params.Smax) - input.used]
    )
end

InfiltrationFlux(input::namedtuple(:S, :in), params::namedtuple(:Smax), ::Val{:3}; infiltration=infiltration) = begin
    HydroFlux(collect(input) => [infiltration], collect(params), exprs=[ifelse(input.S > params.Smax, input.in, 0)])
end

InfiltrationFlux(input::namedtuple(), params::namedtuple(:p1), ::Val{:4}) = begin
    HydroFlux(collect(input) => [infiltration], collect(params), exprs=[params.p1])
end

InfiltrationFlux(input::namedtuple(:S1, :S2), params::namedtuple(:S1max, :S2max, :p1, :p2), ::Val{:5}; infiltration=infiltration) = begin
    HydroFlux(collect(input) => [infiltration], collect(params),
        exprs=[params.p1 * (1 - input.S1 / params.S1max) * abs(input.S2 / params.S2max)^(-params.p2)]
    )
end

InfiltrationFlux(input::namedtuple(:S, :in), params::namedtuple(:Smax, :p1, :p2), ::Val{:6}; infiltration=infiltration) = begin
    HydroFlux(collect(input) => [infiltration], collect(params),
        exprs=[params.p1 * abs(input.S / params.Smax)^params.p2 * input.in]
    )
end

InfiltrationFlux(input::namedtuple(:S), params::namedtuple(:Smax, :p1, :p2), ::Val{:7}; infiltration=infiltration) = begin
    HydroFlux(collect(input) => [infiltration], collect(params),
        exprs=[ifelse(input.S ≥ params.Smax, params.p1 * exp(-params.p2 * input.S / params.Smax), 0)]
    )
end

export InfiltrationFlux