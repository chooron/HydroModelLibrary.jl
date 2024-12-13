@variables saturation

SaturationFlux(input::namedtuple(:S, :in), params::namedtuple(:Smax), ::Val{:1}; saturation=saturation) = begin
    HydroFlux(collect(input) => [saturation], collect(params),
        exprs=[ifelse(input.S > params.Smax, input.in, 0)]
    )
end

SaturationFlux(input::namedtuple(:in), params::namedtuple(:p1, :Smax), ::Val{:2}; saturation=saturation) = begin
    HydroFlux(collect(input) => [saturation], collect(params),
        exprs=[(1 - min(1, (1 - input.S / params.Smax)^params.p1)) * input.in]
    )
end

SaturationFlux(input::namedtuple(:S, :in), params::namedtuple(:p1, :Smax), ::Val{:3}; saturation=saturation) = begin
    HydroFlux(collect(input) => [saturation], collect(params),
        exprs=[(1 - min(1, 1 / (1 + exp((input.S / params.Smax + 0.5) / params.p1)))) * input.in]
    )
end

SaturationFlux(input::namedtuple(:S, :in), params::namedtuple(:p1, :Smax), ::Val{:4}; saturation=saturation) = begin
    HydroFlux(collect(input) => [saturation], collect(params),
        exprs=[(1 - min(1, input.S / params.Smax)^2) * input.in]
    )
end

SaturationFlux(input::namedtuple(:S, :in), params::namedtuple(:p1, :p2), ::Val{:5}; saturation=saturation) = begin
    HydroFlux(collect(input) => [saturation], collect(params),
        exprs=[(1 - min(1, (input.S / params.p1)^params.p2)) * input.in]
    )
end

SaturationFlux(input::namedtuple(:S, :in), params::namedtuple(:p1, :Smax), ::Val{:6}; saturation=saturation) = begin
    HydroFlux(collect(input) => [saturation], collect(params),
        exprs=[params.p1 * input.S / params.Smax * input.in]
    )
end

#TODO Saturation excess from a store with different degrees of saturation (gamma function variant) 

SaturationFlux(input::namedtuple(:S, :in), params::namedtuple(:p1, :p2, :Smax), ::Val{:8}; saturation=saturation) = begin
    HydroFlux(collect(input) => [saturation], collect(params),
        exprs=[(params.p1 + max(0, params.p2 - params.p1) * input.S / params.Smax) * input.in]
    )
end

SaturationFlux(input::namedtuple(:S, :in), params::namedtuple(), ::Val{:9}; saturation=saturation) = begin
    HydroFlux(collect(input) => [saturation], collect(params),
        exprs=[ifelse(input.S == 0, in, 0)]
    )
end

SaturationFlux(input::namedtuple(:S, :in), params::namedtuple(:p1, :p2, :p3), ::Val{:10}; saturation=saturation) = begin
    HydroFlux(collect(input) => [saturation], collect(params),
        exprs=[input.in * min(params.p1, params.p2 + params.p2 * exp(params.p3 * input.S))]
    )
end

SaturationFlux(input::namedtuple(:S, :in), params::namedtuple(:p1, :p2, :Smax, :Smin), ::Val{:11}; saturation=saturation) = begin
    HydroFlux(collect(input) => [saturation], collect(params),
        exprs=[params.p1 * clamp((input.S - params.Smin) / (params.Smax - params.Smin), 0, 1)^params.p2 * input.in]
    )
end

SaturationFlux(input::namedtuple(:in), params::namedtuple(:p1, :p2), ::Val{:12}; saturation=saturation) = begin
    HydroFlux(collect(input) => [saturation], collect(params),
        exprs=[(params.p1 - params.p2) / (1 - params.p2) * input.in]
    )
end

# TODO Saturation excess flow from a store with different degrees of saturation (normal distribution variant) 

SaturationFlux(input::namedtuple(:in), params::namedtuple(:p1, :p2), ::Val{:14}; saturation=saturation) = begin
    HydroFlux(collect(input) => [saturation], collect(params),
        exprs=[ifelse(input.S / params.Smax ≤ 0.5 - params.p1,
            (0.5 - params.p1)^(1 - params.p2) * max(0, input.S / params.Smax)^params.p3,
            1 - (0.5 - params.p1)^(1 - params.p2) * max(0, 1 - input.S / params.Smax)^params.p3) * input.in]
    )
end

export SaturationFlux