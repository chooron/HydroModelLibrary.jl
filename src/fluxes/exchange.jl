@variables exchange

ExchangeFlux(input::namedtuple(:S,), params::namedtuple(:p1, :p2, :p3), ::Val{:1}; exchange=exchange) = begin
    HydroFlux(collect(input) => [exchange], collect(params),
        exprs=[(params.p1 * abs(input.S) + params.p2 * (1 - exp(-params.p3 * abs(input.S)))) * sign(input.S)]
    )
end

ExchangeFlux(input::namedtuple(:S1, :S2), params::namedtuple(:S1max, :S2max), ::Val{:2}; exchange=exchange) = begin
    HydroFlux(collect(input) => [exchange], collect(params),
        exprs=[params.p1 * (input.S1 / params.S1max - input.S2 / params.S2max)]
    )
end

ExchangeFlux(input::namedtuple(:S), params::namedtuple(:p1, :p2), ::Val{:3}; exchange=exchange) = begin
    HydroFlux(collect(input) => [exchange], collect(params),
        exprs=[params.p1 * (input.S - params.p2)]
    )
end

export ExchangeFlux