@variables rainfall

RainfallFlux(input::namedtuple(:P, :T), params::namedtuple(:thres), ::Val{:1}; rainfall=rainfall) = begin
    HydroFlux(collect(input) => [rainfall], collect(params),
        exprs=[ifelse(input.T > params.thres, input.P, 0)]
    )
end

RainfallFlux(input::namedtuple(:P, :T), params::namedtuple(:p1, :p2), ::Val{:2}; rainfall=rainfall) = begin
    HydroFlux(collect(input) => [rainfall], collect(params),
        exprs=[ifelse(input.T ≤ params.p1 - 1 / 2 * params.p2, 0,
            ifelse(input.T ≤ params.p1 + 1 / 2 * params.p2, input.P * (params.p1 + 1 / 2 * params.p2 - input.T) / params.p2, input.P))]
    )
end

export RainfallFlux