@variables snowfall

SnowfallFlux(input::namedtuple(:P, :T), params::namedtuple(:thres), ::Val{:1}; snowfall=snowfall) = begin
    HydroFlux(collect(input) => [snowfall], collect(params),
        exprs=[ifelse(input.T > params.thres, 0, input.P)]
    )
end

SnowfallFlux(input::namedtuple(:P, :T), params::namedtuple(:p1, :p2), ::Val{:2}; snowfall=snowfall) = begin
    HydroFlux(collect(input) => [snowfall], collect(params),
        exprs=[ifelse(input.T ≤ params.p1 - 1 / 2 * params.p2, input.P,
            ifelse(input.T ≤ params.p1 + 1 / 2 * params.p2, input.P * (params.p1 + 1 / 2 * params.p2 - input.T) / params.p2, 0))]
    )
end

export SnowfallFlux