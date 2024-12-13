@variables melt

MeltFlux(input::namedtuple(:T), params::namedtuple(:thres, :p1), ::Val{:1}; melt=melt) = begin
    HydroFlux(collect(input) => [melt], collect(params),
        exprs=[min(0, input.T - params.thres) * params.p1]
    )
end

MeltFlux(input::namedtuple(:S, :T), params::namedtuple(:p1), ::Val{:2}; melt=melt) = begin
    HydroFlux(collect(input) => [melt], collect(params),
        exprs=[ifelse(input.S > 0, params.p1, 0)]
    )
end

MeltFlux(input::namedtuple(:S, :T), params::namedtuple(:p1), ::Val{:2}; melt=melt) = begin
    HydroFlux(collect(input) => [melt], collect(params),
        exprs=[ifelse(input.S == 0, min(0, input.T - params.thres) * params.p1, 0)]
    )
end

export MeltFlux