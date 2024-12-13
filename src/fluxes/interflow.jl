@variables interflow

InterflowFlux(input::namedtuple(:S, :in), params::namedtuple(:p1, :Smax), ::Val{:1}; interflow=interflow) = begin
    HydroFlux(collect(input) => [interflow], collect(params),
        exprs=[params.p1 * input.S / params.Smax * input.in]
    )
end

InterflowFlux(input::namedtuple(:S), params::namedtuple(:p1, :p2), ::Val{:2}; interflow=interflow) = begin
    HydroFlux(collect(input) => [interflow], collect(params),
        exprs=[params.p1 * input.S^(1 + params.p2)]
    )
end

InterflowFlux(input::namedtuple(:S), params::namedtuple(:p1, :p2), ::Val{:3}; interflow=interflow) = begin
    HydroFlux(collect(input) => [interflow], collect(params),
        exprs=[params.p1 * input.S^(params.p2)]
    )
end

InterflowFlux(input::namedtuple(:S), params::namedtuple(:p1, :p2), ::Val{:4}; interflow=interflow) = begin
    HydroFlux(collect(input) => [interflow], collect(params),
        exprs=[params.p1 * input.S + params.p2 * input.S^2]
    )
end

InterflowFlux(input::namedtuple(:S), params::namedtuple(:p1), ::Val{:5}; interflow=interflow) = begin
    HydroFlux(collect(input) => [interflow], collect(params),
        exprs=[params.p1 * input.S]
    )
end

InterflowFlux(input::namedtuple(:S1, :S2), params::namedtuple(:p1, :p2, :S2max), ::Val{:6}; interflow=interflow) = begin
    HydroFlux(collect(input) => [interflow], collect(params),
        exprs=[params.p1 * input.S1 * min(0, input.S2 / params.S2max) / (1 - params.p2)]
    )
end

InterflowFlux(input::namedtuple(:S), params::namedtuple(:p1, :p2, :p3, :Smax), ::Val{:7}; interflow=interflow) = begin
    HydroFlux(collect(input) => [interflow], collect(params),
        exprs=[(min(0, input.S - params.p1 * params.Smax) / params.p2)^(1 / params.p3)]
    )
end

InterflowFlux(input::namedtuple(:S), params::namedtuple(:p1, :p2), ::Val{:8}; interflow=interflow) = begin
    HydroFlux(collect(input) => [interflow], collect(params),
        exprs=[params.p1 * min(0, input.S - params.p2)]
    )
end

InterflowFlux(input::namedtuple(:S), params::namedtuple(:p1, :p2, :p3), ::Val{:9}; interflow=interflow) = begin
    HydroFlux(collect(input) => [interflow], collect(params),
        exprs=[params.p1 * min(0, input.S - params.p2)^params.p3]
    )
end

InterflowFlux(input::namedtuple(:S), params::namedtuple(:p1, :p2, :p3), ::Val{:10}; interflow=interflow) = begin
    HydroFlux(collect(input) => [interflow], collect(params),
        exprs=[params.p1 * min(0, input.S - params.p2) / params.p3]
    )
end

InterflowFlux(input::namedtuple(:S), params::namedtuple(:p1, :p2), ::Val{:11}; interflow=interflow) = begin
    HydroFlux(collect(input) => [interflow], collect(params),
        exprs=[ifelse(input.S > params.p2, params.p1, 0)]
    )
end

export InterflowFlux