@variables percolation

PercolationFlux(input::namedtuple(:S), params::namedtuple(:p1), ::Val{:1}; percolation=percolation) = begin
    HydroFlux(collect(input) => [percolation], collect(params),
        exprs=[ifelse(input.S > 0, params.p1, 0)]
    )
end

PercolationFlux(input::namedtuple(:S), params::namedtuple(:p1, :Smax), ::Val{:2}; percolation=percolation) = begin
    HydroFlux(collect(input) => [percolation], collect(params),
        exprs=[params.p1 * input.S / params.Smax]
    )
end

PercolationFlux(input::namedtuple(:S), params::namedtuple(:p1, :Smax), ::Val{:3}; percolation=percolation) = begin
    HydroFlux(collect(input) => [percolation], collect(params),
        exprs=[params.Smax^(-4) / 4 * (4 / 9)^(-4) * input.S^5]
    )
end

# TODO Demand-based percolation scaled by available moisture 

PercolationFlux(input::namedtuple(:S), params::namedtuple(:p1, :p2, :Smax), ::Val{:5}; percolation=percolation) = begin
    HydroFlux(collect(input) => [percolation], collect(params),
        exprs=[params.p1 * (input.S / params.Smax)^params.p2]
    )
end

PercolationFlux(input::namedtuple(:S), params::namedtuple(:p1, :p2), ::Val{:6}; percolation=percolation) = begin
    HydroFlux(collect(input) => [percolation], collect(params),
        exprs=[min(1, max(0, input.S) / params.p2) * params.p1]
    )
end

export PercolationFlux