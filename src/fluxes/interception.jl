@variables interception

InterceptionFlux(input::namedtuple(:S, :in), params::namedtuple(:Smax), ::Val{:1}; interception=interception) = begin
    HydroFlux(collect(input) => [interception], collect(params), exprs=[ifelse(input.S ≥ params.Smax, input.in, 0)])
end

InterceptionFlux(input::namedtuple(:in), params::namedtuple(:p1), ::Val{:2}; interception=interception) = begin
    HydroFlux(collect(input) => [interception], collect(params), exprs=[max(0, input.in - params.p1)])
end

InterceptionFlux(input::namedtuple(), params::namedtuple(:p1), ::Val{:3}; interception=interception) = begin
    HydroFlux(collect(input) => [interception], collect(params), exprs=[params.p1])
end

# TODO Interception excess after a time-varying fraction is intercepted
# InterceptionFlux(input::namedtuple(:in), params::namedtuple(:p1, :p2), ::Val{:4}) = begin
#     HydroFlux(collect(input) => [interception], collect(params), exprs=[params.p1+(1-params.p1)*cos(2*pi*params.p2)])
# end

InterceptionFlux(input::namedtuple(:in), params::namedtuple(:p1, :p2), ::Val{:5}; interception=interception) = begin
    HydroFlux(collect(input) => [interception], collect(params), exprs=[params.p1 * input.in - params.p2])
end

export InterceptionFlux