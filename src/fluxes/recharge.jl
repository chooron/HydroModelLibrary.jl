@variables recharge

RechargeFlux(input::namedtuple(:S, :in), params::namedtuple(:p1, :Smax), ::Val{:1}; recharge=recharge) = begin
    HydroFlux(collect(input) => [recharge], collect(params),
        exprs=[params.p1 * input.S / params.Smax * input.in]
    )
end

RechargeFlux(input::namedtuple(:S, :in), params::namedtuple(:p1, :Smax), ::Val{:2}; recharge=recharge) = begin
    HydroFlux(collect(input) => [recharge], collect(params),
        exprs=[((input.S / params.Smax)^params.p1) * input.in]
    )
end

RechargeFlux(input::namedtuple(:S), params::namedtuple(:p1), ::Val{:3}; recharge=recharge) = begin
    HydroFlux(collect(input) => [recharge], collect(params),
        exprs=[params.p1 * input.S]
    )
end

RechargeFlux(input::namedtuple(:S), params::namedtuple(:p1), ::Val{:4}; recharge=recharge) = begin
    HydroFlux(collect(input) => [recharge], collect(params),
        exprs=[ifelse(input.S > 0, params.p1, 0)]
    )
end

RechargeFlux(input::namedtuple(:S1, :S2), params::namedtuple(:p1, :p2), ::Val{:5}; recharge=recharge) = begin
    HydroFlux(collect(input) => [recharge], collect(params),
        exprs=[params.p1 * input.S1 * max(0, 1 - input.S2 / params.p2)]
    )
end

RechargeFlux(input::namedtuple(:S), params::namedtuple(:p1, :p2), ::Val{:6}; recharge=recharge) = begin
    HydroFlux(collect(input) => [recharge], collect(params),
        exprs=[params.p1 * input.S^params.p2]
    )
end

RechargeFlux(input::namedtuple(), params::namedtuple(:p1), ::Val{:7}; recharge=recharge) = begin
    HydroFlux(collect(input) => [recharge], collect(params),
        exprs=[params.p1]
    )
end

export RechargeFlux