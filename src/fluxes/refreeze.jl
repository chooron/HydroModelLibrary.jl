@variables refreeze

RefreezeFlux(input::namedtuple(:T), params::namedtuple(:thres, :p1, :p2), ::Val{:1}; refreeze=refreeze) = begin
    HydroFlux(collect(input) => [refreeze], collect(params),
        exprs=[min(0, params.thres - input.T) * params.p1 * params.p2]
    )
end

export RefreezeFlux