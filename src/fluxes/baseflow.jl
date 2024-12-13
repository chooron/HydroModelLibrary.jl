@variables baseflow

"""
    BaseflowFlux

A collection of baseflow functions for hydrological modeling. Each function represents
different types of reservoir outflow processes.

# Available Models

## Linear Reservoir (`:linear`)
```math
q = p_1S
```

## Non-linear Reservoir (`:nonlinear`)
```math
q = \\left(\\frac{S}{p_1}\\right)^{1/p_2}
```

## Empirical Exponential (`:empirical_exponential`)
```math
q = \\frac{S_{max}^{-4}}{4}S^5
```

## Exponential Store (`:exponential`)
```math
q = p_1e^{-p_2S}
```

## Non-linear Scaled (`:nolinear_scaled`)
```math
q = p_1\\left(\\frac{S}{S_{max}}\\right)^{p_2}
```

## Quadratic Threshold (`:quadratic`)
```math
q = \\begin{cases}
p_1|S|^{p_2} & \\text{if } S > p_2 \\\\
0 & \\text{otherwise}
\\end{cases}
```

## Non-linear Reservoir (`:nolinear_reservoir`)
```math
q = p_1|S|^{p_2}
```

## Exponential Scaled (`:exponential_scaled`)
```math
q = p_1(e^{-p_2S/S_{max}} - 1)
```

## Exceeded Linear (`:exceeded_linear`)
```math
q = \\begin{cases}
p_1(S - p_2) & \\text{if } S > p_2 \\\\
0 & \\text{otherwise}
\\end{cases}
```

# Parameters
- `S`: Storage state variable
- `p1`: Primary parameter (coefficient)
- `p2`: Secondary parameter (exponent or threshold)
- `Smax`: Maximum storage capacity

# Notes
- All functions ensure numerical stability by using absolute values where necessary
- Threshold-based functions use smooth transitions for better numerical behavior
- Storage constraints are handled implicitly in the model implementation

# Examples
```julia
# Linear reservoir with single parameter
flux = BaseflowFlux((S=1.0,), (p1=0.1,), Val(:linear))

# Non-linear reservoir with two parameters
flux = BaseflowFlux((S=1.0,), (p1=0.1, p2=2.0), Val(:nonlinear))

# Scaled reservoir with maximum storage
flux = BaseflowFlux((S=1.0,), (p1=0.1, p2=2.0, Smax=100.0), Val(:nolinear_scaled))
```
"""
BaseflowFlux(input::namedtuple(:S,), params::namedtuple(:p1,), ::Val{1}; baseflow=baseflow) = begin
    HydroFlux(collect(input) => [baseflow], collect(params), exprs=[params.p1 * input.S])
end

BaseflowFlux(input::namedtuple(:S,), params::namedtuple(:p1, :p2), ::Val{2}; baseflow=baseflow) = begin
    HydroFlux(collect(input) => [baseflow], collect(params),
        exprs=[abs(1 / params.p1 * input.S)^(1 / params.p2)]
    )
end

BaseflowFlux(input::namedtuple(:S,), params::namedtuple(:Smax,), ::Val{3}; baseflow=baseflow) = begin
    HydroFlux(collect(input) => [baseflow], collect(params),
        exprs=[abs(params.Smax)^(-4) / (4) * (abs(input.S)^5)]
    )
end

BaseflowFlux(input::namedtuple(:S,), params::namedtuple(:p1, :p2), ::Val{4}; baseflow=baseflow) = begin
    HydroFlux(collect(input) => [baseflow], collect(params),
        exprs=[params.p1 * exp(-params.p2 * input.S)]
    )
end

BaseflowFlux(input::namedtuple(:S,), params::namedtuple(:p1, :p2, :Smax), ::Val{5}; baseflow=baseflow) = begin
    HydroFlux(collect(input) => [baseflow], collect(params),
        exprs=[params.p1 * (abs(input.S / params.Smax)^params.p2)]
    )
end

BaseflowFlux(input::namedtuple(:S,), params::namedtuple(:p1, :p2), ::Val{6}; baseflow=baseflow) = begin
    HydroFlux(collect(input) => [baseflow], collect(params),
        exprs=[ifelse(input.S > params.p2, params.p1 * abs(input.S)^params.p2, 0)]
    )
end

BaseflowFlux(input::namedtuple(:S,), params::namedtuple(:p1, :p2), ::Val{7}; baseflow=baseflow) = begin
    HydroFlux(collect(input) => [baseflow], collect(params),
        exprs=[params.p1 * abs(input.S)^params.p2]
    )
end

BaseflowFlux(input::namedtuple(:S,), params::namedtuple(:p1, :p2, :Smax), ::Val{8}; baseflow=baseflow) = begin
    HydroFlux(collect(input) => [baseflow], collect(params),
        exprs=[params.p1 * (exp(-params.p2 * input.S / params.Smax) - 1)]
    )
end

BaseflowFlux(input::namedtuple(:S,), params::namedtuple(:p1, :p2), ::Val{9}; baseflow=baseflow) = begin
    HydroFlux(collect(input) => [baseflow], collect(params),
        exprs=[ifelse(input.S > params.p2, params.p1 * (input.S - params.p2), 0)]
    )
end

export BaseflowFlux