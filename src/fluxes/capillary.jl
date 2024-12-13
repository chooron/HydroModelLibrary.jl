"""
    CapillaryFlux

A collection of capillary flux functions for hydrological modeling.

# Available Models

## Scaled Linear (`:scaled`)
```math
q = p_1\\left(1 - \\frac{S}{S_{max}}\\right)
```
Linear decrease of capillary flux with relative soil moisture.
- `p1`: Maximum capillary flux rate [L/T]
- `Smax`: Maximum soil moisture storage [L]

## Constant (`:constant`)
```math
q = \\begin{cases}
p_1 & \\text{if } S > 0 \\\\
0 & \\text{otherwise}
\\end{cases}
```
Constant capillary flux when soil moisture is available.
- `p1`: Constant capillary flux rate [L/T]

## Threshold (`:threshold`)
```math
q = \\begin{cases}
p_1\\left(1 - \\frac{S}{p_2}\\right) & \\text{if } S < p_2 \\\\
0 & \\text{otherwise}
\\end{cases}
```
Linear decrease with threshold cutoff.
- `p1`: Maximum capillary flux rate [L/T]
- `p2`: Threshold soil moisture storage [L]
"""

@variables capillary

CapillaryFlux(input::namedtuple(:S,), params::namedtuple(:p1,:Smax), ::Val{:scaled}; capillary=capillary) = begin
    HydroFlux(collect(input) => [capillary], collect(params), exprs=[params.p1 * (1 - input.S / params.Smax)])
end

CapillaryFlux(input::namedtuple(:S,), params::namedtuple(:p1,), ::Val{:constant}; capillary=capillary) = begin
    HydroFlux(collect(input) => [capillary], collect(params), exprs=[ifelse(input.S > 0, params.p1, 0)])
end

CapillaryFlux(input::namedtuple(:S,), params::namedtuple(:p1,:p2), ::Val{:threshold}; capillary=capillary) = begin
    HydroFlux(collect(input) => [capillary], collect(params), exprs=[ifelse(input.S > params.p2, 0, params.p1 * (1 - input.S / params.p2))])
end

export CapillaryFlux