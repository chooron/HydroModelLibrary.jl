module SpecialProcesses

using SpecialFunctions: gamma

export CONVOL_GR4J_1,
       CONVOL_GR4J_2,
       CONVOL_GAMMA,
       CONVOL_GAMMA2

"""GR4J unit-hydrograph kernel 1 evaluated at travel time `t`."""
function CONVOL_GR4J_1(;
    t::Real,
    x4::Real,
)
    if t <= 0 || x4 <= 0
        return 0.0
    elseif t <= x4
        return 5 / (2 * x4) * (t / x4)^(3 / 2)
    else
        return 0.0
    end
end

"""GR4J unit-hydrograph kernel 2 evaluated at travel time `t`."""
function CONVOL_GR4J_2(;
    t::Real,
    x4::Real,
)
    if t <= 0 || x4 <= 0
        return 0.0
    elseif t <= x4
        return 5 / x4 * (t / x4)^(3 / 2)
    elseif t <= 2 * x4
        return 5 / x4 * (2 - t / x4)^(3 / 2)
    else
        return 0.0
    end
end

"""Gamma unit-hydrograph kernel evaluated at travel time `t`."""
function CONVOL_GAMMA(;
    t::Real,
    a::Real,
    beta::Real,
)
    if t <= 0 || a <= 0 || beta <= 0
        return 0.0
    end
    return (1 / t) * (beta * t)^a / gamma(a) * exp(-beta * t)
end

"""Second gamma transfer-function alias used in Raven process lists."""
function CONVOL_GAMMA2(; kwargs...)
    return CONVOL_GAMMA(; kwargs...)
end

end
