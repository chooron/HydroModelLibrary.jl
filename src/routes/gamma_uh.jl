module GammaUnitHydro

using ..HydroModels
using ..HydroModelCore
using Symbolics: tosymbol
using NNlib
using ComponentArrays
using SpecialFunctions

struct GammaHydrograph{NT} <: HydroModelCore.AbstractComponent
    name::Symbol
    lenF::Int
    infos::NT

    function GammaHydrograph(
        inputs::Vector{T}, outputs::Vector{T},
        alpha::T, theta::T, lenF::Int=15;
        name::Symbol=:gamma_uh
    ) where {T}
        infos = HydroModelCore.HydroInfos(
            inputs=tosymbol.(inputs),
            outputs=tosymbol.(outputs),
            params=tosymbol.([alpha, theta])
        )
        return new{typeof(infos)}(name, lenF, infos)
    end
end

function _gamma_convolution(uh::GammaHydrograph, x::AbstractArray{T, 2}, pas::ComponentVector) where {T}
    num_grid, time_len = size(x)
    param_names = HydroModels.get_param_names(uh.infos)
    
    α_p = pas.params[param_names[1]]
    θ_p = pas.params[param_names[2]]

    if α_p isa Number
        α_p = fill(α_p, 1, num_grid)
    else
        α_p = reshape(α_p, 1, num_grid)
    end
    if θ_p isa Number
        θ_p = fill(θ_p, 1, num_grid)
    else
        θ_p = reshape(θ_p, 1, num_grid)
    end

    α_eff = α_p .* 60.0 .+ 0.1
    θ_eff = θ_p .* 60.0 .+ 0.5

    t = reshape(T.(0.5:1.0:(uh.lenF-0.5)), uh.lenF, 1)

    log_denom1 = loggamma.(α_eff)
    log_denom2 = α_eff .* log.(θ_eff)
    log_w = (α_eff .- 1) .* log.(t) .- t ./ θ_eff .- (log_denom1 .+ log_denom2)

    w = exp.(log_w)
    w = w ./ sum(w, dims=1)
    
    x_reshaped = reshape(permutedims(x, (2, 1)), time_len, num_grid, 1)
    w_flipped = w[end:-1:1, :]
    w_reshaped = reshape(w_flipped, uh.lenF, 1, num_grid)
    
    padding = (uh.lenF - 1, 0)
    y_conv = conv(x_reshaped, w_reshaped; groups=num_grid, pad=padding)
    
    return dropdims(permutedims(y_conv, (2, 1, 3)), dims=3)
end

(uh::GammaHydrograph)(input::AbstractArray{T, 2}, pas::ComponentVector; kwargs...) where {T} = _gamma_convolution(uh, input, pas)

function (uh::GammaHydrograph)(input::AbstractArray{T, 3}, pas::ComponentVector; kwargs...) where {T}
    y = _gamma_convolution(uh, input[1, :, :], pas)
    return reshape(y, 1, size(y)...)
end

export GammaHydrograph

end