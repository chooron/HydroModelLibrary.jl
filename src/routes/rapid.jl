module rapid

using ..HydroModels
using ..HydroModels: AbstractRoute
using ComponentArrays
using ComponentArrays: getaxes
using LinearAlgebra
"""
    RapidRoute{N} <: AbstractRoute

Implements the RAPID (Routing Application for Parallel computatIon of Discharge) routing model.
"""
struct RapidRoute{N} <: AbstractRoute
    adjacency::AbstractMatrix
    infos::NamedTuple

    function RapidRoute(
        fluxes::Pair{Vector{Num},Vector{Num}};
        network=nothing,
        adjacency::Union{Nothing,AbstractMatrix}=nothing,
        name::Union{Symbol,Nothing}=nothing,
    )
        @parameters rapid_k rapid_x
        inputs, outputs = fluxes[1], fluxes[2]
        @assert length(inputs) == length(outputs) == 1 "The length of inputs and outputs must be 1."
        infos = (; inputs=inputs, outputs=outputs, states=Num[], params=[rapid_k, rapid_x])
        route_name = isnothing(name) ? Symbol("##route#", hash(infos)) : name
        adjacency_mat = if !isnothing(adjacency)
            Matrix(adjacency)
        elseif !isnothing(network)
            Matrix(HydroModels.adjacency_matrix(network))'
        else
            throw(ArgumentError("RapidRoute requires either `adjacency` or `network`."))
        end
        return new{route_name}(adjacency_mat, infos)
    end
end

function (route::RapidRoute)(input::Array, params::ComponentVector; kwargs...)
    ptyidx = get(kwargs, :ptyidx, 1:size(input, 2))
    device = get(kwargs, :device, identity)
    delta_t = get(kwargs, :delta_t, 1.0)
    interp = get(kwargs, :interp, DirectInterpolation)
    solver = get(kwargs, :solver, ManualSolver(mutable=true))
    timeidx = get(kwargs, :timeidx, collect(1:size(input, 3)))
    initstates = zeros(eltype(params), size(input, 2)) |> device

    itpfuncs = interp(input[1, :, :], timeidx)

    expand_params = expand_component_params(params, get_param_names(route), ptyidx)[:params] |> device
    k_ps, x_ps = expand_params[:rapid_k], expand_params[:rapid_x]
    c0 = @. ((delta_t / k_ps) - (2 * x_ps)) / ((2 * (1 - x_ps)) + (delta_t / k_ps))
    c1 = @. ((delta_t / k_ps) + (2 * x_ps)) / ((2 * (1 - x_ps)) + (delta_t / k_ps))
    c2 = @. ((2 * (1 - x_ps)) - (delta_t / k_ps)) / ((2 * (1 - x_ps)) + (delta_t / k_ps))
    new_params = ComponentVector(c0=c0, c1=c1, c2=c2) |> device
    new_params_vec, new_params_axes = Vector(new_params) |> device, getaxes(new_params)
    A = p -> Matrix(I, size(route.adjacency)...) .- diagm(p.c0) * route.adjacency

    function du_func(u, p, t)
        q_out_t1 = u
        q_gen = itpfuncs(t)
        ps = ComponentVector(p, new_params_axes)
        rflux_b = ps.c0 .* q_gen .+ ps.c1 .* (route.adjacency * q_out_t1 .+ q_gen) .+ ps.c2 .* q_out_t1
        A(ps) \ (rflux_b .- A(ps) * u)
    end

    sol_arr = solver(du_func, new_params_vec, initstates, timeidx)
    return reshape(sol_arr, 1, size(sol_arr)...)
end

end

