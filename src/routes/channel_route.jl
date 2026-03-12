module channel_route

using ..HydroModels
using ..HydroModels:
    ChannelRoute,
    RouteIRF,
    SparseRouteKernel,
    SparseRouteConvolution,
    build_irf_kernels,
    aggregate_route_kernel
using ..unithydro:
    ROUTE_DIFFUSIVE_WAVE,
    ROUTE_NONE,
    ROUTE_PLUG_FLOW,
    ROUTE_LINEAR_STORAGE,
    ROUTE_STORAGE_COEFF,
    ROUTE_MUSKINGUM,
    ROUTE_NASH_CASCADE,
    ROUTE_HYDROLOGIC
using ComponentArrays
using SpecialFunctions: loggamma

const _MIN_POSITIVE = 1.0e-6
const _CONVOLUTION_ALGORITHMS = (
    ROUTE_PLUG_FLOW,
    ROUTE_DIFFUSIVE_WAVE,
    ROUTE_LINEAR_STORAGE,
    ROUTE_NASH_CASCADE,
)
const _RECURSIVE_ALGORITHMS = (
    ROUTE_STORAGE_COEFF,
    ROUTE_MUSKINGUM,
    ROUTE_HYDROLOGIC,
)

@inline _is_convolution_algorithm(algorithm::Symbol) = algorithm in _CONVOLUTION_ALGORITHMS
@inline _is_recursive_algorithm(algorithm::Symbol) = algorithm in _RECURSIVE_ALGORITHMS
@inline _channel_route_max_lag(algorithm::Symbol, max_lag::Integer) = _is_convolution_algorithm(algorithm) ? Int(max_lag) : 1

@inline _storage_coeffs(K, delta_t) = begin
    k = min(1 / (max(K, _MIN_POSITIVE) / delta_t + 0.5), 1.0)
    return k / 2, k / 2, 1 - k
end

@inline _muskingum_coeffs(K, X, delta_t) = begin
    denom = 2 * max(K, _MIN_POSITIVE) * (1 - X) + delta_t
    return (
        (delta_t - 2 * K * X) / denom,
        (delta_t + 2 * K * X) / denom,
        (2 * K * (1 - X) - delta_t) / denom,
    )
end

function _plug_kernel(flow_length, wave_celerity, delta_t, max_lag)
    lag = max(flow_length / max(wave_celerity, _MIN_POSITIVE) / delta_t, 0.0)
    weights = zeros(Float64, max(1, max_lag))
    base = floor(Int, lag)
    frac = lag - base
    first_idx = min(base + 1, length(weights))
    weights[first_idx] += 1 - frac
    if frac > 0 && first_idx + 1 <= length(weights)
        weights[first_idx + 1] += frac
    end
    sum(weights) > 0 || (weights[1] = 1.0)
    return weights / sum(weights)
end

function _diffusive_kernel(flow_length, wave_celerity, diffusivity, delta_t, max_lag)
    weights = zeros(Float64, max(1, max_lag))
    c_ref = max(wave_celerity, _MIN_POSITIVE)
    diffusivity_eff = max(diffusivity, _MIN_POSITIVE)

    for lag in 1:length(weights)
        t = max((lag - 0.5) * delta_t, _MIN_POSITIVE)
        weights[lag] = (1 / (2 * sqrt(pi * diffusivity_eff * t))) *
                       exp(-((flow_length - c_ref * t)^2) / (4 * diffusivity_eff * t))
    end

    if sum(weights) <= 0
        return _plug_kernel(flow_length, c_ref, delta_t, max_lag)
    end
    return weights / sum(weights)
end

function _linear_storage_kernel(storage_time, delta_t, max_lag)
    lag_count = max(1, max_lag)
    tau = max(storage_time / max(delta_t, _MIN_POSITIVE), _MIN_POSITIVE)
    x = 1 / (1 + tau)
    weights = zeros(Float64, lag_count)

    @inbounds for lag in 1:lag_count
        weights[lag] = x * (1 - x)^(lag - 1)
    end

    return weights
end

function _nash_cascade_kernel(storage_time, n_reservoirs, delta_t, max_lag)
    lag_count = max(1, max_lag)
    tau = max(storage_time / max(delta_t, _MIN_POSITIVE), _MIN_POSITIVE)
    n = max(Float64(n_reservoirs), 1.0)
    x = 1 / (1 + tau)
    weights = zeros(Float64, lag_count)
    log_x = log(max(x, _MIN_POSITIVE))
    log_one_minus_x = log(max(1 - x, _MIN_POSITIVE))

    @inbounds for lag in 0:(lag_count - 1)
        log_binom = loggamma(lag + n) - loggamma(lag + 1) - loggamma(n)
        weights[lag + 1] = exp(log_binom + n * log_x + lag * log_one_minus_x)
    end

    return weights
end

function _hydrologic_step(q_prev, q_in_prev, q_in_curr, storage_coeff, storage_exponent, delta_t; iterations=12)
    rhs = 0.5 * (q_in_prev + q_in_curr)
    q_new = max(q_prev, rhs, 0.0)
    k = max(storage_coeff, _MIN_POSITIVE)
    m = max(storage_exponent, _MIN_POSITIVE)

    for _ in 1:iterations
        q_eval = max(q_new, 0.0)
        lhs = k * (q_eval^m - max(q_prev, 0.0)^m) / delta_t + 0.5 * (q_prev + q_eval)
        resid = lhs - rhs
        deriv = k * m * max(q_eval, _MIN_POSITIVE)^(m - 1) / delta_t + 0.5
        q_new = max(q_eval - resid / max(deriv, _MIN_POSITIVE), 0.0)
    end

    return q_new
end

function _build_route_irf(::Val{ROUTE_PLUG_FLOW}, route::ChannelRoute)
    return RouteIRF(
        [:flow_length, :wave_celerity],
        (p, dt, horizon) -> _plug_kernel(p[1], p[2], dt, min(route.max_lag, horizon));
        name=route.name,
    )
end

function _build_route_irf(::Val{ROUTE_DIFFUSIVE_WAVE}, route::ChannelRoute)
    return RouteIRF(
        [:flow_length, :wave_celerity, :diffusivity],
        (p, dt, horizon) -> _diffusive_kernel(p[1], p[2], p[3], dt, min(route.max_lag, horizon));
        name=route.name,
    )
end

function _build_route_irf(::Val{ROUTE_LINEAR_STORAGE}, route::ChannelRoute)
    return RouteIRF(
        [:storage_time],
        (p, dt, horizon) -> _linear_storage_kernel(p[1], dt, min(route.max_lag, horizon));
        name=route.name,
    )
end

function _build_route_irf(::Val{ROUTE_NASH_CASCADE}, route::ChannelRoute)
    return RouteIRF(
        [:storage_time, :n_reservoirs],
        (p, dt, horizon) -> _nash_cascade_kernel(p[1], p[2], dt, min(route.max_lag, horizon));
        name=route.name,
    )
end

function _build_route_irf(::Val{algorithm}, route::ChannelRoute) where {algorithm}
    throw(ArgumentError("No RouteIRF available for $(route.algorithm)."))
end

_build_route_irf(route::ChannelRoute) = _build_route_irf(Val(route.algorithm), route)

function _build_sparse_route(route::ChannelRoute, kernel::SparseRouteKernel)
    return SparseRouteConvolution(
        kernel;
        name=route.name,
        inputs=HydroModels.get_input_names(route),
        outputs=HydroModels.get_output_names(route),
        params=HydroModels.get_param_names(route),
    )
end

_simulate_single_none(route::ChannelRoute, inflow::AbstractVector, params::AbstractVector, delta_t::Float64) = Float64.(inflow)

function _simulate_single_convolution(route::ChannelRoute, inflow::AbstractVector, params::AbstractVector, delta_t::Float64)
    horizon = max(1, min(route.max_lag, length(inflow)))
    irf = _build_route_irf(route)
    weights = irf(params; delta_t=delta_t, horizon=horizon)
    kernel = SparseRouteKernel(reshape(Int[1, 1], 1, 2), reshape(weights, 1, :), 1, 1)
    return _build_sparse_route(route, kernel)(Float64.(inflow))
end

function _simulate_single_recursive(route::ChannelRoute, inflow::AbstractVector, params::AbstractVector, delta_t::Float64)
    out = zeros(Float64, length(inflow))
    q_prev = 0.0
    q_in_prev = 0.0

    for tidx in eachindex(inflow)
        q_in_curr = Float64(inflow[tidx])
        if route.algorithm == ROUTE_STORAGE_COEFF
            c1, c2, c3 = _storage_coeffs(params[1], delta_t)
            out[tidx] = max(c1 * q_in_curr + c2 * q_in_prev + c3 * q_prev, 0.0)
        elseif route.algorithm == ROUTE_MUSKINGUM
            c1, c2, c3 = _muskingum_coeffs(params[1], params[2], delta_t)
            out[tidx] = max(c1 * q_in_curr + c2 * q_in_prev + c3 * q_prev, 0.0)
        elseif route.algorithm == ROUTE_HYDROLOGIC
            out[tidx] = _hydrologic_step(q_prev, q_in_prev, q_in_curr, params[1], params[2], delta_t)
        else
            throw(ArgumentError("Unsupported recursive channel routing algorithm: $(route.algorithm)"))
        end
        q_prev = out[tidx]
        q_in_prev = q_in_curr
    end

    return out
end

function _upstream_sum(adjacency, routed, node)
    isnothing(adjacency) && return 0.0
    acc = 0.0
    @inbounds for upstream in 1:length(routed)
        acc += adjacency[node, upstream] * routed[upstream]
    end
    return acc
end

@inline _node_route_params(params::AbstractVector, node::Int) = [param[node] for param in params]

function _simulate_network_none(route::ChannelRoute, input::AbstractMatrix, params::AbstractVector, delta_t::Float64)
    nn, nt = size(input)
    order = isnothing(route.topo_order) ? collect(1:nn) : route.topo_order
    routed = zeros(Float64, nn, nt)

    for tidx in 1:nt
        curr_out = zeros(Float64, nn)
        for node in order
            upstream_curr = _upstream_sum(route.adjacency, curr_out, node)
            curr_out[node] = upstream_curr + Float64(input[node, tidx])
        end
        routed[:, tidx] = curr_out
    end

    return routed
end

function _simulate_network_convolution(route::ChannelRoute, input::AbstractMatrix, params::AbstractVector, delta_t::Float64)
    nn, nt = size(input)
    irf = _build_route_irf(route)
    node_kernels = build_irf_kernels(irf, params; delta_t=delta_t, horizon=min(route.max_lag, nt))
    kernel = aggregate_route_kernel(route.adjacency, node_kernels; topo_order=route.topo_order, horizon=nt)
    return _build_sparse_route(route, kernel)(Float64.(input))
end

function _simulate_network_recursive(route::ChannelRoute, input::AbstractMatrix, params::AbstractVector, delta_t::Float64)
    nn, nt = size(input)
    order = isnothing(route.topo_order) ? collect(1:nn) : route.topo_order
    routed = zeros(Float64, nn, nt)
    prev_out = zeros(Float64, nn)
    prev_local = zeros(Float64, nn)

    for tidx in 1:nt
        curr_out = zeros(Float64, nn)
        for node in order
            upstream_curr = _upstream_sum(route.adjacency, curr_out, node)
            upstream_prev = _upstream_sum(route.adjacency, prev_out, node)
            local_curr = Float64(input[node, tidx])
            node_params = _node_route_params(params, node)

            if route.algorithm == ROUTE_STORAGE_COEFF
                c1, c2, c3 = _storage_coeffs(node_params[1], delta_t)
                curr_out[node] = max(c1 * upstream_curr + c2 * upstream_prev + c3 * prev_out[node] + local_curr, 0.0)
            elseif route.algorithm == ROUTE_MUSKINGUM
                c1, c2, c3 = _muskingum_coeffs(node_params[1], node_params[2], delta_t)
                curr_out[node] = max(c1 * upstream_curr + c2 * upstream_prev + c3 * prev_out[node] + local_curr, 0.0)
            elseif route.algorithm == ROUTE_HYDROLOGIC
                total_prev = upstream_prev + prev_local[node]
                total_curr = upstream_curr + local_curr
                curr_out[node] = _hydrologic_step(prev_out[node], total_prev, total_curr, node_params[1], node_params[2], delta_t)
            else
                throw(ArgumentError("Unsupported recursive channel routing algorithm: $(route.algorithm)"))
            end
        end
        routed[:, tidx] = curr_out
        prev_out = curr_out
        prev_local = Float64.(input[:, tidx])
    end

    return routed
end

function _channel_route_handlers(algorithm::Symbol)
    if algorithm == ROUTE_NONE
        return _simulate_single_none, _simulate_network_none
    elseif _is_convolution_algorithm(algorithm)
        return _simulate_single_convolution, _simulate_network_convolution
    elseif _is_recursive_algorithm(algorithm)
        return _simulate_single_recursive, _simulate_network_recursive
    end

    throw(ArgumentError("Unsupported channel routing algorithm: $algorithm"))
end

function _build_named_channel_route(
    algorithm::Symbol,
    params::AbstractVector;
    input,
    output,
    adjacency=nothing,
    htypes=nothing,
    max_lag::Integer=32,
    name::Union{Nothing,Symbol}=nothing,
)
    simulate_single, simulate_network = _channel_route_handlers(algorithm)
    return @channelroute begin
        name = name
        algorithm = algorithm
        simulate_single = simulate_single
        simulate_network = simulate_network
        inputs = [input]
        outputs = [output]
        params = params
        adjacency = adjacency
        max_lag = _channel_route_max_lag(algorithm, max_lag)
        htypes = htypes
    end
end

_required_route_kw(kwargs, key::Symbol, algorithm::Symbol) = get(kwargs, key, nothing)

function _required_route_param(kwargs, key::Symbol, algorithm::Symbol)
    value = _required_route_kw(kwargs, key, algorithm)
    isnothing(value) && throw(ArgumentError("$(algorithm) requires `$(key)`."))
    return value
end

_channel_route_params(::Val{ROUTE_NONE}; kwargs...) = HydroModels.Num[]
_channel_route_params(::Val{ROUTE_PLUG_FLOW}; kwargs...) = [
    _required_route_param(kwargs, :flow_length, ROUTE_PLUG_FLOW),
    _required_route_param(kwargs, :wave_celerity, ROUTE_PLUG_FLOW),
]
_channel_route_params(::Val{ROUTE_DIFFUSIVE_WAVE}; kwargs...) = [
    _required_route_param(kwargs, :flow_length, ROUTE_DIFFUSIVE_WAVE),
    _required_route_param(kwargs, :wave_celerity, ROUTE_DIFFUSIVE_WAVE),
    _required_route_param(kwargs, :diffusivity, ROUTE_DIFFUSIVE_WAVE),
]
_channel_route_params(::Val{ROUTE_LINEAR_STORAGE}; kwargs...) = [
    _required_route_param(kwargs, :storage_time, ROUTE_LINEAR_STORAGE),
]
_channel_route_params(::Val{ROUTE_STORAGE_COEFF}; kwargs...) = [
    _required_route_param(kwargs, :storage_time, ROUTE_STORAGE_COEFF),
]
_channel_route_params(::Val{ROUTE_MUSKINGUM}; kwargs...) = [
    _required_route_param(kwargs, :muskingum_k, ROUTE_MUSKINGUM),
    _required_route_param(kwargs, :muskingum_x, ROUTE_MUSKINGUM),
]
_channel_route_params(::Val{ROUTE_NASH_CASCADE}; kwargs...) = [
    _required_route_param(kwargs, :storage_time, ROUTE_NASH_CASCADE),
    _required_route_param(kwargs, :n_reservoirs, ROUTE_NASH_CASCADE),
]
_channel_route_params(::Val{ROUTE_HYDROLOGIC}; kwargs...) = [
    _required_route_param(kwargs, :storage_coeff, ROUTE_HYDROLOGIC),
    _required_route_param(kwargs, :storage_exponent, ROUTE_HYDROLOGIC),
]

function _channel_route_params(::Val{algorithm}; kwargs...) where {algorithm}
    throw(ArgumentError("Unsupported channel routing algorithm: $(algorithm)"))
end

function build_channel_route(
    algorithm::Symbol;
    input,
    output,
    adjacency=nothing,
    htypes=nothing,
    max_lag::Integer=32,
    name::Union{Nothing,Symbol}=nothing,
    kwargs...,
)
    params = _channel_route_params(Val(algorithm); kwargs...)
    return _build_named_channel_route(
        algorithm,
        params;
        input=input,
        output=output,
        adjacency=adjacency,
        htypes=htypes,
        max_lag=max_lag,
        name=name,
    )
end

no_channel_route(; kwargs...) = build_channel_route(ROUTE_NONE; kwargs...)
plug_flow_route(; kwargs...) = build_channel_route(ROUTE_PLUG_FLOW; kwargs...)
diffusive_wave_route(; kwargs...) = build_channel_route(ROUTE_DIFFUSIVE_WAVE; kwargs...)
linear_storage_route(; kwargs...) = build_channel_route(ROUTE_LINEAR_STORAGE; kwargs...)
storage_coeff_route(; kwargs...) = build_channel_route(ROUTE_STORAGE_COEFF; kwargs...)
muskingum_route(; kwargs...) = build_channel_route(ROUTE_MUSKINGUM; kwargs...)
nash_cascade_route(; kwargs...) = build_channel_route(ROUTE_NASH_CASCADE; kwargs...)
hydrologic_route(; kwargs...) = build_channel_route(ROUTE_HYDROLOGIC; kwargs...)

export ChannelRoute,
       build_channel_route,
       no_channel_route,
       plug_flow_route,
       diffusive_wave_route,
       linear_storage_route,
       storage_coeff_route,
       muskingum_route,
       nash_cascade_route,
       hydrologic_route

end


