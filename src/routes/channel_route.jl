module channel_route

using ..HydroModels
using ..HydroModels: _as_componentvector,
                     tosymbol,
                     RouteIRF,
                     SparseRouteKernel,
                     SparseRouteConvolution,
                     build_irf_kernels,
                     aggregate_route_kernel
using ..unithydro: ROUTE_DIFFUSIVE_WAVE,
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

struct ChannelRoute{AT,OT,HT,NT} <: HydroModels.AbstractComponent
    name::Symbol
    algorithm::Symbol
    adjacency::AT
    topo_order::OT
    max_lag::Int
    htypes::HT
    infos::NT
end

function ChannelRoute(
    inputs::AbstractVector,
    outputs::AbstractVector,
    params::AbstractVector;
    algorithm::Symbol,
    adjacency::Union{Nothing,AbstractMatrix}=nothing,
    max_lag::Integer=32,
    htypes=nothing,
    name::Union{Nothing,Symbol}=nothing,
)
    length(inputs) == length(outputs) == 1 || throw(ArgumentError("ChannelRoute supports one input and one output."))

    htypes = htypes isa Vector{Int} && isempty(htypes) ? nothing : htypes
    adjacency_mat = isnothing(adjacency) ? nothing : Matrix{Float64}(adjacency)
    topo_order = isnothing(adjacency_mat) ? nothing : _topological_order(adjacency_mat)
    param_names = isempty(params) ? Symbol[] : tosymbol.(params)
    infos = HydroModels.HydroInfos(
        inputs=tosymbol.(inputs),
        outputs=tosymbol.(outputs),
        params=param_names,
    )
    route_name = isnothing(name) ? algorithm : name

    return ChannelRoute{typeof(adjacency_mat),typeof(topo_order),typeof(htypes),typeof(infos)}(
        route_name,
        algorithm,
        adjacency_mat,
        topo_order,
        Int(max_lag),
        htypes,
        infos,
    )
end

function _topological_order(adjacency::AbstractMatrix)
    size(adjacency, 1) == size(adjacency, 2) || throw(ArgumentError("adjacency must be square."))
    nn = size(adjacency, 1)
    indegree = [count(!iszero, @view adjacency[row, :]) for row in 1:nn]
    pending = collect(findall(==(0), indegree))
    order = Int[]

    while !isempty(pending)
        node = popfirst!(pending)
        push!(order, node)
        for downstream in 1:nn
            if !iszero(adjacency[downstream, node])
                indegree[downstream] -= 1
                if indegree[downstream] == 0
                    push!(pending, downstream)
                end
            end
        end
    end

    length(order) == nn || throw(ArgumentError("adjacency must define an acyclic upstream-to-downstream network."))
    return order
end

function _scalar_param(value, name::Symbol)
    if value isa Number
        return Float64(value)
    elseif value isa AbstractVector && length(value) == 1
        return Float64(first(value))
    end
    throw(ArgumentError("Parameter $name must be scalar for single-node routing."))
end

function _node_param_array(value, name::Symbol, nn::Int, htypes)
    if !isnothing(htypes)
        if value isa Number
            return fill(Float64(value), length(htypes))
        end
        return Float64.(collect(value[htypes]))
    elseif value isa Number
        return fill(Float64(value), nn)
    elseif value isa AbstractVector && length(value) == nn
        return Float64.(collect(value))
    elseif value isa AbstractVector && length(value) == 1
        return fill(Float64(first(value)), nn)
    end
    throw(ArgumentError("Parameter $name must be scalar, length 1, or length $nn."))
end

function _extract_scalar_params(route::ChannelRoute, params::ComponentVector)
    param_names = collect(HydroModels.get_param_names(route))
    return [_scalar_param(params[:params][name], name) for name in param_names]
end

function _extract_node_params(route::ChannelRoute, params::ComponentVector, nn::Int)
    param_names = collect(HydroModels.get_param_names(route))
    return [_node_param_array(params[:params][name], name, nn, route.htypes) for name in param_names]
end

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
    D = max(diffusivity, _MIN_POSITIVE)

    for lag in 1:length(weights)
        t = max((lag - 0.5) * delta_t, _MIN_POSITIVE)
        weights[lag] = (1 / (2 * sqrt(pi * D * t))) * exp(-((flow_length - c_ref * t)^2) / (4 * D * t))
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

@inline _is_convolution_route(route::ChannelRoute) = route.algorithm == ROUTE_PLUG_FLOW ||
    route.algorithm == ROUTE_DIFFUSIVE_WAVE ||
    route.algorithm == ROUTE_LINEAR_STORAGE ||
    route.algorithm == ROUTE_NASH_CASCADE

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

function _build_route_irf(route::ChannelRoute)
    if route.algorithm == ROUTE_PLUG_FLOW
        return RouteIRF([:flow_length, :wave_celerity], (p, dt, horizon) -> _plug_kernel(p[1], p[2], dt, min(route.max_lag, horizon)); name=route.name)
    elseif route.algorithm == ROUTE_DIFFUSIVE_WAVE
        return RouteIRF([:flow_length, :wave_celerity, :diffusivity], (p, dt, horizon) -> _diffusive_kernel(p[1], p[2], p[3], dt, min(route.max_lag, horizon)); name=route.name)
    elseif route.algorithm == ROUTE_LINEAR_STORAGE
        return RouteIRF([:storage_time], (p, dt, horizon) -> _linear_storage_kernel(p[1], dt, min(route.max_lag, horizon)); name=route.name)
    elseif route.algorithm == ROUTE_NASH_CASCADE
        return RouteIRF([:storage_time, :n_reservoirs], (p, dt, horizon) -> _nash_cascade_kernel(p[1], p[2], dt, min(route.max_lag, horizon)); name=route.name)
    end
    throw(ArgumentError("No RouteIRF available for $(route.algorithm)."))
end

function _simulate_single_convolution(route::ChannelRoute, inflow::AbstractVector, params::AbstractVector, delta_t::Float64)
    horizon = max(1, min(route.max_lag, length(inflow)))
    irf = _build_route_irf(route)
    weights = irf(params; delta_t=delta_t, horizon=horizon)
    kernel = SparseRouteKernel(reshape(Int[1, 1], 1, 2), reshape(weights, 1, :), 1, 1)
    return SparseRouteConvolution(kernel)(Float64.(inflow))
end

function _simulate_single(route::ChannelRoute, inflow::AbstractVector, params::AbstractVector, delta_t::Float64)
    if route.algorithm == ROUTE_NONE
        return Float64.(inflow)
    elseif _is_convolution_route(route)
        return _simulate_single_convolution(route, inflow, params, delta_t)
    elseif route.algorithm == ROUTE_STORAGE_COEFF
        out = zeros(Float64, length(inflow))
        q_prev = 0.0
        q_in_prev = 0.0
        for tidx in eachindex(inflow)
            c1, c2, c3 = _storage_coeffs(params[1], delta_t)
            q_in_curr = Float64(inflow[tidx])
            out[tidx] = max(c1 * q_in_curr + c2 * q_in_prev + c3 * q_prev, 0.0)
            q_prev = out[tidx]
            q_in_prev = q_in_curr
        end
        return out
    elseif route.algorithm == ROUTE_MUSKINGUM
        out = zeros(Float64, length(inflow))
        q_prev = 0.0
        q_in_prev = 0.0
        for tidx in eachindex(inflow)
            c1, c2, c3 = _muskingum_coeffs(params[1], params[2], delta_t)
            q_in_curr = Float64(inflow[tidx])
            out[tidx] = max(c1 * q_in_curr + c2 * q_in_prev + c3 * q_prev, 0.0)
            q_prev = out[tidx]
            q_in_prev = q_in_curr
        end
        return out
    elseif route.algorithm == ROUTE_HYDROLOGIC
        out = zeros(Float64, length(inflow))
        q_prev = 0.0
        q_in_prev = 0.0
        for tidx in eachindex(inflow)
            q_in_curr = Float64(inflow[tidx])
            out[tidx] = _hydrologic_step(q_prev, q_in_prev, q_in_curr, params[1], params[2], delta_t)
            q_prev = out[tidx]
            q_in_prev = q_in_curr
        end
        return out
    end

    throw(ArgumentError("Unsupported channel routing algorithm: $(route.algorithm)"))
end

function _upstream_sum(adjacency, routed, node)
    isnothing(adjacency) && return 0.0
    acc = 0.0
    @inbounds for upstream in 1:length(routed)
        acc += adjacency[node, upstream] * routed[upstream]
    end
    return acc
end

function _simulate_network_convolution(route::ChannelRoute, input::AbstractMatrix, params::AbstractVector, delta_t::Float64)
    nn, nt = size(input)
    irf = _build_route_irf(route)
    node_kernels = build_irf_kernels(irf, params; delta_t=delta_t, horizon=min(route.max_lag, nt))
    kernel = aggregate_route_kernel(route.adjacency, node_kernels; topo_order=route.topo_order, horizon=nt)
    return SparseRouteConvolution(kernel)(Float64.(input))
end

function _simulate_network(route::ChannelRoute, input::AbstractMatrix, params::AbstractVector, delta_t::Float64)
    nn, nt = size(input)
    order = isnothing(route.topo_order) ? collect(1:nn) : route.topo_order
    routed = zeros(Float64, nn, nt)
    prev_out = zeros(Float64, nn)
    prev_local = zeros(Float64, nn)

    if _is_convolution_route(route)
        return _simulate_network_convolution(route, input, params, delta_t)
    end

    for tidx in 1:nt
        curr_out = zeros(Float64, nn)
        for node in order
            upstream_curr = _upstream_sum(route.adjacency, curr_out, node)
            upstream_prev = _upstream_sum(route.adjacency, prev_out, node)
            local_curr = Float64(input[node, tidx])

            if route.algorithm == ROUTE_NONE
                curr_out[node] = upstream_curr + local_curr
            elseif route.algorithm == ROUTE_STORAGE_COEFF
                c1, c2, c3 = _storage_coeffs(params[1][node], delta_t)
                curr_out[node] = max(c1 * upstream_curr + c2 * upstream_prev + c3 * prev_out[node] + local_curr, 0.0)
            elseif route.algorithm == ROUTE_MUSKINGUM
                c1, c2, c3 = _muskingum_coeffs(params[1][node], params[2][node], delta_t)
                curr_out[node] = max(c1 * upstream_curr + c2 * upstream_prev + c3 * prev_out[node] + local_curr, 0.0)
            elseif route.algorithm == ROUTE_HYDROLOGIC
                total_prev = upstream_prev + prev_local[node]
                total_curr = upstream_curr + local_curr
                curr_out[node] = _hydrologic_step(prev_out[node], total_prev, total_curr, params[1][node], params[2][node], delta_t)
            else
                throw(ArgumentError("Unsupported channel routing algorithm: $(route.algorithm)"))
            end
        end
        routed[:, tidx] = curr_out
        prev_out = curr_out
        prev_local = Float64.(input[:, tidx])
    end

    return routed
end

function (route::ChannelRoute)(input::AbstractArray{T,2}, params::AbstractVector, config=HydroModels.default_config(); kwargs...) where {T}
    size(input, 1) == 1 || throw(ArgumentError("ChannelRoute expects a single input variable."))
    params_cv = _as_componentvector(params)
    delta_t = Float64(get(kwargs, :delta_t, 1.0))
    routed = _simulate_single(route, vec(input[1, :]), _extract_scalar_params(route, params_cv), delta_t)
    return reshape(convert.(T, routed), 1, :)
end

function (route::ChannelRoute)(input::AbstractArray{T,3}, params::AbstractVector, config=HydroModels.default_config(); kwargs...) where {T}
    size(input, 1) == 1 || throw(ArgumentError("ChannelRoute expects a single input variable."))
    params_cv = _as_componentvector(params)
    delta_t = Float64(get(kwargs, :delta_t, 1.0))
    routed = _simulate_network(route, input[1, :, :], _extract_node_params(route, params_cv, size(input, 2)), delta_t)
    return reshape(convert.(T, routed), 1, size(routed, 1), size(routed, 2))
end

function build_channel_route(algorithm::Symbol; input, output, adjacency=nothing, htypes=nothing, max_lag::Integer=32, name::Union{Nothing,Symbol}=nothing, kwargs...)
    if algorithm == ROUTE_NONE
        return ChannelRoute([input], [output], HydroModels.Num[]; algorithm=algorithm, adjacency=adjacency, htypes=htypes, max_lag=1, name=name)
    elseif algorithm == ROUTE_PLUG_FLOW
        flow_length = get(kwargs, :flow_length, nothing)
        wave_celerity = get(kwargs, :wave_celerity, nothing)
        isnothing(flow_length) && throw(ArgumentError("ROUTE_PLUG_FLOW requires `flow_length`."))
        isnothing(wave_celerity) && throw(ArgumentError("ROUTE_PLUG_FLOW requires `wave_celerity`."))
        return ChannelRoute([input], [output], [flow_length, wave_celerity]; algorithm=algorithm, adjacency=adjacency, htypes=htypes, max_lag=max_lag, name=name)
    elseif algorithm == ROUTE_DIFFUSIVE_WAVE
        flow_length = get(kwargs, :flow_length, nothing)
        wave_celerity = get(kwargs, :wave_celerity, nothing)
        diffusivity = get(kwargs, :diffusivity, nothing)
        isnothing(flow_length) && throw(ArgumentError("ROUTE_DIFFUSIVE_WAVE requires `flow_length`."))
        isnothing(wave_celerity) && throw(ArgumentError("ROUTE_DIFFUSIVE_WAVE requires `wave_celerity`."))
        isnothing(diffusivity) && throw(ArgumentError("ROUTE_DIFFUSIVE_WAVE requires `diffusivity`."))
        return ChannelRoute([input], [output], [flow_length, wave_celerity, diffusivity]; algorithm=algorithm, adjacency=adjacency, htypes=htypes, max_lag=max_lag, name=name)
    elseif algorithm == ROUTE_LINEAR_STORAGE
        storage_time = get(kwargs, :storage_time, nothing)
        isnothing(storage_time) && throw(ArgumentError("ROUTE_LINEAR_STORAGE requires `storage_time`."))
        return ChannelRoute([input], [output], [storage_time]; algorithm=algorithm, adjacency=adjacency, htypes=htypes, max_lag=max_lag, name=name)
    elseif algorithm == ROUTE_STORAGE_COEFF
        storage_time = get(kwargs, :storage_time, nothing)
        isnothing(storage_time) && throw(ArgumentError("ROUTE_STORAGE_COEFF requires `storage_time`."))
        return ChannelRoute([input], [output], [storage_time]; algorithm=algorithm, adjacency=adjacency, htypes=htypes, max_lag=1, name=name)
    elseif algorithm == ROUTE_MUSKINGUM
        muskingum_k = get(kwargs, :muskingum_k, nothing)
        muskingum_x = get(kwargs, :muskingum_x, nothing)
        isnothing(muskingum_k) && throw(ArgumentError("ROUTE_MUSKINGUM requires `muskingum_k`."))
        isnothing(muskingum_x) && throw(ArgumentError("ROUTE_MUSKINGUM requires `muskingum_x`."))
        return ChannelRoute([input], [output], [muskingum_k, muskingum_x]; algorithm=algorithm, adjacency=adjacency, htypes=htypes, max_lag=1, name=name)
    elseif algorithm == ROUTE_NASH_CASCADE
        storage_time = get(kwargs, :storage_time, nothing)
        n_reservoirs = get(kwargs, :n_reservoirs, nothing)
        isnothing(storage_time) && throw(ArgumentError("ROUTE_NASH_CASCADE requires `storage_time`."))
        isnothing(n_reservoirs) && throw(ArgumentError("ROUTE_NASH_CASCADE requires `n_reservoirs`."))
        return ChannelRoute([input], [output], [storage_time, n_reservoirs]; algorithm=algorithm, adjacency=adjacency, htypes=htypes, max_lag=max_lag, name=name)
    elseif algorithm == ROUTE_HYDROLOGIC
        storage_coeff = get(kwargs, :storage_coeff, nothing)
        storage_exponent = get(kwargs, :storage_exponent, nothing)
        isnothing(storage_coeff) && throw(ArgumentError("ROUTE_HYDROLOGIC requires `storage_coeff`."))
        isnothing(storage_exponent) && throw(ArgumentError("ROUTE_HYDROLOGIC requires `storage_exponent`."))
        return ChannelRoute([input], [output], [storage_coeff, storage_exponent]; algorithm=algorithm, adjacency=adjacency, htypes=htypes, max_lag=1, name=name)
    end

    throw(ArgumentError("Unsupported channel routing algorithm: $algorithm"))
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
