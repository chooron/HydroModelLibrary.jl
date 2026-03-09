module unithydro

using ..HydroModels
using ..HydroModels: tosymbol
using SpecialFunctions

const ROUTE_DUMP = :ROUTE_DUMP
const ROUTE_GAMMA_CONVOLUTION = :ROUTE_GAMMA_CONVOLUTION
const ROUTE_TRI_CONVOLUTION = :ROUTE_TRI_CONVOLUTION
const ROUTE_RESERVOIR_SERIES = :ROUTE_RESERVOIR_SERIES
const ROUTE_DIFFUSIVE_WAVE = :ROUTE_DIFFUSIVE_WAVE
const ROUTE_NONE = :ROUTE_NONE
const ROUTE_PLUG_FLOW = :ROUTE_PLUG_FLOW
const ROUTE_LINEAR_STORAGE = :ROUTE_LINEAR_STORAGE
const ROUTE_STORAGE_COEFF = :ROUTE_STORAGE_COEFF
const ROUTE_MUSKINGUM = :ROUTE_MUSKINGUM
const ROUTE_NASH_CASCADE = :ROUTE_NASH_CASCADE
const ROUTE_HYDROLOGIC = :ROUTE_HYDROLOGIC

const TOC_MCDERMOTT_PILGRIM = :TOC_MCDERMOTT_PILGRIM
const TOC_BRANSBY_WILLIAMS = :TOC_BRANSBY_WILLIAMS
const TOC_WILLIAMS_1922 = :TOC_WILLIAMS_1922

@inline _normalize_htypes(htypes) = htypes isa Vector{Int} && isempty(htypes) ? nothing : htypes

@inline function _as_inputs_outputs(input, output)
    return [input], [output]
end

function _build_weighted_unit_hydrograph(
    input,
    output,
    params,
    weight_at::Function,
    max_lag_func::Function;
    name::Union{Nothing,Symbol}=nothing,
    htypes=nothing,
)
    inputs, outputs = _as_inputs_outputs(input, output)
    htypes = _normalize_htypes(htypes)

    function cumulative_uh(step, pas)
        max_lag = max(1, Int(ceil(max_lag_func(pas))))
        upper = min(max(Int(floor(step)), 0), max_lag)
        value = 0.0
        @inbounds for idx in 1:upper
            value += weight_at(idx, pas)
        end
        return value
    end

    return HydroModels.UnitHydrograph(
        inputs,
        outputs,
        params,
        cumulative_uh,
        pas -> max(1, Int(ceil(max_lag_func(pas))));
        name=isnothing(name) ? Symbol("##uh#", hash((inputs, outputs, tosymbol.(params)))) : name,
        htypes=htypes,
    )
end

function dump_unit_hydrograph(; input, output, name::Symbol=ROUTE_DUMP, htypes=nothing)
    return _build_weighted_unit_hydrograph(
        input,
        output,
        HydroModels.Num[],
        (step, _) -> step == 1 ? 1.0 : 0.0,
        _ -> 1;
        name=name,
        htypes=htypes,
    )
end

function gamma_unit_hydrograph(;
    input,
    output,
    alpha,
    beta=nothing,
    tp=nothing,
    max_lag::Integer=32,
    name::Symbol=ROUTE_GAMMA_CONVOLUTION,
    htypes=nothing,
)
    alpha_name = tosymbol(alpha)
    beta_name = isnothing(beta) ? nothing : tosymbol(beta)
    tp_name = isnothing(tp) ? nothing : tosymbol(tp)

    if isnothing(beta_name) && isnothing(tp_name)
        throw(ArgumentError("gamma_unit_hydrograph requires either `beta` or `tp`."))
    end

    params = isnothing(beta_name) ? [alpha, tp] : [alpha, beta]

    function weight_at(step, pas)
        a = max(Float64(pas.params[alpha_name]), 1.0e-6)
        b = if isnothing(beta_name)
            max((a - 1.0) / max(Float64(pas.params[tp_name]), 1.0e-6), 1.0e-6)
        else
            max(Float64(pas.params[beta_name]), 1.0e-6)
        end
        t_mid = max(step - 0.5, 1.0e-6)
        return (1 / t_mid) * (b * t_mid)^a / gamma(a) * exp(-b * t_mid)
    end

    return _build_weighted_unit_hydrograph(
        input,
        output,
        params,
        weight_at,
        _ -> max_lag;
        name=name,
        htypes=htypes,
    )
end

function triangular_unit_hydrograph(;
    input,
    output,
    tp,
    te,
    name::Symbol=ROUTE_TRI_CONVOLUTION,
    htypes=nothing,
)
    tp_name = tosymbol(tp)
    te_name = tosymbol(te)

    function weight_at(step, pas)
        tp_val = max(Float64(pas.params[tp_name]), 1.0e-6)
        te_val = max(Float64(pas.params[te_name]), tp_val + 1.0e-6)
        t_mid = step - 0.5
        if t_mid < tp_val
            return (2.0 / te_val) * (t_mid / tp_val)
        elseif t_mid <= te_val
            return (2.0 / te_val) * ((te_val - t_mid) / (te_val - tp_val))
        end
        return 0.0
    end

    return _build_weighted_unit_hydrograph(
        input,
        output,
        [tp, te],
        weight_at,
        pas -> max(1, Int(ceil(Float64(pas.params[te_name]))));
        name=name,
        htypes=htypes,
    )
end

function nash_unit_hydrograph(;
    input,
    output,
    n,
    k,
    max_lag::Integer=32,
    name::Symbol=ROUTE_RESERVOIR_SERIES,
    htypes=nothing,
)
    n_name = tosymbol(n)
    k_name = tosymbol(k)

    function weight_at(step, pas)
        n_val = max(Float64(pas.params[n_name]), 1.0e-6)
        k_val = max(Float64(pas.params[k_name]), 1.0e-6)
        t_mid = max(step - 0.5, 1.0e-6)
        return t_mid^(n_val - 1.0) * k_val^n_val * exp(-k_val * t_mid) / gamma(n_val)
    end

    return _build_weighted_unit_hydrograph(
        input,
        output,
        [n, k],
        weight_at,
        _ -> max_lag;
        name=name,
        htypes=htypes,
    )
end

function build_unit_hydrograph(algorithm::Symbol; kwargs...)
    if algorithm == ROUTE_DUMP
        return dump_unit_hydrograph(; kwargs...)
    elseif algorithm == ROUTE_GAMMA_CONVOLUTION
        return gamma_unit_hydrograph(; kwargs...)
    elseif algorithm == ROUTE_TRI_CONVOLUTION
        return triangular_unit_hydrograph(; kwargs...)
    elseif algorithm == ROUTE_RESERVOIR_SERIES
        return nash_unit_hydrograph(; kwargs...)
    end
    throw(ArgumentError("Unsupported in-catchment routing algorithm: $algorithm"))
end

@inline toc_mcdermott_pilgrim(A) = 0.031667 * A^0.38
@inline toc_bransby_williams(A; L, S) = 0.0359 * L * S^(-0.2) * A^(-0.1)
@inline toc_williams_1922(A; L, S) = 0.02539 * L * S^(-0.2) * A^(-0.1)

function time_of_concentration(algorithm::Symbol, A; L=nothing, S=nothing)
    if algorithm == TOC_MCDERMOTT_PILGRIM
        return toc_mcdermott_pilgrim(A)
    elseif algorithm == TOC_BRANSBY_WILLIAMS
        isnothing(L) && throw(ArgumentError("TOC_BRANSBY_WILLIAMS requires `L`."))
        isnothing(S) && throw(ArgumentError("TOC_BRANSBY_WILLIAMS requires `S`."))
        return toc_bransby_williams(A; L=L, S=S)
    elseif algorithm == TOC_WILLIAMS_1922
        isnothing(L) && throw(ArgumentError("TOC_WILLIAMS_1922 requires `L`."))
        isnothing(S) && throw(ArgumentError("TOC_WILLIAMS_1922 requires `S`."))
        return toc_williams_1922(A; L=L, S=S)
    end
    throw(ArgumentError("Unsupported time-of-concentration algorithm: $algorithm"))
end

export ROUTE_DUMP,
       ROUTE_GAMMA_CONVOLUTION,
       ROUTE_TRI_CONVOLUTION,
       ROUTE_RESERVOIR_SERIES,
       ROUTE_DIFFUSIVE_WAVE,
       ROUTE_NONE,
       ROUTE_PLUG_FLOW,
       ROUTE_LINEAR_STORAGE,
       ROUTE_STORAGE_COEFF,
       ROUTE_MUSKINGUM,
       ROUTE_NASH_CASCADE,
       ROUTE_HYDROLOGIC,
       TOC_MCDERMOTT_PILGRIM,
       TOC_BRANSBY_WILLIAMS,
       TOC_WILLIAMS_1922,
       build_unit_hydrograph,
       dump_unit_hydrograph,
       gamma_unit_hydrograph,
       triangular_unit_hydrograph,
       nash_unit_hydrograph,
       toc_mcdermott_pilgrim,
       toc_bransby_williams,
       toc_williams_1922,
       time_of_concentration

end





