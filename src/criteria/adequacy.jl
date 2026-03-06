using Dates
using Statistics

const ADEQUACY_PE_LIMIT = 50.0
const WATER_YEAR_START_MONTH = 10
const WATER_YEAR_START_DAY = 1

"""
    percentage_error(sim, obs)

Compute percentage error:
`PE = 100 * (sim - obs) / obs`.
"""
percentage_error(sim::Real, obs::Real) = 100.0 * (float(sim) - float(obs)) / float(obs)

function percentage_error(sim::AbstractVector{<:Real}, obs::AbstractVector{<:Real})
    length(sim) == length(obs) || throw(ArgumentError("sim and obs must have the same length."))
    out = Vector{Float64}(undef, length(sim))
    @inbounds for i in eachindex(sim, obs)
        out[i] = percentage_error(sim[i], obs[i])
    end
    return out
end

"""
    is_adequate(pe; limit=ADEQUACY_PE_LIMIT)

Return `true` when `pe` is finite and inside `+/-limit`.
"""
is_adequate(pe::Real; limit::Real=ADEQUACY_PE_LIMIT) = isfinite(pe) && abs(float(pe)) <= float(limit)

function is_adequate(pe::AbstractVector{<:Real}; limit::Real=ADEQUACY_PE_LIMIT)
    out = Vector{Bool}(undef, length(pe))
    @inbounds for i in eachindex(pe)
        out[i] = is_adequate(pe[i]; limit=limit)
    end
    return out
end

"""
    monthly_mean_daily_flow(Q, dates)

Return a 12-element vector with mean daily flow for each calendar month
aggregated over all years in `dates`. Missing months are `NaN`.
"""
function monthly_mean_daily_flow(Q::AbstractVector{<:Real}, dates::AbstractVector{<:Date})
    sums = zeros(Float64, 12)
    counts = zeros(Int, 12)

    @inbounds for i in eachindex(Q, dates)
        m = month(dates[i])
        sums[m] += float(Q[i])
        counts[m] += 1
    end

    means = fill(NaN, 12)
    @inbounds for m in 1:12
        if counts[m] > 0
            means[m] = sums[m] / counts[m]
        end
    end
    return means
end

"""
    yearly_mean_daily_flow(Q, dates)

Return `(years, values)` where `values` are mean daily flows per calendar year.
"""
function yearly_mean_daily_flow(Q::AbstractVector{<:Real}, dates::AbstractVector{<:Date})
    first_year = year(first(dates))
    last_year = year(last(dates))
    years = collect(first_year:last_year)

    sums = zeros(Float64, length(years))
    counts = zeros(Int, length(years))
    @inbounds for i in eachindex(Q, dates)
        idx = year(dates[i]) - first_year + 1
        sums[idx] += float(Q[i])
        counts[idx] += 1
    end

    values = fill(NaN, length(years))
    @inbounds for i in eachindex(values)
        if counts[i] > 0
            values[i] = sums[i] / counts[i]
        end
    end

    return (years=years, values=values)
end

"""
    runoff_ratio(Q, P)

Runoff ratio `RR = sum(Q) / sum(P)`.
"""
function runoff_ratio(Q::AbstractVector{<:Real}, P::AbstractVector{<:Real})
    return sum(float, Q) / sum(float, P)
end

function _event_lengths(mask::AbstractVector{Bool})
    lengths = Int[]
    current = 0

    @inbounds for flag in mask
        if flag
            current += 1
        elseif current > 0
            push!(lengths, current)
            current = 0
        end
    end

    if current > 0
        push!(lengths, current)
    end

    return lengths
end

function _events_per_year(mask::AbstractVector{Bool}, dates::AbstractVector{<:Date})
    lengths = _event_lengths(mask)
    event_count = length(lengths)
    year_span = (Dates.value(last(dates) - first(dates)) + 1) / 365.25
    year_span = max(year_span, eps(Float64))
    return event_count / year_span
end

function _mean_event_duration(mask::AbstractVector{Bool})
    lengths = _event_lengths(mask)
    return isempty(lengths) ? 0.0 : mean(lengths)
end

"""
    bfi_ioh1980(Q)

Baseflow index (BFI) from the IoH (1980) continuous baseflow separation
based on non-overlapping 5-day minima and turning points.
"""
function bfi_ioh1980(Q::AbstractVector{<:Real})
    n = length(Q)
    n == 0 && return NaN

    flow = Float64.(Q)
    total_flow = sum(flow)
    total_flow == 0.0 && return NaN

    block_size = 5
    nblocks = cld(n, block_size)
    block_min = Vector{Float64}(undef, nblocks)
    block_idx = Vector{Int}(undef, nblocks)

    @inbounds for b in 1:nblocks
        s = (b - 1) * block_size + 1
        e = min(b * block_size, n)

        min_val = Inf
        min_pos = s
        for i in s:e
            qi = flow[i]
            if qi < min_val
                min_val = qi
                min_pos = i
            end
        end

        block_min[b] = min_val
        block_idx[b] = min_pos
    end

    turning_blocks = Int[1]
    @inbounds for b in 2:(nblocks - 1)
        c = block_min[b]
        if c <= 0.9 * block_min[b - 1] && c <= 0.9 * block_min[b + 1]
            push!(turning_blocks, b)
        end
    end
    if nblocks > 1
        push!(turning_blocks, nblocks)
    end

    sort!(unique!(turning_blocks))
    tp_idx = block_idx[turning_blocks]
    tp_val = block_min[turning_blocks]

    baseflow = zeros(Float64, n)

    if length(tp_idx) == 1
        baseflow .= tp_val[1]
    else
        @inbounds for seg in 1:(length(tp_idx) - 1)
            i1 = tp_idx[seg]
            i2 = tp_idx[seg + 1]
            q1 = tp_val[seg]
            q2 = tp_val[seg + 1]
            span = i2 - i1

            if span == 0
                baseflow[i1] = q1
                continue
            end

            for i in i1:i2
                w = (i - i1) / span
                baseflow[i] = q1 + w * (q2 - q1)
            end
        end
    end

    first_tp = tp_idx[1]
    @inbounds for i in 1:first_tp
        baseflow[i] = tp_val[1]
    end

    last_tp = tp_idx[end]
    @inbounds for i in last_tp:n
        baseflow[i] = tp_val[end]
    end

    @inbounds for i in eachindex(baseflow, flow)
        baseflow[i] = clamp(baseflow[i], 0.0, flow[i])
    end

    return sum(baseflow) / total_flow
end

function _water_year_label(d::Date)
    return month(d) >= WATER_YEAR_START_MONTH ? year(d) + 1 : year(d)
end

"""
    mhfd(Q, dates)

Mean half-flow date in day-of-water-year (water year starts on Oct 1).
"""
function mhfd(Q::AbstractVector{<:Real}, dates::AbstractVector{<:Date})
    n = length(Q)
    n == 0 && return NaN

    all_days = Float64[]
    complete_days = Float64[]

    i = 1
    while i <= n
        wy = _water_year_label(dates[i])
        j = i
        while j < n && _water_year_label(dates[j + 1]) == wy
            j += 1
        end

        total_q = sum(float, @view Q[i:j])
        if total_q > 0.0
            half_q = 0.5 * total_q
            cumulative = 0.0
            half_idx = i

            @inbounds for k in i:j
                cumulative += float(Q[k])
                if cumulative >= half_q
                    half_idx = k
                    break
                end
            end

            wy_start = Date(wy - 1, WATER_YEAR_START_MONTH, WATER_YEAR_START_DAY)
            day_of_wy = Dates.value(dates[half_idx] - wy_start) + 1
            push!(all_days, day_of_wy)

            is_complete = dates[i] == wy_start && dates[j] == Date(wy, 9, 30)
            if is_complete
                push!(complete_days, day_of_wy)
            end
        end

        i = j + 1
    end

    selected = isempty(complete_days) ? all_days : complete_days
    return isempty(selected) ? NaN : mean(selected)
end

"""
    slope_fdc(Q)

Slope of the flow duration curve between the 33rd and 66th percentiles
in log space.
"""
function slope_fdc(Q::AbstractVector{<:Real})
    flow = Float64.(Q)
    q33 = max(quantile(flow, 0.33), eps(Float64))
    q66 = max(quantile(flow, 0.66), eps(Float64))
    return (log(q66) - log(q33)) / (0.66 - 0.33)
end

function _compute_signatures(Q::AbstractVector{<:Real}, P::AbstractVector{<:Real}, dates::AbstractVector{<:Date})
    qmonth = monthly_mean_daily_flow(Q, dates)
    ymean = yearly_mean_daily_flow(Q, dates)

    qmean = mean(float, Q)
    qmedian = median(Float64.(Q))

    low_thr = 0.2 * qmean
    high_thr = 9.0 * qmedian

    low_mask = Float64.(Q) .< low_thr
    high_mask = Float64.(Q) .> high_thr

    return (
        Qmonth=qmonth,
        Qyear=ymean.values,
        RR=runoff_ratio(Q, P),
        BFI=bfi_ioh1980(Q),
        Q5=quantile(Float64.(Q), 0.05),
        Q95=quantile(Float64.(Q), 0.95),
        LFfreq=_events_per_year(low_mask, dates),
        HFfreq=_events_per_year(high_mask, dates),
        LFdur=_mean_event_duration(low_mask),
        HFdur=_mean_event_duration(high_mask),
        MHFD=mhfd(Q, dates),
        slopeFDC=slope_fdc(Q),
    )
end

function _map_namedtuple_binary(f, a::NamedTuple, b::NamedTuple)
    keys = propertynames(a)
    vals = ntuple(i -> f(getfield(a, i), getfield(b, i)), length(keys))
    return NamedTuple{keys}(vals)
end

function _map_namedtuple_unary(f, a::NamedTuple)
    keys = propertynames(a)
    vals = ntuple(i -> f(getfield(a, i)), length(keys))
    return NamedTuple{keys}(vals)
end

function _validate_daily_inputs(
    Q_sim::AbstractVector{<:Real},
    Q_obs::AbstractVector{<:Real},
    P::AbstractVector{<:Real},
    dates::AbstractVector{<:Date},
)
    n = length(Q_sim)
    n > 0 || throw(ArgumentError("Input arrays must be non-empty."))
    length(Q_obs) == n || throw(ArgumentError("Q_obs must have the same length as Q_sim."))
    length(P) == n || throw(ArgumentError("P must have the same length as Q_sim."))
    length(dates) == n || throw(ArgumentError("dates must have the same length as Q_sim."))

    @inbounds for i in 2:n
        dates[i] == dates[i - 1] + Day(1) || throw(
            ArgumentError("Input dates must be strictly daily with 1-day increments."),
        )
    end

    return nothing
end

"""
    evaluate_adequacy_assessment(Q_sim, Q_obs, P, dates; pe_limit=ADEQUACY_PE_LIMIT)

Evaluate adequacy for daily hydrological series using:

Part A:
- Qmonth: mean daily flow by calendar month over the tested period.
- Qyear: mean daily flow by calendar year.

Part B (10 signatures):
- RR, BFI (IoH 1980), Q5, Q95, LFfreq, HFfreq, LFdur, HFdur, MHFD, slopeFDC.

Core metric:
- PE = 100 * (sim - obs) / obs

Adequacy threshold:
- true if PE is within +/-50%.

Returns a `NamedTuple` with simulated values, observed values,
percentage errors, and adequacy flags.
"""
function evaluate_adequacy_assessment(
    Q_sim::AbstractVector{<:Real},
    Q_obs::AbstractVector{<:Real},
    P::AbstractVector{<:Real},
    dates::AbstractVector{<:Date};
    pe_limit::Real=ADEQUACY_PE_LIMIT,
)
    _validate_daily_inputs(Q_sim, Q_obs, P, dates)

    sim = _compute_signatures(Q_sim, P, dates)
    obs = _compute_signatures(Q_obs, P, dates)

    pe = _map_namedtuple_binary(percentage_error, sim, obs)
    adequacy = _map_namedtuple_unary(x -> is_adequate(x; limit=pe_limit), pe)

    qyear_meta = yearly_mean_daily_flow(Q_obs, dates)

    return (
        metadata=(
            frequency=:daily,
            months=collect(1:12),
            years=qyear_meta.years,
        ),
        simulated=sim,
        observed=obs,
        percentage_error=pe,
        adequate=adequacy,
    )
end
