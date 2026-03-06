using LinearAlgebra
using Statistics

function _validate_accuracy_inputs(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real})
    n = length(simulated_array)
    n > 0 || throw(ArgumentError("Input arrays must be non-empty."))
    length(observed_array) == n || throw(
        ArgumentError("simulated_array and observed_array must have the same length."),
    )
    return nothing
end

function _paired_series(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real})
    _validate_accuracy_inputs(simulated_array, observed_array)
    return Float64.(simulated_array), Float64.(observed_array)
end

_safe_div(num::Real, den::Real) = den == 0 ? NaN : float(num) / float(den)

function _safe_acos(num::Real, den::Real)
    den == 0 && return NaN
    return acos(clamp(float(num) / float(den), -1.0, 1.0))
end

function _safe_asin(value::Real)
    return asin(clamp(float(value), -1.0, 1.0))
end

function _rankdata(values::AbstractVector{<:Real})
    n = length(values)
    perm = sortperm(values)
    ranks = zeros(Float64, n)

    i = 1
    while i <= n
        j = i
        while j < n && values[perm[j + 1]] == values[perm[i]]
            j += 1
        end

        avg_rank = (i + j) / 2
        @inbounds for k in i:j
            ranks[perm[k]] = avg_rank
        end

        i = j + 1
    end

    return ranks
end

function _geometric_mean(values::AbstractVector{<:Real})
    arr = Float64.(values)
    any(x -> x <= 0, arr) && return NaN
    return exp(mean(log.(arr)))
end

function acc(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    a = sim .- mean(sim)
    b = obs .- mean(obs)
    c = std(obs) * std(sim) * length(sim)
    return _safe_div(dot(a, b), c)
end

function d(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    a = (obs .- sim) .^ 2
    b = abs.(sim .- mean(obs))
    c = abs.(obs .- mean(obs))
    return 1 - _safe_div(sum(a), sum((b .+ c) .^ 2))
end

function d1(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    obs_mean = mean(obs)
    a = abs.(sim .- obs)
    b = abs.(sim .- obs_mean)
    c = abs.(obs .- obs_mean)
    return 1 - _safe_div(sum(a), sum(b .+ c))
end

function d1_p(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    obs_bar_p = get(kwargs, :obs_bar_p, nothing)
    ref = isnothing(obs_bar_p) ? mean(obs) : float(obs_bar_p)

    a = abs.(obs .- sim)
    b = abs.(sim .- ref) .+ abs.(obs .- ref)
    return 1 - _safe_div(sum(a), sum(b))
end

function dmod(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    j = float(get(kwargs, :j, 1.0))

    a = abs.(sim .- obs) .^ j
    b = abs.(sim .- mean(obs))
    c = abs.(obs .- mean(obs))
    e = (b .+ c) .^ j
    return 1 - _safe_div(sum(a), sum(e))
end

function dr(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    a = sum(abs.(sim .- obs))
    b = 2 * sum(abs.(obs .- mean(obs)))

    if a <= b
        return 1 - _safe_div(a, b)
    else
        return _safe_div(b, a) - 1
    end
end

function drel(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    obs_mean = mean(obs)

    a = ((sim .- obs) ./ obs) .^ 2
    b = abs.(sim .- obs_mean)
    c = abs.(obs .- obs_mean)
    e = ((b .+ c) ./ obs_mean) .^ 2
    return 1 - _safe_div(sum(a), sum(e))
end

function ed(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    return norm(obs .- sim)
end

function g_mean_diff(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    sim_log = log1p.(sim)
    obs_log = log1p.(obs)
    return exp(_geometric_mean(sim_log) - _geometric_mean(obs_log))
end

function h1_mhe(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    h = (sim .- obs) ./ obs
    return mean(h)
end

function h1_mahe(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    h = (sim .- obs) ./ obs
    return mean(abs.(h))
end

function h1_rmshe(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    h = (sim .- obs) ./ obs
    return sqrt(mean(h .^ 2))
end

function h10_mhe(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    h = log1p.(sim) .- log1p.(obs)
    return mean(h)
end

function h10_mahe(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    h = log1p.(sim) .- log1p.(obs)
    return mean(abs.(h))
end

function h10_rmshe(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    h = log1p.(sim) .- log1p.(obs)
    return sqrt(mean(h .^ 2))
end

function h2_mhe(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    h = (sim .- obs) ./ sim
    return mean(h)
end

function h2_mahe(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    h = (sim .- obs) ./ sim
    return mean(abs.(h))
end

function h2_rmshe(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    h = (sim .- obs) ./ sim
    return sqrt(mean(h .^ 2))
end

function h3_mhe(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    h = (sim .- obs) ./ (0.5 .* (sim .+ obs))
    return mean(h)
end

function h3_mahe(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    h = (sim .- obs) ./ (0.5 .* (sim .+ obs))
    return mean(abs.(h))
end

function h3_rmshe(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    h = (sim .- obs) ./ (0.5 .* (sim .+ obs))
    return sqrt(mean(h .^ 2))
end

function h4_mhe(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    h = (sim .- obs) ./ sqrt.(sim .* obs)
    return mean(h)
end

function h4_mahe(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    h = (sim .- obs) ./ sqrt.(sim .* obs)
    return mean(abs.(h))
end

function h4_rmshe(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    h = (sim .- obs) ./ sqrt.(sim .* obs)
    return sqrt(mean(h .^ 2))
end

function h5_mhe(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    top = sim .- obs
    bot = 1.0 ./ (0.5 .* ((1.0 ./ obs) .+ (1.0 ./ sim)))
    h = top ./ bot
    return mean(h)
end

function h5_mahe(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    top = sim .- obs
    bot = 1.0 ./ (0.5 .* ((1.0 ./ obs) .+ (1.0 ./ sim)))
    h = top ./ bot
    return mean(abs.(h))
end

function h5_rmshe(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    top = sim .- obs
    bot = 1.0 ./ (0.5 .* ((1.0 ./ obs) .+ (1.0 ./ sim)))
    h = top ./ bot
    return sqrt(mean(h .^ 2))
end

function h6_mhe(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    k = float(get(kwargs, :k, 1.0))

    ratio = sim ./ obs
    top = ratio .- 1
    bot = (0.5 .* (1 .+ ratio .^ k)) .^ (1 / k)
    h = top ./ bot
    return mean(h)
end

function h6_mahe(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    k = float(get(kwargs, :k, 1.0))

    ratio = sim ./ obs
    top = ratio .- 1
    bot = (0.5 .* (1 .+ ratio .^ k)) .^ (1 / k)
    h = top ./ bot
    return mean(abs.(h))
end

function h6_rmshe(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    k = float(get(kwargs, :k, 1.0))

    ratio = sim ./ obs
    top = ratio .- 1
    bot = (0.5 .* (1 .+ ratio .^ k)) .^ (1 / k)
    h = top ./ bot
    return sqrt(mean(h .^ 2))
end

function h7_mhe(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    ratio = sim ./ obs
    h = (ratio .- 1) ./ minimum(ratio)
    return mean(h)
end

function h7_mahe(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    ratio = sim ./ obs
    h = (ratio .- 1) ./ minimum(ratio)
    return mean(abs.(h))
end

function h7_rmshe(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    ratio = sim ./ obs
    h = (ratio .- 1) ./ minimum(ratio)
    return sqrt(mean(h .^ 2))
end

function h8_mhe(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    ratio = sim ./ obs
    h = (ratio .- 1) ./ maximum(ratio)
    return mean(h)
end

function h8_mahe(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    ratio = sim ./ obs
    h = (ratio .- 1) ./ maximum(ratio)
    return mean(abs.(h))
end

function h8_rmshe(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    ratio = sim ./ obs
    h = (ratio .- 1) ./ maximum(ratio)
    return sqrt(mean(h .^ 2))
end

function irmse(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    length(obs) < 2 && return NaN

    obs_grad = diff(obs)
    obs_grad_std = std(obs_grad)

    rmse_value = sqrt(mean((sim .- obs) .^ 2))
    return _safe_div(rmse_value, obs_grad_std)
end

function kge_2009(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    s = get(kwargs, :s, (1.0, 1.0, 1.0))

    pr = person_r(sim, obs)
    beta = _safe_div(mean(sim), mean(obs))
    alpha = _safe_div(std(sim), std(obs))

    any(isnan, (pr, beta, alpha)) && return NaN

    return 1 - sqrt((s[1] * (pr - 1))^2 + (s[2] * (alpha - 1))^2 + (s[3] * (beta - 1))^2)
end

function kge_2012(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    s = get(kwargs, :s, (1.0, 1.0, 1.0))

    pr = person_r(sim, obs)
    beta = _safe_div(mean(sim), mean(obs))

    sim_cv = _safe_div(std(sim), mean(sim))
    obs_cv = _safe_div(std(obs), mean(obs))
    gam = _safe_div(sim_cv, obs_cv)

    any(isnan, (pr, beta, gam)) && return NaN

    return 1 - sqrt((s[1] * (pr - 1))^2 + (s[2] * (gam - 1))^2 + (s[3] * (beta - 1))^2)
end

function lm_index(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    obs_bar_p = get(kwargs, :obs_bar_p, nothing)
    ref = isnothing(obs_bar_p) ? mean(obs) : float(obs_bar_p)

    a = abs.(sim .- obs)
    b = abs.(obs .- ref)
    return 1 - _safe_div(sum(a), sum(b))
end

function maappe(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    return mean(atan.(abs.((sim .- obs) ./ obs)))
end

function mae(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    return mean(abs.(sim .- obs))
end

function male(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    sim_log = log1p.(sim)
    obs_log = log1p.(obs)
    return mean(abs.(sim_log .- obs_log))
end

function mapd(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    return _safe_div(sum(abs.(sim .- obs)), sum(abs.(obs)))
end

function mape(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    return 100 / length(sim) * sum(abs.((sim .- obs) ./ obs))
end

function mase(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    m = Int(get(kwargs, :m, 1))

    n = length(obs)
    (m < 1 || n <= m) && return NaN

    mae_model = mean(abs.(sim .- obs))
    mae_naive = mean(abs.(obs[(m + 1):end] .- obs[1:(end - m)]))
    return _safe_div(mae_model, mae_naive)
end

function mdae(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    return median(abs.(sim .- obs))
end

function mde(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    return median(sim .- obs)
end

function mdse(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    return median((sim .- obs) .^ 2)
end

function me(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    return mean(sim .- obs)
end

function mean_var(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    return var(log1p.(obs) .- log1p.(sim))
end

function mhe_h1 end
function mhe_h2 end
function mhe_h3 end
function mhe_h4 end

function mle(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    sim_log = log1p.(sim)
    obs_log = log1p.(obs)
    return mean(sim_log .- obs_log)
end

function mse(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    return mean((sim .- obs) .^ 2)
end

function msle(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    sim_log = log1p.(sim)
    obs_log = log1p.(obs)
    return mean((sim_log .- obs_log) .^ 2)
end

function ned(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    a = obs ./ mean(obs)
    b = sim ./ mean(sim)
    return norm(a .- b)
end

function nrmse_range(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    rmse_value = rmse(sim, obs)
    obs_range = maximum(obs) - minimum(obs)
    return _safe_div(rmse_value, obs_range)
end

function nrmse_mean(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    rmse_value = rmse(sim, obs)
    return _safe_div(rmse_value, mean(obs))
end

function nrmse_iqr(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    rmse_value = rmse(sim, obs)
    q = quantile(obs, [0.25, 0.75])
    iqr = q[2] - q[1]
    return _safe_div(rmse_value, iqr)
end

function nse(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    a = (sim .- obs) .^ 2
    b = (obs .- mean(obs)) .^ 2
    return 1 - _safe_div(sum(a), sum(b))
end

function nse_mod(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    j = float(get(kwargs, :j, 1.0))

    a = abs.(sim .- obs) .^ j
    b = abs.(obs .- mean(obs)) .^ j
    return 1 - _safe_div(sum(a), sum(b))
end

function nse_rel(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    obs_mean = mean(obs)

    a = abs.((sim .- obs) ./ obs) .^ 2
    b = abs.((obs .- obs_mean) ./ obs_mean) .^ 2
    return 1 - _safe_div(sum(a), sum(b))
end

function person_r(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    sim_centered = sim .- mean(sim)
    obs_centered = obs .- mean(obs)

    top = sum(obs_centered .* sim_centered)
    bot = sqrt(sum(obs_centered .^ 2) * sum(sim_centered .^ 2))
    return _safe_div(top, bot)
end

function spearman_r(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    rank_sim = _rankdata(sim)
    rank_obs = _rankdata(obs)
    return person_r(rank_sim, rank_obs)
end

function r_squared(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    r = person_r(simulated_array, observed_array)
    return r^2
end

function mb_r(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    n = length(sim)

    tot = 0.0
    for i in 1:n
        tot += sum(abs.(sim .- obs[i]))
    end

    mae_val = mean(abs.(sim .- obs))
    return 1 - _safe_div((n^2) * mae_val, tot)
end

function r2(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    return r_squared(simulated_array, observed_array)
end

function rmse(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    return sqrt(mean((sim .- obs) .^ 2))
end

function rmsle(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    return sqrt(mean((log1p.(sim) .- log1p.(obs)) .^ 2))
end

function sa(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    return _safe_acos(dot(sim, obs), norm(sim) * norm(obs))
end

function sc(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    sim_centered = sim .- mean(sim)
    obs_centered = obs .- mean(obs)

    return _safe_acos(dot(obs_centered, sim_centered), norm(obs_centered) * norm(sim_centered))
end

function sga(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    length(sim) < 2 && return NaN

    sgx = diff(obs)
    sgy = diff(sim)
    return _safe_acos(dot(sgx, sgy), norm(sgx) * norm(sgy))
end

function sid(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)

    first = (obs ./ mean(obs)) .- (sim ./ mean(sim))
    second = (log10.(obs) .- log10(mean(obs))) .- (log10.(sim) .- log10(mean(sim)))
    return dot(first, second)
end

function smape1(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    return 100 / length(sim) * sum(abs.(sim .- obs) ./ (abs.(sim) .+ abs.(obs)))
end

function smape2(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    return 100 / length(sim) * sum(abs.((sim .- obs) ./ ((sim .+ obs) ./ 2)))
end

function ve(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)
    return 1 - _safe_div(sum(abs.(sim .- obs)), sum(obs))
end

function watt_m(simulated_array::AbstractVector{<:Real}, observed_array::AbstractVector{<:Real}; kwargs...)
    sim, obs = _paired_series(simulated_array, observed_array)

    a = 2 / pi
    b = mean((sim .- obs) .^ 2)
    c = std(obs)^2 + std(sim)^2
    e = (mean(sim) - mean(obs))^2
    f = c + e

    f == 0 && return NaN
    return a * _safe_asin(1 - b / f)
end