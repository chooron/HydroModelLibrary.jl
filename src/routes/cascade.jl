module nash_cascade

using .HydroModelCore
using .HydroModelLibrary: tosymbol

struct NashCascade <: HydroModelCore.AbstractComponent

    function NashCascade(
        input::T, output::T, k::T, n::T,
        name::Symbol=:nash_cascade
    ) where T
        infos = HydroModelCore.HydroInfos(
            inputs=tosymbol(input),
            outputs=tosymbol(output),
            params=tosymbol.([k, n])
        )
        return new{typeof(infos)}(name, lenF, infos)
    end

end

function (cascade::NashCascade)(input::AbstractArray{T, 2}, pas::ComponentVector; config::NamedTuple=NamedTuple()) where T
    n = Int(pas.params.n)
    init_states = zeros(n)
    input_len = size(input)[2]
    input_itp = LinearInterpolation(input[1, :], collect(1:input_len))
    config = get(config, :solver, HydroModels.ManualSolver())

    function nash_unithydro!(du, u, p, t)
        k = p
        du[1] = input_itp(t) - u[1] / k
        for i in 2:n
            du[i] = (u[i-1] - u[i]) / k
        end
    end

    prob = ODEProblem(nash_unithydro!, init_states, (1, length(input_vec)), (pas[:params].k,))
    sol = solve(prob, Tsit5())
    sol_vec = sol.u[:, end] .* pas[:params].k
    reshape(sol_vec, 1, input_len)
end

end
