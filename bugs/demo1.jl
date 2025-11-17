using DifferentiationInterface, Enzyme
using DifferentialEquations, DataInterpolations, SciMLSensitivity
using ComponentArrays
using Zygote

t_data = collect(0.0:1.0:11.0)
alpha_values = rand(length(t_data))

function loss_function(params, config)

    interpolator = config.interpolator_type(
        alpha_values,
        config.t_data,
    )

    function my_ode!(du, u, p, t)
        alpha_t = interpolator(t)
        du[1] = -alpha_t * u[1] * p[1] + p[2]
        return nothing
    end

    u0 = [10.0]
    tspan = (1.0, 10.0)

    prob = ODEProblem(my_ode!, u0, tspan, params, sensealg=InterpolatingAdjoint(autojacvec=ZygoteVJP()))
    sol = solve(prob, Tsit5(), save_everystep=false, save_start=false)
    return Array(sol) |> sum
end

config = (
    interpolator_type=ConstantInterpolation,
    t_data=t_data,
)

initial_alpha_values = rand(length(t_data))
params = [1.0, 2.0]
target_func(p) = loss_function(p, config)
target_func(params)
ad_backend = AutoZygote()
val, grad = value_and_gradient(target_func, ad_backend, params)