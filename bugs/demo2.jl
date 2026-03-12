using Pkg
Pkg.activate(".")

using DifferentiationInterface
using Enzyme
using DifferentialEquations, SciMLSensitivity
using ComponentArrays

# --- 1. 数据准备 ---
# 既然是均匀时间步，我们只需要知道起始时间和步长
const t_start = 0.0
const t_step = 1.0
# 注意：alpha_values 放在这里，稍后传入闭包
alpha_values_data = rand(12) # 对应 0.0 到 11.0

# --- 2. 核心：手动实现高性能插值 (替代 Interpolations.jl) ---
# 这是一个纯数学函数，Enzyme 可以完美求导，绝不崩溃
function manual_constant_interpolation(t, values, t0, dt)
    # 计算当前时间对应第几个格点
    # 逻辑：(当前时间 - 起始时间) / 步长，向下取整，然后 +1 得到 Julia 的 1-based 索引
    # 加上 1e-9 是为了防止浮点数精度问题导致 t=1.0 变成 0.99999
    idx = floor(Int, (t - t0 + 1e-9) / dt) + 1
    
    # 边界保护：防止求解器尝试计算范围外的时间（例如 t=11.00001）
    idx = clamp(idx, 1, length(values))
    
    # 直接返回数组中的值
    return values[idx]
end

# --- 3. 损失函数 ---
function loss_function(params, alpha_values)
    
    # 定义 ODE
    function my_ode!(du, u, p, t)
        # 解包参数
        param_1, param_2 = p
        
        # 【关键修改】调用手动插值函数，而不是库函数
        # 这行代码等价于 ConstantInterpolation，但对 AD 友好
        alpha_t = manual_constant_interpolation(t, alpha_values, t_start, t_step)
        
        du[1] = -alpha_t * u[1] * param_1 + param_2
        return nothing
    end

    u0 = [10.0]
    tspan = (1.0, 10.0)

    # 定义问题
    # 注意：EnzymeVJP 处理闭包（捕获 alpha_values）比处理复杂的 Struct 稳定得多
    prob = ODEProblem(
        my_ode!, 
        u0, 
        tspan, 
        params, 
        sensealg=GaussAdjoint(autojacvec=EnzymeVJP()) # 你要求的最高效后端
    )
    
    sol = solve(prob, Tsit5(), save_everystep=false, save_start=false)
    
    # 简单的防崩溃检查
    if sol.retcode != ReturnCode.Success
        return 1e9
    end
    
    return sum(sol.u[end])
end

# --- 4. 运行测试 ---
params = [1.0, 2.0]

# 封装目标函数
target_func(p) = loss_function(p, alpha_values_data)

println("Forward pass value: ", target_func(params))

# 设置 Enzyme
ad_backend = AutoEnzyme(mode=Enzyme.Reverse)

println("Calculating Gradient with Enzyme...")

# 这次应该可以顺利跑通，不会有 UndefVarError，也不会有 Segfault
val, grad = value_and_gradient(target_func, ad_backend, params)

println("Value: ", val)
println("Gradient: ", grad)