# HydroModelLibrary.jl

HydroModelLibrary.jl 是一个基于 [HydroModels.jl](https://github.com/chooron/HydroModels.jl) 构建的概念性水文模型库（v0.2.0），提供 40+ 个经典水文模型、多种路由算法、通量过程模块以及模型评估指标。

---

## 安装

```julia
# 从本地路径开发安装
using Pkg
Pkg.develop(path="path/to/HydroModelLibrary")

# 或直接 add（发布后）
Pkg.add("HydroModelLibrary")
```

---

## 快速开始

```julia
using HydroModelLibrary
using ComponentArrays
using HydroModels

# 1. 加载模型
m = load_model(:gr4j)
model = m.model

# 2. 加载示例数据
data = load_sample_data(:gr4j)

# 3. 查看参数范围与随机参数
bounds = get_params_bounds(:gr4j)
params = get_random_params(:gr4j)

# 4. 准备输入（以 gr4j 为例）
input = (P=data[!, "prec"], Ep=data[!, "pet"])
input_arr = stack(input, dims=1)

# 5. 初始化状态
state_names = HydroModels.get_state_names(model)
init_states = ComponentVector(NamedTuple{Tuple(state_names)}(zeros(length(state_names))))

# 6. 运行模型
result = model(input_arr, params, initstates=init_states)
```

---

## 核心 API

### 模型管理

| 函数 | 说明 |
|------|------|
| `load_model(name::Symbol)` | 加载指定模型，返回模型模块（含 `model`、`model_parameters` 等字段） |
| `get_params_bounds(name::Symbol)` | 获取模型参数的上下界，返回 `NamedTuple` |
| `get_random_params(name::Symbol)` | 在参数范围内随机采样，返回 `ComponentVector` |
| `AVAILABLE_MODELS` | 所有可用模型名称的列表（`Vector{Symbol}`） |

```julia
# 查看所有可用模型
println(AVAILABLE_MODELS)

# 加载模型并查看参数
m = load_model(:hymod)
bounds = get_params_bounds(:hymod)   # => (Smax=(1.0, 900.0), b=(...), ...)
params = get_random_params(:hymod)   # => ComponentVector(params=(...))
```

### 示例数据加载

| 函数 | 说明 |
|------|------|
| `load_sample_data(name::Symbol)` | 加载指定数据集，返回 `DataFrame` |
| `load_sample_data_for_model(name::Symbol)` | 按模型名自动匹配数据集 |
| `list_sample_data()` | 打印所有可用数据集及描述 |
| `get_sample_data_info(name::Symbol)` | 获取数据集元信息（不加载数据） |
| `AVAILABLE_SAMPLE_DATA` | 所有可用数据集名称列表 |

```julia
data = load_sample_data(:exphydro)   # => DataFrame
info = get_sample_data_info(:gr4j)   # => (name, description, path, full_path)
list_sample_data()                   # 打印所有数据集信息
```

### 路由算法

**单元线路由（集总式）**

```julia
@variables q_in q_out
@parameters alpha beta

uh = build_unit_hydrograph(ROUTE_GAMMA_CONVOLUTION;
    input=q_in, output=q_out, alpha=alpha, beta=beta, max_lag=24)
result = uh(input_matrix, ComponentVector(params=(alpha=3.0, beta=1.5)))
```

可用路由类型常量：`ROUTE_DUMP`, `ROUTE_GAMMA_CONVOLUTION`, `ROUTE_TRI_CONVOLUTION`,
`ROUTE_RESERVOIR_SERIES`, `ROUTE_DIFFUSIVE_WAVE`, `ROUTE_NONE`, `ROUTE_PLUG_FLOW`,
`ROUTE_LINEAR_STORAGE`, `ROUTE_STORAGE_COEFF`, `ROUTE_MUSKINGUM`, `ROUTE_NASH_CASCADE`, `ROUTE_HYDROLOGIC`

**河道路由（分布式）**

```julia
route = build_channel_route(ROUTE_MUSKINGUM;
    input=q_in, output=q_out,
    muskingum_k=muskingum_k, muskingum_x=muskingum_x)
result = route(input_matrix, ComponentVector(params=(muskingum_k=1.0, muskingum_x=0.0)))
```

**集水时间**

```julia
tc = time_of_concentration(TOC_MCDERMOTT_PILGRIM, 100.0)
```

### 模型评估

**精度指标（accuracy）**

```julia
# Nash-Sutcliffe 效率系数
nse(q_obs, q_sim)

# Kling-Gupta 效率系数
kge_2009(q_obs, q_sim)
kge_2012(q_obs, q_sim)

# 均方根误差
rmse(q_obs, q_sim)

# 其他：mae, r_squared, person_r, nrmse_mean, ...
```

**充分性评估（adequacy）**

```julia
using Dates
dates = collect(Date(2000,1,1):Day(1):Date(2005,12,31))
result = evaluate_adequacy_assessment(q_sim, q_obs, precip, dates)
# result.adequate.RR  => 径流比是否充分
# result.adequate.Qmonth => 各月是否充分
```

### Raven 模板模型

```julia
model = build_raven_gr4j()
model = build_raven_hymod()
model = build_raven_hbv_light()

# 查看所有 Raven 模型
println(AVAILABLE_RAVEN_MODELS)
println(RAVEN_TEMPLATE_MODEL_STATUS)
```

---

## 可用模型列表

```
alpine1, alpine2, australia,
collie1, collie2, collie3,
cemaneigegr4j, echo, exphydro,
flexb, gr4j, gsfb, gsmsocont,
hbv_edu, hbv, hbv96,
hillslope, hmets, hymod,
ihacres, ihm19, lascam, mcrm,
modhydrolog, mopex1, mopex2, mopex3, mopex4, mopex5,
nam, newzealand1, newzealand2,
penman, plateau, prms, sacramento,
susannah1, susannah2, tank, tcm,
unitedstates, wetland, xaj
```

---

## 可用示例数据集

| 数据集 | 描述 |
|--------|------|
| `:cemaneige` | CemaNeige 融雪模型数据（温度、降水、融雪） |
| `:gr4j` | GR4J 模型数据（降水、蒸发、流量） |
| `:hbv_edu` | HBV 教育版模型数据 |
| `:hymod` | HYMOD 模型数据 |
| `:exphydro` | ExpHydro 数据（USGS 站点 01013500） |
| `:exphydro_alt` | ExpHydro 备用数据（USGS 站点 03604000） |
| `:marrmot` | MARRMOT 框架数据 |
| `:symhyd` | SYMHYD 模型数据 |
| `:m50` | M50 模型数据（USGS 站点 01013500） |

---

## 运行测试

```bash
# 在项目根目录
julia --project=. -e "using Pkg; Pkg.test()"
```

测试覆盖：示例数据加载、模型加载与参数接口、Raven 模型构建、充分性评估、路由算法。

---

## 项目结构

```
src/
├── HydroModelLibrary.jl   # 主模块
├── sample_data.jl         # 示例数据加载
├── models/                # 40+ 水文模型实现
├── fluxes/                # 35+ 通量过程模块
├── routes/                # 路由算法
├── criteria/              # 精度与充分性评估
└── raven-model/           # Raven 模板模型
data/                      # 示例数据（CSV）
test/                      # 测试套件
docs/                      # 文档
```
