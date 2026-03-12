# 示例数据加载指南

HydroModelLibrary 提供了便捷的示例数据加载功能，帮助用户快速上手和测试各种水文模型。

## 快速开始

```julia
using HydroModelLibrary

# 查看所有可用数据集
list_sample_data()

# 加载 GR4J 示例数据
data = load_sample_data(:gr4j)
first(data, 5)
```

---

## 可用数据集

`AVAILABLE_SAMPLE_DATA` 常量列出所有可用数据集名称：

```julia
println(AVAILABLE_SAMPLE_DATA)
# [:cemaneige, :gr4j, :hbv_edu, :hymod, :exphydro, :exphydro_alt, :marrmot, :symhyd, :m50]
```

### 1. CemaNeige (`:cemaneige`)

CemaNeige 融雪模型示例数据，包含温度、降水和融雪观测。

```julia
data = load_sample_data(:cemaneige)
```

主要列：`date`, `min_temp`, `mean_temp`, `max_temp`, `precipitation`, `liquid_outflow`, `qsim`, `snowwater`, `thermal`

### 2. GR4J (`:gr4j`)

GR4J 模型示例数据，包含降水、蒸发和流量观测。

```julia
data = load_sample_data(:gr4j)
```

主要列：`time`, `pet`, `prec`, `qobs`, `qsim` 及中间计算变量

### 3. HBV-EDU (`:hbv_edu`)

HBV 教育版模型示例数据。

```julia
data = load_sample_data(:hbv_edu)
```

### 4. HYMOD (`:hymod`)

HYMOD 模型示例数据，包含气象和水文观测。

```julia
data = load_sample_data(:hymod)
```

主要列：`isodate`, `precip`, `tmax`, `tmin`, `wind`, `pet`, `q`

### 5. ExpHydro (`:exphydro`, `:exphydro_alt`)

ExpHydro 模型示例数据，来自 USGS 站点。

```julia
data     = load_sample_data(:exphydro)      # 站点 01013500
data_alt = load_sample_data(:exphydro_alt)  # 站点 03604000
```

主要列：`date`, `prcp(mm/day)`, `tmean(C)`, `dayl(day)`, `srad(W/m2)`, `vp(Pa)`, `flow(mm)`

### 6. MARRMOT (`:marrmot`)

MARRMOT 模型框架示例数据。

```julia
data = load_sample_data(:marrmot)
```

### 7. SYMHYD (`:symhyd`)

SYMHYD 模型示例数据。

```julia
data = load_sample_data(:symhyd)
```

### 8. M50 (`:m50`)

M50 模型示例数据（USGS 站点 01013500）。

```julia
data = load_sample_data(:m50)
```

---

## API 参考

### `list_sample_data()`

打印所有可用数据集及其描述和路径。

### `load_sample_data(dataset_name; data_dir=nothing)`

加载指定数据集，返回 `DataFrame`。

- `dataset_name`：`Symbol` 或 `String`，须在 `AVAILABLE_SAMPLE_DATA` 中
- `data_dir`：可选，自定义数据目录路径（默认使用包内 `data/` 目录）

```julia
data = load_sample_data(:gr4j)
data = load_sample_data("gr4j")                        # 字符串也可以
data = load_sample_data(:gr4j, data_dir="/my/data")    # 自定义路径
```

### `load_sample_data_for_model(model_name; data_dir=nothing)`

按模型名自动匹配并加载对应数据集。支持以下映射：

| 模型名 | 映射数据集 |
|--------|-----------|
| `:gr4j` | `:gr4j` |
| `:hymod` | `:hymod` |
| `:hbv` | `:hbv_edu` |
| `:exphydro1` | `:exphydro` |
| `:exphydro2` | `:exphydro_alt` |

```julia
data = load_sample_data_for_model(:gr4j)
data = load_sample_data_for_model(:hbv)   # 自动加载 :hbv_edu 数据
```

### `get_sample_data_info(dataset_name)`

获取数据集元信息，不加载数据，返回 `NamedTuple`。

```julia
info = get_sample_data_info(:gr4j)
# (name=:gr4j, description="...", path="gr4j/sample.csv", full_path="...")
```

---

## 完整示例：加载数据并运行模型

```julia
using HydroModelLibrary
using ComponentArrays
using HydroModels

# 加载模型和数据
m = load_model(:exphydro)
data = load_sample_data(:exphydro)

# 准备输入
input_names = HydroModels.get_input_names(m.model)
input = (prcp=data[!, "prcp(mm/day)"], temp=data[!, "tmean(C)"])
input_arr = stack(input[input_names], dims=1)

# 初始化状态和参数
state_names = HydroModels.get_state_names(m.model)
init_states = ComponentVector(NamedTuple{Tuple(state_names)}(zeros(length(state_names))))
params = get_random_params(:exphydro)

# 运行模型
result = m.model(input_arr, params, initstates=init_states)
```

---

## 错误处理

```julia
# 数据集不存在时抛出 ArgumentError
try
    load_sample_data(:nonexistent)
catch e
    println(e.msg)
end
```

---

## 添加新数据集

1. 将 CSV 文件放入 `data/<dataset_name>/` 目录
2. 在 `src/sample_data.jl` 中更新 `AVAILABLE_SAMPLE_DATA`、`DATA_PATHS`、`DATA_DESCRIPTIONS`
3. 更新本文档
