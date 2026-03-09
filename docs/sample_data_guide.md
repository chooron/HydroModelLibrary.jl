# 示例数据加载指南

HydroModelLibrary 提供了便捷的示例数据加载功能，帮助用户快速上手和测试各种水文模型。

## 快速开始

### 1. 查看所有可用数据集

```julia
using HydroModelLibrary

# 列出所有可用的示例数据集
list_sample_data()
```

这将显示所有可用的数据集及其描述信息。

### 2. 加载示例数据

```julia
# 加载 GR4J 模型的示例数据
data = load_sample_data(:gr4j)

# 查看数据
first(data, 5)
```

### 3. 查看数据结构

```julia
# 获取数据维度
size(data)

# 查看列名
names(data)

# 查看数据摘要
describe(data)
```

## 可用的数据集

### 1. CemaNeige (`cemaneige`)
- **描述**: CemaNeige 融雪模型示例数据
- **包含**: 温度、降水、融雪观测数据
- **用途**: 测试融雪过程模拟

```julia
data = load_sample_data(:cemaneige)
```

**数据列**:
- `date`: 日期
- `min_temp`, `mean_temp`, `max_temp`: 最低、平均、最高温度
- `precipitation`: 降水量
- `liquid_outflow`: 液态径流
- `qsim`: 模拟流量
- `snowwater`: 雪水当量
- `thermal`: 热量

### 2. GR4J (`gr4j`)
- **描述**: GR4J 模型示例数据
- **包含**: 降水、蒸发、流量观测
- **用途**: 测试经典的 GR4J 降雨径流模型

```julia
data = load_sample_data(:gr4j)
```

**数据列**:
- `time`: 时间
- `pet`: 潜在蒸散发
- `prec`: 降水
- `qobs`: 观测流量
- `qsim`: 模拟流量
- 以及其他中间计算变量

### 3. HBV-EDU (`hbv_edu`)
- **描述**: HBV 教育版模型示例数据
- **包含**: HBV 模型所需的气象和水文数据
- **用途**: 学习和测试 HBV 模型

```julia
data = load_sample_data(:hbv_edu)
```

### 4. HYMOD (`hymod`)
- **描述**: HYMOD 模型示例数据
- **包含**: 气象和水文观测数据
- **用途**: 测试 HYMOD 概念性水文模型

```julia
data = load_sample_data(:hymod)
```

**数据列**:
- `isodate`: ISO 格式日期
- `precip`: 降水
- `tmax`, `tmin`: 最高、最低温度
- `wind`: 风速
- `pet`: 潜在蒸散发
- `q`: 观测流量
- 以及其他变量

### 5. ExpHydro (`exphydro`, `exphydro_alt`)
- **描述**: ExpHydro 模型示例数据
- **包含**: USGS 站点的气象和流量数据
- **用途**: 测试 ExpHydro 可解释性水文模型

```julia
# 站点 01013500
data = load_sample_data(:exphydro)

# 站点 03604000（备用）
data_alt = load_sample_data(:exphydro_alt)
```

**数据列**:
- `date`: 日期
- `prcp(mm/day)`: 降水量
- `tmean(C)`: 平均温度
- `dayl(day)`: 日长
- `srad(W/m2)`: 太阳辐射
- `vp(Pa)`: 水汽压
- `flow(mm)`: 流量

### 6. MARRMOT (`marrmot`)
- **描述**: MARRMOT 模型框架示例数据
- **包含**: 通用的水文模型输入数据

```julia
data = load_sample_data(:marrmot)
```

### 7. SYMHYD (`symhyd`)
- **描述**: SYMHYD 模型示例数据

```julia
data = load_sample_data(:symhyd)
```

### 8. M50 (`m50`)
- **描述**: M50 模型示例数据

```julia
data = load_sample_data(:m50)
```

## 高级用法

### 根据模型名称自动加载数据

```julia
# 自动为特定模型加载对应的示例数据
data = load_sample_data_for_model(:gr4j)
data = load_sample_data_for_model(:hbv)  # 自动加载 hbv_edu 数据
```

### 获取数据集信息（不加载数据）

```julia
info = get_sample_data_info(:gr4j)
println("描述: ", info.description)
println("路径: ", info.path)
```

### 使用自定义数据目录

```julia
# 如果数据存储在自定义位置
data = load_sample_data(:gr4j, data_dir="/path/to/your/data")
```

### 数据预处理示例

```julia
# 加载数据
data = load_sample_data(:exphydro)

# 提取特定列
precipitation = data[!, "prcp(mm/day)"]
temperature = data[!, "tmean(C)"]
flow = data[!, "flow(mm)"]

# 计算统计信息
using Statistics
println("平均降水: ", mean(precipitation), " mm/day")
println("平均温度: ", mean(temperature), " °C")
println("平均流量: ", mean(flow), " mm")

# 数据可视化
using Plots
plot(data[!, "date"], flow, label="流量", xlabel="日期", ylabel="流量 (mm)")
```

## 完整示例：运行模型

```julia
using HydroModelLibrary

# 1. 加载模型
model = load_model(:gr4j)

# 2. 加载示例数据
data = load_sample_data(:gr4j)

# 3. 获取参数范围
params_bounds = get_params_bounds(:gr4j)

# 4. 获取随机参数（或使用您自己的参数）
params = get_random_params(:gr4j)

# 5. 准备输入数据
# (根据具体模型需求提取相应的列)

# 6. 运行模型
# (根据模型的具体接口运行)
```

## 错误处理

```julia
# 如果数据集不存在，会抛出 ArgumentError
try
    data = load_sample_data(:nonexistent)
catch e
    println("错误: ", e.msg)
end

# 显示所有可用数据集
println("可用数据集: ", AVAILABLE_SAMPLE_DATA)
```

## API 参考

### `list_sample_data()`
列出所有可用的示例数据集及其描述。

### `load_sample_data(dataset_name::Symbol; data_dir=nothing)`
加载指定的示例数据集。

**参数**:
- `dataset_name`: 数据集名称（Symbol）
- `data_dir`: 可选的自定义数据目录

**返回**: DataFrame

### `get_sample_data_info(dataset_name::Symbol)`
获取数据集的详细信息（不加载数据）。

**返回**: NamedTuple (name, description, path, full_path)

### `load_sample_data_for_model(model_name::Symbol; data_dir=nothing)`
根据模型名称自动加载对应的示例数据。

### `AVAILABLE_SAMPLE_DATA`
常量：包含所有可用数据集名称的数组。

## 贡献

如果您想添加新的示例数据集，请：

1. 将数据文件放在 `data/` 目录下的适当子目录中
2. 在 `src/sample_data.jl` 中更新 `AVAILABLE_SAMPLE_DATA`、`DATA_PATHS` 和 `DATA_DESCRIPTIONS`
3. 更新此文档

## 许可

示例数据来自各种公开数据源，请注意各数据集的具体使用许可。

## 支持

如有问题，请提交 Issue 或查看项目文档。

