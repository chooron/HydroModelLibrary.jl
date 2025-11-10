"""
    示例数据加载模块

提供快速加载示例数据集的功能，帮助用户快速上手使用 HydroModelLibrary。
"""

using CSV
using DataFrames

"""
    AVAILABLE_SAMPLE_DATA

可用的示例数据集列表。每个数据集对应一个特定的水文模型。

# 可用数据集
- `:cemaneige` - CemaNeige 模型示例数据（温度、降水、融雪）
- `:gr4j` - GR4J 模型示例数据（降水、蒸发、流量）
- `:hbv_edu` - HBV 教育版模型示例数据
- `:hymod` - HYMOD 模型示例数据
- `:exphydro` - ExpHydro 模型示例数据（站点：01013500）
- `:exphydro_alt` - ExpHydro 模型备用数据（站点：03604000）
- `:marrmot` - MARRMOT 模型示例数据
- `:symhyd` - SYMHYD 模型示例数据
- `:m50` - M50 模型示例数据

# 示例
```julia
# 查看所有可用数据集
println(AVAILABLE_SAMPLE_DATA)

# 加载 GR4J 示例数据
data = load_sample_data(:gr4j)
```
"""
const AVAILABLE_SAMPLE_DATA = [
    :cemaneige,
    :gr4j,
    :hbv_edu,
    :hymod,
    :exphydro,
    :exphydro_alt,
    :marrmot,
    :symhyd,
    :m50
]

"""
    DATA_PATHS

内部使用：存储每个数据集的文件路径映射
"""
const DATA_PATHS = Dict(
    :cemaneige => joinpath("cemaneige", "sample.csv"),
    :gr4j => joinpath("gr4j", "sample.csv"),
    :hbv_edu => joinpath("hbv_edu", "hbv_sample.csv"),
    :hymod => joinpath("hymod", "sample.csv"),
    :exphydro => joinpath("exphydro", "01013500.csv"),
    :exphydro_alt => joinpath("exphydro", "03604000.csv"),
    :marrmot => joinpath("marrmot", "3604000.csv"),
    :symhyd => joinpath("symhyd", "sample.csv"),
    :m50 => joinpath("m50", "01013500.csv")
)

"""
    DATA_DESCRIPTIONS

内部使用：每个数据集的描述信息
"""
const DATA_DESCRIPTIONS = Dict(
    :cemaneige => "CemaNeige 融雪模型示例数据，包含温度、降水和融雪观测",
    :gr4j => "GR4J 模型示例数据，包含降水、蒸发和流量观测",
    :hbv_edu => "HBV 教育版模型示例数据",
    :hymod => "HYMOD 模型示例数据，包含气象和水文观测",
    :exphydro => "ExpHydro 模型示例数据（USGS 站点 01013500）",
    :exphydro_alt => "ExpHydro 模型备用示例数据（USGS 站点 03604000）",
    :marrmot => "MARRMOT 模型框架示例数据",
    :symhyd => "SYMHYD 模型示例数据",
    :m50 => "M50 模型示例数据（USGS 站点 01013500）"
)

"""
    load_sample_data(dataset_name::Symbol; data_dir=nothing)

加载指定的示例数据集。

# 参数
- `dataset_name::Symbol`: 数据集名称，可选值见 `AVAILABLE_SAMPLE_DATA`
- `data_dir`: 数据目录路径（可选），默认为包安装目录下的 data 文件夹

# 返回
- `DataFrame`: 加载的数据框

# 示例
```julia
# 加载 GR4J 示例数据
data = load_sample_data(:gr4j)

# 查看数据前几行
first(data, 5)

# 加载 CemaNeige 示例数据
cemaneige_data = load_sample_data(:cemaneige)

# 使用自定义数据目录
data = load_sample_data(:gr4j, data_dir="/path/to/data")
```

# 异常
- 如果数据集名称不存在，将抛出 `ArgumentError`
- 如果数据文件不存在，将抛出 `SystemError`
"""
function load_sample_data(dataset_name::Symbol; data_dir=nothing)
    # 检查数据集是否可用
    if !(dataset_name in AVAILABLE_SAMPLE_DATA)
        available_list = join(string.(AVAILABLE_SAMPLE_DATA), ", ")
        throw(ArgumentError(
            "数据集 '$dataset_name' 不可用。可用的数据集有: $available_list"
        ))
    end
    
    # 确定数据目录
    if isnothing(data_dir)
        # 默认使用包目录下的 data 文件夹
        package_dir = dirname(dirname(@__FILE__))
        data_dir = joinpath(package_dir, "data")
    end
    
    # 构建完整的文件路径
    relative_path = DATA_PATHS[dataset_name]
    full_path = joinpath(data_dir, relative_path)
    
    # 检查文件是否存在
    if !isfile(full_path)
        throw(SystemError(
            "数据文件不存在: $full_path"
        ))
    end
    
    # 加载 CSV 数据
    df = CSV.read(full_path, DataFrame)
    
    return df
end

# 重载：支持字符串输入
load_sample_data(dataset_name::String; data_dir=nothing) = 
    load_sample_data(Symbol(dataset_name); data_dir=data_dir)

"""
    list_sample_data()

列出所有可用的示例数据集及其描述信息。

# 示例
```julia
list_sample_data()
```
"""
function list_sample_data()
    println("="^70)
    println("HydroModelLibrary - 可用示例数据集")
    println("="^70)
    println()
    
    for dataset in AVAILABLE_SAMPLE_DATA
        desc = DATA_DESCRIPTIONS[dataset]
        path = DATA_PATHS[dataset]
        println("📊 :$dataset")
        println("   描述: $desc")
        println("   路径: data/$path")
        println()
    end
    
    println("="^70)
    println("使用方法: data = load_sample_data(:dataset_name)")
    println("示例: data = load_sample_data(:gr4j)")
    println("="^70)
end

"""
    get_sample_data_info(dataset_name::Symbol)

获取指定数据集的详细信息（不加载数据）。

# 参数
- `dataset_name::Symbol`: 数据集名称

# 返回
- `NamedTuple`: 包含数据集的描述、路径等信息

# 示例
```julia
info = get_sample_data_info(:gr4j)
println(info.description)
```
"""
function get_sample_data_info(dataset_name::Symbol)
    if !(dataset_name in AVAILABLE_SAMPLE_DATA)
        available_list = join(string.(AVAILABLE_SAMPLE_DATA), ", ")
        throw(ArgumentError(
            "数据集 '$dataset_name' 不可用。可用的数据集有: $available_list"
        ))
    end
    
    return (
        name = dataset_name,
        description = DATA_DESCRIPTIONS[dataset_name],
        path = DATA_PATHS[dataset_name],
        full_path = joinpath(dirname(dirname(@__FILE__)), "data", DATA_PATHS[dataset_name])
    )
end

get_sample_data_info(dataset_name::String) = get_sample_data_info(Symbol(dataset_name))

"""
    load_sample_data_for_model(model_name::Symbol; data_dir=nothing)

根据模型名称自动加载对应的示例数据。

# 参数
- `model_name::Symbol`: 模型名称（如 :gr4j, :hymod 等）
- `data_dir`: 数据目录路径（可选）

# 返回
- `DataFrame`: 加载的数据框

# 示例
```julia
# 自动加载 GR4J 模型的示例数据
data = load_sample_data_for_model(:gr4j)

# 自动加载 HYMOD 模型的示例数据
data = load_sample_data_for_model(:hymod)
```

# 注意
如果模型名称没有对应的示例数据，将尝试使用通用数据集。
"""
function load_sample_data_for_model(model_name::Symbol; data_dir=nothing)
    # 尝试直接匹配模型名称
    if model_name in AVAILABLE_SAMPLE_DATA
        return load_sample_data(model_name; data_dir=data_dir)
    end
    
    # 尝试常见的模型名称变体
    model_map = Dict(
        :hbv => :hbv_edu,
        :exphydro1 => :exphydro,
        :exphydro2 => :exphydro_alt
    )
    
    if haskey(model_map, model_name)
        dataset = model_map[model_name]
        @info "使用 :$dataset 作为 :$model_name 的示例数据"
        return load_sample_data(dataset; data_dir=data_dir)
    end
    
    # 如果没有找到对应的数据集，提供建议
    throw(ArgumentError(
        "模型 '$model_name' 没有专门的示例数据集。\n" *
        "请从以下可用数据集中选择: $(join(string.(AVAILABLE_SAMPLE_DATA), ", "))\n" *
        "使用方法: load_sample_data(:dataset_name)"
    ))
end

load_sample_data_for_model(model_name::String; data_dir=nothing) = 
    load_sample_data_for_model(Symbol(model_name); data_dir=data_dir)

export load_sample_data, list_sample_data, get_sample_data_info, 
       load_sample_data_for_model, AVAILABLE_SAMPLE_DATA

