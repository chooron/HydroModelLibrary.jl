"""
示例数据加载使用示例

本文件展示如何使用 HydroModelLibrary 中的示例数据加载功能。
"""

using HydroModelLibrary

# ============================================================
# 1. 查看所有可用的示例数据集
# ============================================================
println("=" ^ 70)
println("1. 列出所有可用的示例数据集")
println("=" ^ 70)
list_sample_data()

# ============================================================
# 2. 加载特定的示例数据集
# ============================================================
println("\n" * "=" ^ 70)
println("2. 加载 GR4J 模型的示例数据")
println("=" ^ 70)

# 加载 GR4J 示例数据
gr4j_data = load_sample_data(:gr4j)
println("✓ 成功加载 GR4J 数据")
println("数据维度: ", size(gr4j_data))
println("列名: ", names(gr4j_data))
println("\n前 5 行数据:")
println(first(gr4j_data, 5))

# ============================================================
# 3. 加载其他模型的示例数据
# ============================================================
println("\n" * "=" ^ 70)
println("3. 加载 CemaNeige 模型的示例数据")
println("=" ^ 70)

cemaneige_data = load_sample_data(:cemaneige)
println("✓ 成功加载 CemaNeige 数据")
println("数据维度: ", size(cemaneige_data))
println("列名: ", names(cemaneige_data))
println("\n前 5 行数据:")
println(first(cemaneige_data, 5))

# ============================================================
# 4. 加载 HYMOD 模型的示例数据
# ============================================================
println("\n" * "=" ^ 70)
println("4. 加载 HYMOD 模型的示例数据")
println("=" ^ 70)

hymod_data = load_sample_data(:hymod)
println("✓ 成功加载 HYMOD 数据")
println("数据维度: ", size(hymod_data))
println("列名: ", names(hymod_data))
println("\n前 5 行数据:")
println(first(hymod_data, 5))

# ============================================================
# 5. 获取数据集的详细信息（不加载数据）
# ============================================================
println("\n" * "=" ^ 70)
println("5. 获取数据集详细信息")
println("=" ^ 70)

info = get_sample_data_info(:gr4j)
println("数据集名称: ", info.name)
println("描述: ", info.description)
println("相对路径: ", info.path)
println("完整路径: ", info.full_path)

# ============================================================
# 6. 根据模型名称自动加载数据
# ============================================================
println("\n" * "=" ^ 70)
println("6. 根据模型名称自动加载对应数据")
println("=" ^ 70)

# 自动加载 HBV 模型的数据
hbv_data = load_sample_data_for_model(:hbv)
println("✓ 为 HBV 模型加载了对应的示例数据")
println("数据维度: ", size(hbv_data))

# ============================================================
# 7. 实际应用示例：准备数据用于模型运行
# ============================================================
println("\n" * "=" ^ 70)
println("7. 实际应用示例：准备数据用于模型")
println("=" ^ 70)

# 加载数据
data = load_sample_data(:exphydro)
println("✓ 加载 ExpHydro 数据")

# 提取关键变量
if "prcp(mm/day)" in names(data)
    prcp = data[!, "prcp(mm/day)"]
    println("降水数据点数: ", length(prcp))
    println("平均降水: ", round(mean(prcp), digits=2), " mm/day")
end

if "tmean(C)" in names(data)
    tmean = data[!, "tmean(C)"]
    println("温度数据点数: ", length(tmean))
    println("平均温度: ", round(mean(tmean), digits=2), " °C")
end

if "flow(mm)" in names(data)
    flow = data[!, "flow(mm)"]
    println("流量数据点数: ", length(flow))
    println("平均流量: ", round(mean(flow), digits=2), " mm")
end

# ============================================================
# 8. 错误处理示例
# ============================================================
println("\n" * "=" ^ 70)
println("8. 错误处理示例")
println("=" ^ 70)

try
    # 尝试加载不存在的数据集
    data = load_sample_data(:nonexistent_data)
catch e
    println("✓ 捕获预期的错误:")
    println("错误类型: ", typeof(e))
    println("错误信息: ", e.msg)
end

# ============================================================
# 完成
# ============================================================
println("\n" * "=" ^ 70)
println("✓ 示例完成！")
println("=" ^ 70)
println("\n使用提示：")
println("- 使用 list_sample_data() 查看所有可用数据集")
println("- 使用 load_sample_data(:dataset_name) 加载数据")
println("- 使用 get_sample_data_info(:dataset_name) 获取数据集信息")
println("- 使用 load_sample_data_for_model(:model_name) 自动加载模型对应的数据")
println("=" ^ 70)

