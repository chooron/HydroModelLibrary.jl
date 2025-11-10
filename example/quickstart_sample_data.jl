"""
快速入门：加载和使用示例数据

这是一个最简单的示例，展示如何快速开始使用 HydroModelLibrary 的示例数据。
"""

using HydroModelLibrary

println("="^70)
println("HydroModelLibrary - 示例数据快速入门")
println("="^70)

# 步骤 1：查看所有可用的数据集
println("\n【步骤 1】查看所有可用的示例数据集：")
println("-"^70)
list_sample_data()

# 步骤 2：加载一个示例数据集
println("\n【步骤 2】加载 GR4J 示例数据：")
println("-"^70)
data = load_sample_data(:gr4j)
println("✓ 数据加载成功！")
println("  数据维度: ", size(data))
println("  列名: ", names(data))

# 步骤 3：查看数据的前几行
println("\n【步骤 3】查看数据前 5 行：")
println("-"^70)
println(first(data, 5))

# 步骤 4：提取和使用数据
println("\n【步骤 4】数据分析示例：")
println("-"^70)
if "prec" in names(data)
    using Statistics
    prec = data[!, "prec"]
    println("降水统计：")
    println("  总数据点: ", length(prec))
    println("  平均值: ", round(mean(prec), digits=2))
    println("  最大值: ", round(maximum(prec), digits=2))
    println("  最小值: ", round(minimum(prec), digits=2))
end

# 完成
println("\n" * "="^70)
println("✓ 快速入门完成！")
println("="^70)
println("\n接下来您可以：")
println("  1. 尝试加载其他数据集: load_sample_data(:hymod)")
println("  2. 查看完整文档: docs/sample_data_guide.md")
println("  3. 运行完整示例: examples/load_sample_data_example.jl")
println("  4. 运行测试: include(\"test/test_sample_data.jl\")")
println("="^70)

