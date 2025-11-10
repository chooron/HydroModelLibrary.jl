"""
测试示例数据加载功能
"""

using Test
using HydroModelLibrary
using DataFrames

@testset "示例数据加载测试" begin
    
    @testset "AVAILABLE_SAMPLE_DATA 常量" begin
        @test isdefined(HydroModelLibrary, :AVAILABLE_SAMPLE_DATA)
        @test length(AVAILABLE_SAMPLE_DATA) > 0
        @test :gr4j in AVAILABLE_SAMPLE_DATA
        @test :hymod in AVAILABLE_SAMPLE_DATA
    end
    
    @testset "load_sample_data 基本功能" begin
        # 测试加载 GR4J 数据
        data = load_sample_data(:gr4j)
        @test data isa DataFrame
        @test nrow(data) > 0
        @test ncol(data) > 0
        
        # 测试加载 HYMOD 数据
        data = load_sample_data(:hymod)
        @test data isa DataFrame
        @test nrow(data) > 0
        
        # 测试字符串输入
        data = load_sample_data("gr4j")
        @test data isa DataFrame
    end
    
    @testset "load_sample_data 错误处理" begin
        # 测试不存在的数据集
        @test_throws ArgumentError load_sample_data(:nonexistent_dataset)
    end
    
    @testset "get_sample_data_info" begin
        info = get_sample_data_info(:gr4j)
        @test info isa NamedTuple
        @test haskey(info, :name)
        @test haskey(info, :description)
        @test haskey(info, :path)
        @test haskey(info, :full_path)
        @test info.name == :gr4j
        
        # 测试字符串输入
        info = get_sample_data_info("gr4j")
        @test info isa NamedTuple
        
        # 测试不存在的数据集
        @test_throws ArgumentError get_sample_data_info(:nonexistent)
    end
    
    @testset "load_sample_data_for_model" begin
        # 测试直接匹配的模型
        data = load_sample_data_for_model(:gr4j)
        @test data isa DataFrame
        
        # 测试映射的模型名称
        data = load_sample_data_for_model(:hbv)
        @test data isa DataFrame
        
        # 测试字符串输入
        data = load_sample_data_for_model("gr4j")
        @test data isa DataFrame
        
        # 测试不存在映射的模型
        @test_throws ArgumentError load_sample_data_for_model(:unknown_model)
    end
    
    @testset "list_sample_data" begin
        # 测试函数不会抛出错误
        @test begin
            list_sample_data()
            true
        end
    end
    
    @testset "加载所有数据集" begin
        # 确保所有列出的数据集都可以加载
        for dataset in AVAILABLE_SAMPLE_DATA
            @testset "加载 $dataset" begin
                data = load_sample_data(dataset)
                @test data isa DataFrame
                @test nrow(data) > 0
                println("✓ $dataset: $(nrow(data)) 行, $(ncol(data)) 列")
            end
        end
    end
    
    @testset "数据质量检查" begin
        # 检查 ExpHydro 数据的关键列
        data = load_sample_data(:exphydro)
        expected_cols = ["date", "prcp(mm/day)", "tmean(C)", "flow(mm)"]
        for col in expected_cols
            @test col in names(data)
        end
        
        # 检查数据没有全为缺失值的列
        for col in names(data)
            if eltype(data[!, col]) <: Number
                @test !all(ismissing.(data[!, col]))
            end
        end
    end
end

println("\n" * "="^70)
println("✓ 所有测试通过！")
println("="^70)

