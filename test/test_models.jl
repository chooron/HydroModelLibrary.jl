using Test
using HydroModelLibrary
using ComponentArrays

# 导入需要的函数
using HydroModelLibrary: load_model, get_params_bounds, get_random_params,
                         AVAILABLE_MODELS, AVAILABLE_RAVEN_MODELS,
                         build_raven_gr4j, build_raven_hymod, build_raven_hbv_light

@testset "Model Loading" begin
    @testset "AVAILABLE_MODELS" begin
        @test length(AVAILABLE_MODELS) > 0
        @test :gr4j in AVAILABLE_MODELS
        @test :hymod in AVAILABLE_MODELS
        @test :exphydro in AVAILABLE_MODELS
        @test :hbv_edu in AVAILABLE_MODELS
    end

    @testset "load_model" begin
        m = load_model(:gr4j)
        @test m isa Module
        @test isdefined(m, :model)
        @test isdefined(m, :model_parameters)

        # 字符串输入
        m2 = load_model("gr4j")
        @test m2 isa Module

        # 代表性模型都能加载（只测试确定存在的模型）
        # 注意：由于 world age 问题，我们只测试模块是否成功加载
        for nm in [:hymod, :exphydro, :hbv_edu]
            m = load_model(nm)
            @test m isa Module
        end
    end

    @testset "get_params_bounds" begin
        bounds = get_params_bounds(:gr4j)
        @test bounds isa NamedTuple
        @test length(bounds) > 0
        # 每个参数的 bounds 应是长度为 2 的元组
        for (k, v) in pairs(bounds)
            @test length(v) == 2
            @test v[1] <= v[2]
        end

        # 字符串输入
        bounds2 = get_params_bounds("gr4j")
        @test bounds2 isa NamedTuple
    end

    @testset "get_random_params" begin
        params = get_random_params(:gr4j)
        @test params isa ComponentVector
        @test haskey(params, :params)

        bounds = get_params_bounds(:gr4j)
        for k in keys(bounds)
            v = params.params[k]
            @test bounds[k][1] <= v <= bounds[k][2]
        end

        # 字符串输入
        params2 = get_random_params("hymod")
        @test params2 isa ComponentVector
    end
end

@testset "Raven Models" begin
    @test length(AVAILABLE_RAVEN_MODELS) > 0
    @test :gr4j in AVAILABLE_RAVEN_MODELS
    @test :hymod in AVAILABLE_RAVEN_MODELS
    @test :hbv_light in AVAILABLE_RAVEN_MODELS

    m = build_raven_gr4j()
    @test m !== nothing

    m = build_raven_hymod()
    @test m !== nothing

    m = build_raven_hbv_light()
    @test m !== nothing
end
