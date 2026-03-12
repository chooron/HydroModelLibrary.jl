using Test
using Dates
using HydroModelLibrary

@testset "Adequacy Assessment" begin
    @test percentage_error(15.0, 10.0) == 50.0
    @test is_adequate(50.0)
    @test !is_adequate(50.01)

    dates = collect(Date(2020, 1, 1):Day(1):Date(2021, 12, 31))
    n = length(dates)

    q_obs = [10.0 + 2.0 * sin(2.0 * pi * i / 365.0) + 0.5 * cos(2.0 * pi * i / 30.0) for i in 1:n]
    q_sim = 1.2 .* q_obs
    p = fill(20.0, n)

    result = evaluate_adequacy_assessment(q_sim, q_obs, p, dates)

    @test haskey(result, :metadata)
    @test haskey(result, :simulated)
    @test haskey(result, :observed)
    @test haskey(result, :percentage_error)
    @test haskey(result, :adequate)

    @test result.metadata.frequency == :daily
    @test length(result.metadata.months) == 12
    @test length(result.metadata.years) == 2

    @test length(result.simulated.Qmonth) == 12
    @test length(result.observed.Qmonth) == 12
    @test length(result.simulated.Qyear) == length(result.metadata.years)
    @test length(result.observed.Qyear) == length(result.metadata.years)

    @test result.percentage_error.RR ≈ 20.0 atol = 1e-8
    @test result.adequate.RR
    @test all(result.adequate.Qmonth)
    @test all(result.adequate.Qyear)

    result_bad = evaluate_adequacy_assessment(2.0 .* q_obs, q_obs, p, dates)
    @test result_bad.percentage_error.RR ≈ 100.0 atol = 1e-8
    @test !result_bad.adequate.RR

    bad_dates = [Date(2020, 1, 1), Date(2020, 1, 3)]
    @test_throws ArgumentError evaluate_adequacy_assessment([1.0, 2.0], [1.0, 2.0], [1.0, 1.0], bad_dates)
end
