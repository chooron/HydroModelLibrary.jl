using Test
using ComponentArrays
using HydroModels
using HydroModels: RouteIRF, SparseRouteConvolution, aggregate_route_kernel, build_irf_kernels

module RouteHarness
using HydroModels
using HydroModels: RouteIRF, SparseRouteConvolution, aggregate_route_kernel, build_irf_kernels
using ComponentArrays

include(joinpath(@__DIR__, "..", "src", "routes", "unithydro.jl"))
include(joinpath(@__DIR__, "..", "src", "routes", "channel_route.jl"))

using .unithydro: ROUTE_DUMP,
                  ROUTE_GAMMA_CONVOLUTION,
                  ROUTE_TRI_CONVOLUTION,
                  ROUTE_RESERVOIR_SERIES,
                  ROUTE_DIFFUSIVE_WAVE,
                  ROUTE_NONE,
                  ROUTE_PLUG_FLOW,
                  ROUTE_LINEAR_STORAGE,
                  ROUTE_STORAGE_COEFF,
                  ROUTE_MUSKINGUM,
                  ROUTE_NASH_CASCADE,
                  ROUTE_HYDROLOGIC,
                  TOC_MCDERMOTT_PILGRIM,
                  TOC_BRANSBY_WILLIAMS,
                  TOC_WILLIAMS_1922,
                  build_unit_hydrograph,
                  toc_mcdermott_pilgrim,
                  time_of_concentration
using .channel_route: build_channel_route,
                      ChannelRoute

export ROUTE_DUMP,
       ROUTE_GAMMA_CONVOLUTION,
       ROUTE_TRI_CONVOLUTION,
       ROUTE_RESERVOIR_SERIES,
       ROUTE_DIFFUSIVE_WAVE,
       ROUTE_NONE,
       ROUTE_PLUG_FLOW,
       ROUTE_LINEAR_STORAGE,
       ROUTE_STORAGE_COEFF,
       ROUTE_MUSKINGUM,
       ROUTE_NASH_CASCADE,
       ROUTE_HYDROLOGIC,
       TOC_MCDERMOTT_PILGRIM,
       TOC_BRANSBY_WILLIAMS,
       TOC_WILLIAMS_1922,
       build_unit_hydrograph,
       build_channel_route,
       ChannelRoute,
       toc_mcdermott_pilgrim,
       time_of_concentration
end

using .RouteHarness: ROUTE_DUMP,
                     ROUTE_GAMMA_CONVOLUTION,
                     ROUTE_TRI_CONVOLUTION,
                     ROUTE_RESERVOIR_SERIES,
                     ROUTE_DIFFUSIVE_WAVE,
                     ROUTE_NONE,
                     ROUTE_PLUG_FLOW,
                     ROUTE_LINEAR_STORAGE,
                     ROUTE_STORAGE_COEFF,
                     ROUTE_MUSKINGUM,
                     ROUTE_NASH_CASCADE,
                     ROUTE_HYDROLOGIC,
                     TOC_MCDERMOTT_PILGRIM,
                     TOC_BRANSBY_WILLIAMS,
                     TOC_WILLIAMS_1922,
                     build_unit_hydrograph,
                     build_channel_route,
                     ChannelRoute,
                     toc_mcdermott_pilgrim,
                     time_of_concentration

@testset "Routing" begin
    @variables q_in q_out
    @parameters alpha beta tp te n k
    @parameters reach_length wave_celerity diffusivity
    @parameters storage_time muskingum_k muskingum_x n_reservoirs
    @parameters hydro_k hydro_m

    impulse_short = reshape(Float64[1, 0, 0, 0, 0, 0], 1, :)
    impulse_long = zeros(Float64, 1, 40)
    impulse_long[1, 1] = 1.0

    @testset "In-Catchment Unit Hydrographs" begin
        dump = build_unit_hydrograph(ROUTE_DUMP; input=q_in, output=q_out)
        @test dump(impulse_short, ComponentVector(params=NamedTuple())) == impulse_short

        gamma_uh = build_unit_hydrograph(
            ROUTE_GAMMA_CONVOLUTION;
            input=q_in,
            output=q_out,
            alpha=alpha,
            beta=beta,
            max_lag=24,
        )
        gamma_out = gamma_uh(impulse_long, ComponentVector(params=(alpha=3.0, beta=1.5)))
        @test size(gamma_out) == size(impulse_long)
        @test all(gamma_out .>= 0)
        @test isapprox(sum(gamma_out), 1.0; atol=1.0e-3)

        gamma_tp = build_unit_hydrograph(
            ROUTE_GAMMA_CONVOLUTION;
            input=q_in,
            output=q_out,
            alpha=alpha,
            tp=tp,
            max_lag=24,
        )
        gamma_tp_out = gamma_tp(impulse_long, ComponentVector(params=(alpha=3.0, tp=2.0)))
        @test size(gamma_tp_out) == size(impulse_long)
        @test all(gamma_tp_out .>= 0)
        @test isapprox(sum(gamma_tp_out), 1.0; atol=1.0e-2)

        tri_uh = build_unit_hydrograph(
            ROUTE_TRI_CONVOLUTION;
            input=q_in,
            output=q_out,
            tp=tp,
            te=te,
        )
        tri_out = tri_uh(impulse_long, ComponentVector(params=(tp=2.0, te=5.0)))
        @test size(tri_out) == size(impulse_long)
        @test all(tri_out .>= 0)
        @test isapprox(sum(tri_out), 1.0; atol=1.0e-6)
        @test maximum(tri_out) > 0.0
        @test findmax(vec(tri_out))[2] in 2:3

        nash_uh = build_unit_hydrograph(
            ROUTE_RESERVOIR_SERIES;
            input=q_in,
            output=q_out,
            n=n,
            k=k,
            max_lag=24,
        )
        nash_out = nash_uh(impulse_long, ComponentVector(params=(n=3.0, k=0.9)))
        @test size(nash_out) == size(impulse_long)
        @test all(nash_out .>= 0)
        @test isapprox(sum(nash_out), 1.0; atol=2.0e-2)

        multi_input = zeros(Float64, 1, 2, 40)
        multi_input[1, 1, 1] = 1.0
        multi_input[1, 2, 1] = 1.0
        multi_gamma = build_unit_hydrograph(
            ROUTE_GAMMA_CONVOLUTION;
            input=q_in,
            output=q_out,
            alpha=alpha,
            beta=beta,
            max_lag=24,
            htypes=[1, 2],
        )
        multi_gamma_out = multi_gamma(
            multi_input,
            ComponentVector(params=(alpha=[2.0, 4.0], beta=[1.0, 2.0])),
        )
        @test size(multi_gamma_out) == size(multi_input)
        @test all(multi_gamma_out .>= 0)
        @test isapprox(sum(multi_gamma_out[1, 1, :]), 1.0; atol=2.0e-2)
        @test isapprox(sum(multi_gamma_out[1, 2, :]), 1.0; atol=2.0e-2)
        @test collect(multi_gamma_out[1, 1, :]) != collect(multi_gamma_out[1, 2, :])

        @test_throws ArgumentError build_unit_hydrograph(:UNKNOWN_ROUTE; input=q_in, output=q_out)
        @test_throws ArgumentError build_unit_hydrograph(ROUTE_GAMMA_CONVOLUTION; input=q_in, output=q_out, alpha=alpha)
    end


    @testset "Routing Kernel Layers" begin
        lag_irf = RouteIRF([:delay], (p, dt, horizon) -> begin
            weights = zeros(Float64, horizon)
            if horizon >= 2
                weights[2] = 1.0
            else
                weights[1] = 1.0
            end
            weights
        end; name=:lag1_irf)
        lag_kernels = build_irf_kernels(lag_irf, reshape([1.0, 1.0], 2, 1); delta_t=1.0, horizon=4)
        sparse_kernel = aggregate_route_kernel([0.0 0.0; 1.0 0.0], lag_kernels; topo_order=[1, 2], horizon=6)
        runoff = zeros(Float64, 2, 6)
        runoff[1, 1] = 1.0
        routed = SparseRouteConvolution(sparse_kernel)(runoff)
        @test isapprox(routed[1, 2], 1.0; atol=1.0e-8)
        @test isapprox(routed[2, 3], 1.0; atol=1.0e-8)
        @test all(routed .>= 0)
    end

    @testset "Time of Concentration" begin
        @test isapprox(toc_mcdermott_pilgrim(100.0), 0.031667 * 100.0^0.38; atol=1.0e-10)
        @test isapprox(time_of_concentration(TOC_MCDERMOTT_PILGRIM, 100.0), toc_mcdermott_pilgrim(100.0); atol=1.0e-10)
        @test isapprox(time_of_concentration(TOC_BRANSBY_WILLIAMS, 120.0; L=15.0, S=0.02), 0.0359 * 15.0 * 0.02^(-0.2) * 120.0^(-0.1); atol=1.0e-10)
        @test isapprox(time_of_concentration(TOC_WILLIAMS_1922, 120.0; L=15.0, S=0.02), 0.02539 * 15.0 * 0.02^(-0.2) * 120.0^(-0.1); atol=1.0e-10)
        @test_throws ArgumentError time_of_concentration(TOC_BRANSBY_WILLIAMS, 120.0; S=0.02)
        @test_throws ArgumentError time_of_concentration(TOC_WILLIAMS_1922, 120.0; L=15.0)
        @test_throws ArgumentError time_of_concentration(:UNKNOWN_TOC, 120.0)
    end

    @testset "Single-Reach Channel Routing" begin
        no_route = build_channel_route(ROUTE_NONE; input=q_in, output=q_out)
        @test no_route isa ChannelRoute
        @test no_route isa ChannelRoute
        @test no_route(impulse_short, ComponentVector(params=NamedTuple())) == impulse_short

        plug_route = build_channel_route(
            ROUTE_PLUG_FLOW;
            input=q_in,
            output=q_out,
            flow_length=reach_length,
            wave_celerity=wave_celerity,
            max_lag=6,
        )
        @test plug_route isa ChannelRoute
        plug_out = plug_route(impulse_short, ComponentVector(params=(reach_length=2.0, wave_celerity=1.0)))
        @test plug_out[1, 1] == 0.0
        @test plug_out[1, 2] == 0.0
        @test isapprox(plug_out[1, 3], 1.0; atol=1.0e-8)

        diff_route = build_channel_route(
            ROUTE_DIFFUSIVE_WAVE;
            input=q_in,
            output=q_out,
            flow_length=reach_length,
            wave_celerity=wave_celerity,
            diffusivity=diffusivity,
            max_lag=8,
        )
        @test diff_route isa ChannelRoute
        diff_out = diff_route(impulse_short, ComponentVector(params=(reach_length=2.0, wave_celerity=1.0, diffusivity=0.5)))
        @test size(diff_out) == size(impulse_short)
        @test all(diff_out .>= 0)
        @test 0.0 < sum(diff_out) <= 1.0

        linear_storage_route = build_channel_route(
            ROUTE_LINEAR_STORAGE;
            input=q_in,
            output=q_out,
            storage_time=storage_time,
            max_lag=8,
        )
        @test linear_storage_route isa ChannelRoute
        linear_storage_out = linear_storage_route(impulse_short, ComponentVector(params=(storage_time=1.0,)))
        @test isapprox(linear_storage_out[1, 1], 0.5; atol=1.0e-8)
        @test isapprox(linear_storage_out[1, 2], 0.25; atol=1.0e-8)
        @test all(linear_storage_out .>= 0)

        nash_cascade = build_channel_route(
            ROUTE_NASH_CASCADE;
            input=q_in,
            output=q_out,
            storage_time=storage_time,
            n_reservoirs=n_reservoirs,
            max_lag=12,
        )
        @test nash_cascade isa ChannelRoute
        nash_cascade_out = nash_cascade(impulse_short, ComponentVector(params=(storage_time=1.0, n_reservoirs=3.0)))
        @test size(nash_cascade_out) == size(impulse_short)
        @test all(nash_cascade_out .>= 0)
        @test nash_cascade_out[1, 1] < linear_storage_out[1, 1]
        @test sum(nash_cascade_out) < 1.0

        storage_route = build_channel_route(
            ROUTE_STORAGE_COEFF;
            input=q_in,
            output=q_out,
            storage_time=storage_time,
        )
        @test storage_route isa ChannelRoute
        storage_out = storage_route(impulse_short, ComponentVector(params=(storage_time=1.0,)))
        @test isapprox(storage_out[1, 1], 1 / 3; atol=1.0e-8)
        @test isapprox(storage_out[1, 2], 4 / 9; atol=1.0e-8)
        @test all(storage_out .>= 0)

        muskingum = build_channel_route(
            ROUTE_MUSKINGUM;
            input=q_in,
            output=q_out,
            muskingum_k=muskingum_k,
            muskingum_x=muskingum_x,
        )
        @test muskingum isa ChannelRoute
        muskingum_out = muskingum(impulse_short, ComponentVector(params=(muskingum_k=1.0, muskingum_x=0.0)))
        @test isapprox(muskingum_out, storage_out; atol=1.0e-8)

        hydrologic = build_channel_route(
            ROUTE_HYDROLOGIC;
            input=q_in,
            output=q_out,
            storage_coeff=hydro_k,
            storage_exponent=hydro_m,
        )
        @test hydrologic isa ChannelRoute
        hydrologic_out = hydrologic(impulse_short, ComponentVector(params=(hydro_k=1.0, hydro_m=1.2)))
        @test size(hydrologic_out) == size(impulse_short)
        @test all(isfinite, hydrologic_out)
        @test all(hydrologic_out .>= 0)

        @test_throws ArgumentError build_channel_route(ROUTE_PLUG_FLOW; input=q_in, output=q_out, flow_length=reach_length)
        @test_throws ArgumentError build_channel_route(ROUTE_DIFFUSIVE_WAVE; input=q_in, output=q_out, flow_length=reach_length, wave_celerity=wave_celerity)
        @test_throws ArgumentError build_channel_route(ROUTE_NASH_CASCADE; input=q_in, output=q_out, storage_time=storage_time)
        @test_throws ArgumentError build_channel_route(:UNKNOWN_ROUTE; input=q_in, output=q_out)

        bad_scalar_params = ComponentVector(params=(muskingum_k=[1.0, 2.0], muskingum_x=0.0))
        @test_throws ArgumentError muskingum(impulse_short, bad_scalar_params)
    end

    @testset "Network Channel Routing" begin
        adjacency = [0.0 0.0; 1.0 0.0]
        network_input = zeros(Float64, 1, 2, 6)
        network_input[1, 1, 1] = 1.0

        no_route = build_channel_route(ROUTE_NONE; input=q_in, output=q_out, adjacency=adjacency)
        no_route_out = no_route(network_input, ComponentVector(params=NamedTuple()))
        @test size(no_route_out) == size(network_input)
        @test no_route_out[1, 1, 1] == 1.0
        @test no_route_out[1, 2, 1] == 1.0

        network_plug = build_channel_route(
            ROUTE_PLUG_FLOW;
            input=q_in,
            output=q_out,
            adjacency=adjacency,
            flow_length=reach_length,
            wave_celerity=wave_celerity,
            max_lag=4,
        )
        network_plug_out = network_plug(
            network_input,
            ComponentVector(params=(reach_length=[1.0, 1.0], wave_celerity=[1.0, 1.0])),
        )
        @test size(network_plug_out) == size(network_input)
        @test isapprox(network_plug_out[1, 1, 2], 1.0; atol=1.0e-8)
        @test isapprox(network_plug_out[1, 2, 3], 1.0; atol=1.0e-8)

        network_linear_storage = build_channel_route(
            ROUTE_LINEAR_STORAGE;
            input=q_in,
            output=q_out,
            adjacency=adjacency,
            storage_time=storage_time,
            max_lag=8,
        )
        network_linear_storage_out = network_linear_storage(network_input, ComponentVector(params=(storage_time=[1.0, 1.0],)))
        @test size(network_linear_storage_out) == size(network_input)
        @test all(network_linear_storage_out .>= 0)
        @test network_linear_storage_out[1, 2, 2] > 0.0

        network_nash = build_channel_route(
            ROUTE_NASH_CASCADE;
            input=q_in,
            output=q_out,
            adjacency=adjacency,
            storage_time=storage_time,
            n_reservoirs=n_reservoirs,
            max_lag=12,
        )
        network_nash_out = network_nash(network_input, ComponentVector(params=(storage_time=[1.0, 1.0], n_reservoirs=[3.0, 3.0])))
        @test size(network_nash_out) == size(network_input)
        @test all(network_nash_out .>= 0)
        @test network_nash_out[1, 2, 3] > 0.0

        network_storage = build_channel_route(
            ROUTE_STORAGE_COEFF;
            input=q_in,
            output=q_out,
            adjacency=adjacency,
            storage_time=storage_time,
        )
        network_storage_out = network_storage(network_input, ComponentVector(params=(storage_time=[1.0, 1.0],)))
        @test size(network_storage_out) == size(network_input)
        @test all(network_storage_out .>= 0)
        @test network_storage_out[1, 2, 2] > 0.0

        network_hydrologic = build_channel_route(
            ROUTE_HYDROLOGIC;
            input=q_in,
            output=q_out,
            adjacency=adjacency,
            storage_coeff=hydro_k,
            storage_exponent=hydro_m,
        )
        network_hydrologic_out = network_hydrologic(network_input, ComponentVector(params=(hydro_k=[1.0, 1.0], hydro_m=[1.2, 1.2])))
        @test size(network_hydrologic_out) == size(network_input)
        @test all(isfinite, network_hydrologic_out)
        @test all(network_hydrologic_out .>= 0)

        @test_throws ArgumentError build_channel_route(ROUTE_NONE; input=q_in, output=q_out, adjacency=[0.0 1.0; 1.0 0.0])
    end
end











