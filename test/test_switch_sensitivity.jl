# test/test_switch_sensitivity.jl
# 3-step sensitivity integration tests verifying that each active physics
# switch in SimulationConfig changes at least one output field when toggled,
# that all arrays remain finite, and that monotonic responses follow expected signs.

using Test
using Erebus
using Erebus.Config

@testset "Physics Switch Sensitivity" begin
    base_toml = normpath(
        joinpath(@__DIR__, "..", "configs", "magma_ocean_cooling_turb_on_32.toml")
    )
    base_cfg = load_config(base_toml)

    created_dirs = String[]
    function make_cfg(overrides::AbstractDict)
        out_dir = mktempdir()
        push!(created_dirs, out_dir)
        merged = Dict{String,Any}(
            "grid.Nx" => 17,
            "grid.Ny" => 17,
            "time.n_steps" => 3,
            "time.dt_initial" => 0.05,
            "time.dt_longest" => 0.05,
            "solver.p2m_mode" => :tiled,
            "solver.hydromech_solver" => :direct,
            "solver.seed" => 42,
            "output.mode" => :telemetry,
            "output.save_final" => false,
            "output.output_dir" => out_dir,
        )
        for (k, v) in overrides
            merged[k] = v
        end
        return Erebus.override_config(base_cfg, merged)
    end

    @testset "Magma Degassing Switch" begin
        cfg_on = make_cfg(
            Dict(
                "melting.active" => true,
                "volatiles.active" => true,
                "magma_degassing.active" => true,
                "magma_degassing.mode" => :dynamic_flux,
                "magma_degassing.degas_depth_fraction" => 0.90,
                "atmosphere.active" => true,
                "escape.active" => false,
            ),
        )
        cfg_off = make_cfg(
            Dict(
                "melting.active" => true,
                "volatiles.active" => true,
                "magma_degassing.active" => false,
                "atmosphere.active" => true,
                "escape.active" => false,
            ),
        )

        res_on = simulation_loop(cfg_on)
        res_off = simulation_loop(cfg_off)

        # Degassing on increases atmospheric elemental H mass
        @test res_on.atm.elem.H > res_off.atm.elem.H
        @test any(res_on.markers.XH2Om .!= res_off.markers.XH2Om)
        @test all(isfinite, res_on.markers.XH2Om)
        @test all(isfinite, res_on.grids.tk1)
    end

    @testset "Atmospheric Escape Switch" begin
        cfg_on = make_cfg(
            Dict(
                "melting.active" => true,
                "volatiles.active" => true,
                "magma_degassing.active" => true,
                "magma_degassing.mode" => :dynamic_flux,
                "atmosphere.active" => true,
                "escape.active" => true,
                "escape.multi_species" => true,
            ),
        )
        cfg_off = make_cfg(
            Dict(
                "melting.active" => true,
                "volatiles.active" => true,
                "magma_degassing.active" => true,
                "magma_degassing.mode" => :dynamic_flux,
                "atmosphere.active" => true,
                "escape.active" => false,
            ),
        )

        res_on = simulation_loop(cfg_on)
        res_off = simulation_loop(cfg_off)

        # Escape on reduces retained atmospheric H mass and increases escaped H
        @test res_on.atm.elem.H < res_off.atm.elem.H
        @test res_on.atm.escaped.H > res_off.atm.escaped.H
        @test all(isfinite, res_on.grids.tk1)
        @test all(isfinite, res_on.markers.tkm)
    end

    @testset "Radiogenic Heating Switch" begin
        cfg_on = make_cfg(
            Dict(
                "thermodynamics.hr_al" => true,
                "thermodynamics.ratio_al" => 1.0e-4,
                "melting.active" => true,
            ),
        )
        cfg_off = make_cfg(Dict("thermodynamics.hr_al" => false, "melting.active" => true))

        res_on = simulation_loop(cfg_on)
        res_off = simulation_loop(cfg_off)

        # Radiogenic heating on raises mean rock temperature
        @test sum(res_on.markers.tkm) > sum(res_off.markers.tkm)
        @test any(res_on.grids.tk1 .!= res_off.grids.tk1)
        @test all(isfinite, res_on.grids.tk1)
        @test all(isfinite, res_on.markers.tkm)
    end

    @testset "Melting Switch" begin
        cfg_on = make_cfg(Dict("melting.active" => true))
        cfg_off = make_cfg(
            Dict("melting.active" => false, "melting.soft_turbulence" => false)
        )

        res_on = simulation_loop(cfg_on)
        res_off = simulation_loop(cfg_off)

        @test any(res_on.markers.Fm .!= res_off.markers.Fm)
        @test all(isfinite, res_on.markers.Fm)
        @test all(isfinite, res_on.grids.tk1)
    end

    @testset "Volatiles Switch" begin
        cfg_on = make_cfg(
            Dict(
                "volatiles.active" => true,
                "volatiles.initial_water_wtpct" => 1.0,
                "melting.active" => true,
            ),
        )
        cfg_off = make_cfg(Dict("volatiles.active" => false, "melting.active" => true))

        res_on = simulation_loop(cfg_on)
        res_off = simulation_loop(cfg_off)

        @test haskey(res_on.markers, :XH2Om)
        @test !haskey(res_off.markers, :XH2Om)
        @test all(isfinite, res_on.markers.XH2Om)
        @test all(isfinite, res_on.grids.tk1)
    end

    @testset "Venting Switch" begin
        cfg_on = make_cfg(
            Dict(
                "melting.active" => true,
                "volatiles.active" => true,
                "venting.active" => true,
                "poroelasticity.hydrofracture" => true,
                "retention.active" => true,
                "retention.venting_drainage_active" => true,
            ),
        )
        cfg_off = make_cfg(
            Dict(
                "melting.active" => true,
                "volatiles.active" => true,
                "venting.active" => false,
                "poroelasticity.hydrofracture" => true,
                "retention.active" => false,
                "retention.venting_drainage_active" => false,
            ),
        )

        res_on = simulation_loop(cfg_on)
        res_off = simulation_loop(cfg_off)

        @test any(res_on.grids.pf .!= res_off.grids.pf)
        @test all(isfinite, res_on.grids.pf)
        @test all(isfinite, res_on.markers.phim)
    end

    @testset "Metal Partition Switch" begin
        cfg_on = make_cfg(
            Dict(
                "melting.active" => true,
                "volatiles.active" => true,
                "metal_partition.active" => true,
            ),
        )
        cfg_off = make_cfg(
            Dict(
                "melting.active" => true,
                "volatiles.active" => true,
                "metal_partition.active" => false,
            ),
        )

        res_on = simulation_loop(cfg_on)
        res_off = simulation_loop(cfg_off)

        @test haskey(res_on.markers, :Xfe_H_m)
        @test !haskey(res_off.markers, :Xfe_H_m)
        @test all(isfinite, res_on.markers.Xfe_H_m)
        @test all(isfinite, res_on.grids.tk1)
    end

    @testset "Redox Switch" begin
        cfg_on = make_cfg(
            Dict(
                "melting.active" => true, "volatiles.active" => true, "redox.active" => true
            ),
        )
        cfg_off = make_cfg(
            Dict(
                "melting.active" => true,
                "volatiles.active" => true,
                "redox.active" => false,
            ),
        )

        res_on = simulation_loop(cfg_on)
        res_off = simulation_loop(cfg_off)

        @test haskey(res_on.markers, :deltaIW_m)
        @test !haskey(res_off.markers, :deltaIW_m)
        @test all(isfinite, res_on.markers.deltaIW_m)
        @test all(isfinite, res_on.grids.tk1)
    end

    for d in created_dirs
        rm(d; recursive=true, force=true)
    end
end
