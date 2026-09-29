using Test
using Random
using Erebus

make_test_coordinates(Nx=33, Ny=33) = default_grid_coordinates()

@testset "Marker Container Architecture" begin
    @testset "Authoritative Field Set Equality" begin
        # 69 authoritative field names per PR 3c specification
        expected_fields = Set([
            :F_extract_m,
            :Fm,
            :XCm,
            :XH2Om,
            :XNm,
            :XSm,
            :XWsolidm,
            :XWsolidm0,
            :X_graphite_m,
            :X_ice_CH4_m,
            :X_ice_CO2_m,
            :X_ice_CO_m,
            :X_ice_H2O_m,
            :X_ice_H2S_m,
            :X_ice_N2_m,
            :X_ice_NH3_m,
            :X_ice_PH3_m,
            :X_refr_C_m,
            :X_refr_H_m,
            :X_refr_N_m,
            :X_refr_P_m,
            :X_refr_S_m,
            :Xfe_C_m,
            :Xfe_H_m,
            :Xfe_N_m,
            :Xfe_S_m,
            :Xfe_bulk,
            :Xfem,
            :Xfem0,
            :Xmin_cohenite_m,
            :Xmin_graphite_m,
            :Xmin_metal_matrix_m,
            :Xmin_nitride_m,
            :Xmin_schreibersite_m,
            :Xmin_troilite_m,
            :alphafluidcur,
            :alphasolidcur,
            :cohestotalm,
            :deltaIW_m,
            :etafluidcur_inv_kphim,
            :etatotalm,
            :etavpm,
            :fricttotalm,
            :hrtotalm,
            :inv_gggtotalm,
            :ktotalm,
            :nCH4_m,
            :nCO2_m,
            :nCO_m,
            :nC_graphite_m,
            :nFe0_m,
            :nFe2_m,
            :nFe3_m,
            :pfm0,
            :phim,
            :phinewm,
            :rhocptotalm,
            :rhofluidcur,
            :rhototalm,
            :sxxm,
            :sxym,
            :t_accreted,
            :tenstotalm,
            :tkm,
            :tkm_rhocptotalm,
            :tm,
            :w3d_m,
            :xm,
            :ym,
        ])

        # Instantiate full MarkerArrays containing all groups via override_config
        cfg = Erebus.override_config(
            default_config(),
            Dict(
                "coreformation.percolation_active" => true,
                "metal_partition.active" => true,
                "volatiles.active" => true,
                "redox.active" => true,
                "volatile_mixture.active" => true,
                "refractory.active" => true,
                "phase_tracking.active" => true,
                "accretion.active" => true,
            ),
        )

        coords = make_test_coordinates(33, 33)
        markers = init_marker_arrays(10, cfg, coords)

        container_fields = Set(all_marker_array_names(markers))
        @test container_fields == expected_fields
        @test length(container_fields) == 69
    end

    @testset "Deep Copy Isolation" begin
        cfg = Erebus.override_config(
            default_config(),
            Dict("coreformation.percolation_active" => true, "volatiles.active" => true),
        )
        coords = make_test_coordinates(33, 33)
        markers = init_marker_arrays(10, cfg, coords)

        markers.xm[1] = 123.456
        markers.Xfem[1] = 0.789
        markers_copied = copy(markers)

        # Mutate copied arrays
        markers_copied.xm[1] = 999.999
        markers_copied.Xfem[1] = 0.111

        # Assert original remains completely unchanged
        @test markers.xm[1] ≈ 123.456
        @test markers.Xfem[1] ≈ 0.789
        @test markers_copied.xm[1] ≈ 999.999
        @test markers_copied.Xfem[1] ≈ 0.111
    end

    @testset "Absent Group Preservation and Push Marker" begin
        # Inactive optional groups must be absent from groups NamedTuple
        cfg = Erebus.override_config(
            default_config(),
            Dict(
                "coreformation.percolation_active" => false,
                "coreformation.settling_active" => false,
                "metal_partition.active" => false,
                "thermodynamics.hr_fe" => false,
                "volatiles.active" => false,
                "magma_degassing.active" => false,
                "magma_transport.active" => false,
                "redox.active" => false,
                "volatile_mixture.active" => false,
                "refractory.active" => false,
                "phase_tracking.active" => false,
                "accretion.active" => false,
            ),
        )

        coords = make_test_coordinates(33, 33)
        markers = init_marker_arrays(5, cfg, coords)

        @test !haskey(markers.groups, :metal)
        @test !haskey(markers.groups, :volatiles)
        @test !haskey(markers.groups, :redox)
        @test !haskey(markers.groups, :hcnspo)
        @test !haskey(markers.groups, :phase)
        @test !haskey(markers.groups, :accretion)

        orig_len = length(markers)
        push_marker!(markers; xm=5000.0, ym=6000.0, tm=1, tkm=300.0, phim=0.1)

        @test length(markers) == orig_len + 1
        @test markers.xm[end] ≈ 5000.0
        @test markers.ym[end] ≈ 6000.0
        # Optional groups must remain absent
        @test !haskey(markers.groups, :metal)
        @test !haskey(markers.groups, :volatiles)
    end

    @testset "Type Stability and Group Access" begin
        cfg = Erebus.override_config(
            default_config(), Dict("coreformation.percolation_active" => true)
        )
        coords = make_test_coordinates(33, 33)
        markers = init_marker_arrays(10, cfg, coords)

        get_metal_group(m::MarkerArrays) = m.groups.metal
        @test @inferred(get_metal_group(markers)) isa MetalGroup
        @test length(markers.groups.metal.Xfem) == 10
    end

    @testset "A4 Defect Fix: Redox Array Length on Reseeding" begin
        cfg = Erebus.override_config(default_config(), Dict("redox.active" => true))
        coords = make_test_coordinates(33, 33)
        markers = init_marker_arrays(10, cfg, coords)

        @test haskey(markers.groups, :redox)
        @test length(markers.groups.redox.deltaIW_m) == 10
        @test length(markers.groups.redox.nFe0_m) == 10

        # Reseed marker array with expansion
        mdis = fill(1.0, coords.Ny, coords.Nx)
        mnum = zeros(Int, coords.Ny, coords.Nx)
        # Call container-based replenish
        marknum_new = replenish_markers!(markers, mdis, mnum; coords=coords, cfg=cfg)

        @test length(markers.groups.redox.deltaIW_m) == length(markers)
        @test length(markers.groups.redox.nFe0_m) == length(markers)
        @test length(markers.groups.redox.nC_graphite_m) == length(markers)
    end

    @testset "A5 Defect Fix: Accretion Retyping Initialization" begin
        cfg = Erebus.override_config(
            default_config(),
            Dict(
                "accretion.active" => true,
                "accretion.Xfe_bulk_accreted" => 0.18,
                "volatiles.active" => true,
                "metal_partition.active" => true,
            ),
        )
        coords = make_test_coordinates(33, 33)
        markers = init_marker_arrays(10, cfg, coords)

        # Set marker 1 as air
        markers.tm[1] = 3
        markers.xm[1] = coords.xcenter + 1000.0
        markers.ym[1] = coords.ycenter + 1000.0

        n_conv = advance_accretion_boundary!(
            0.0, 5000.0, markers; cfg=cfg, current_time=1.0e6
        )

        @test n_conv == 1
        @test markers.tm[1] == 2
        @test markers.Xfe_bulk[1] ≈ 0.18
        @test markers.Xfe_H_m[1] ≈ 0.0
        @test markers.Xfe_C_m[1] ≈ 0.0
    end

    @testset "Multi-Step Length Invariant" begin
        out_dir = mktempdir()
        cfg = Erebus.override_config(
            default_config(),
            Dict(
                "grid.Nx" => 17,
                "grid.Ny" => 17,
                "time.n_steps" => 5,
                "time.dt_initial" => 0.05,
                "time.dt_longest" => 0.05,
                "solver.p2m_mode" => :tiled,
                "solver.hydromech_solver" => :direct,
                "output.savematstep" => 100,
                "output.save_final" => false,
                "output.output_dir" => out_dir,
                "coreformation.percolation_active" => true,
                "volatiles.active" => true,
                "redox.active" => true,
                "accretion.active" => true,
            ),
        )

        sim_res = simulation_loop(cfg)
        markers = sim_res.markers
        expected_len = length(markers.xm)

        @test length(markers.ym) == expected_len
        @test length(markers.tm) == expected_len
        if haskey(markers.groups, :metal)
            @test length(markers.groups.metal.Xfem) == expected_len
        end
        if haskey(markers.groups, :volatiles)
            @test length(markers.groups.volatiles.XH2Om) == expected_len
        end
        if haskey(markers.groups, :redox)
            @test length(markers.groups.redox.deltaIW_m) == expected_len
        end
    end

    @testset "Container Iteration and Pairs" begin
        cfg = Erebus.override_config(
            default_config(),
            Dict(
                "volatiles.active" => true,
                "volatile_mixture.active" => true,
                "refractory.active" => true,
                "redox.active" => true,
            ),
        )
        coords = make_test_coordinates(33, 33)
        markers = init_marker_arrays(10, cfg, coords)

        # Test Base.values and Base.iterate on groups
        @test length(values(markers.groups.hcnspo)) == 13
        for prop in values(markers.groups.hcnspo)
            @test length(prop) == 10
        end
        for prop in markers.groups.hcnspo
            @test length(prop) == 10
        end

        # Test Base.pairs and Base.iterate on MarkerArrays
        pair_dict = Dict(pairs(markers))
        @test haskey(pair_dict, :xm)
        @test haskey(pair_dict, :ym)
        @test haskey(pair_dict, :XH2Om)
        @test length(pair_dict[:xm]) == 10

        # Test keys(markers)
        k_tuple = keys(markers)
        @test :xm in k_tuple
        @test :ym in k_tuple
        @test !(:core in k_tuple)
        @test !(:groups in k_tuple)
    end

    @testset "Accretion Group with track_accretion_time = false" begin
        cfg = Erebus.override_config(
            default_config(),
            Dict("accretion.active" => true, "accretion.track_accretion_time" => false),
        )
        coords = make_test_coordinates(33, 33)
        markers = init_marker_arrays(10, cfg, coords; initial_time=5.0e5)
        @test !haskey(markers.groups, :accretion)

        # With track_accretion_time = true, verify initial_time propagation
        cfg_tracked = Erebus.override_config(
            default_config(),
            Dict("accretion.active" => true, "accretion.track_accretion_time" => true),
        )
        markers_tracked = init_marker_arrays(10, cfg_tracked, coords; initial_time=5.0e5)
        @test haskey(markers_tracked.groups, :accretion)
        @test all(markers_tracked.groups.accretion.t_accreted .≈ 5.0e5)
    end

    @testset "Replenishment Redox and Graphite Inheritance" begin
        cfg = Erebus.override_config(
            default_config(), Dict("volatiles.active" => true, "redox.active" => true)
        )
        coords = make_test_coordinates(33, 33)
        markers = init_marker_arrays(10, cfg, coords)

        # Set specific donor values
        markers.groups.volatiles.X_graphite_m .= 0.05
        markers.groups.redox.deltaIW_m .= -1.5
        markers.groups.redox.nFe0_m .= 2.5

        mdis = fill(1.0, coords.Ny, coords.Nx)
        mnum = zeros(Int, coords.Ny, coords.Nx)
        replenish_markers!(markers, mdis, mnum; coords=coords, cfg=cfg)

        @test all(markers.groups.volatiles.X_graphite_m .≈ 0.05)
        @test all(markers.groups.redox.deltaIW_m .≈ -1.5)
        @test all(markers.groups.redox.nFe0_m .≈ 2.5)
    end
end
