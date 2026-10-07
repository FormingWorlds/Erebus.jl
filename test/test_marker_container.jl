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

    @testset "Container Operations, Indexing, and Serialization" begin
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

        # Indexing and key checks
        @test markers[:xm] === markers.core.xm
        @test markers["ym"] === markers.core.ym
        @test haskey(markers, "xm")
        @test haskey(markers, :XH2Om)
        @test !haskey(markers, "nonexistent")
        @test !haskey(markers, :nonexistent)
        @test_throws ErrorException markers.nonexistent_field

        # propertynames
        pnames = propertynames(markers)
        @test :core in pnames
        @test :groups in pnames
        @test :xm in pnames
        @test :deltaIW_m in pnames

        # Group conversions and operations
        for grp in values(markers.groups)
            nt = NamedTuple(grp)
            @test nt isa NamedTuple
            @test length(nt) == length(fieldnames(typeof(grp)))
            merged = merge((; test_entry=1), grp)
            @test haskey(merged, :test_entry)
            @test keys(grp) == fieldnames(typeof(grp))
        end

        # Base.resize!
        resize!(markers, 20)
        @test length(markers) == 20
        @test length(markers.core.xm) == 20
        @test length(markers.groups.metal.Xfem) == 20
        @test length(markers.groups.volatiles.XH2Om) == 20
        @test length(markers.groups.redox.deltaIW_m) == 20

        # Iteration, keys, and pairs
        @test :xm in keys(markers)
        @test :deltaIW_m in keys(markers)
        pairs_list = collect(pairs(markers))
        @test any(p -> p.first === :xm, pairs_list)
        iter_count = 0
        for (k, v) in markers
            iter_count += 1
            @test length(v) == 20
        end
        @test iter_count == length(keys(markers))

        # Base.copy
        copied = copy(markers)
        @test length(copied) == 20
        @test copied.core.xm ≈ markers.core.xm
        @test copied.groups.metal.Xfem ≈ markers.groups.metal.Xfem
        copied.core.xm[1] += 999.0
        @test markers.core.xm[1] != copied.core.xm[1]

        # push_marker! with groups present
        push_marker!(markers; xm=1.0, ym=2.0, tm=2, XH2Om=0.5, deltaIW_m=-1.0)
        @test length(markers) == 21
        @test markers.xm[21] ≈ 1.0
        @test markers.ym[21] ≈ 2.0
        @test markers.tm[21] == 2
        @test markers.groups.volatiles.XH2Om[21] ≈ 0.5
        @test markers.groups.redox.deltaIW_m[21] ≈ -1.0

        # serialize_marker_arrays
        serialized = serialize_marker_arrays(markers)
        @test serialized isa Dict{String,Any}
        @test haskey(serialized, "xm")
        @test haskey(serialized, "deltaIW_m")
        @test length(serialized["xm"]) == 21

        # restore_marker_arrays!
        target = init_marker_arrays(5, cfg, coords)
        @test length(target) == 5
        restore_marker_arrays!(target, markers)
        @test length(target) == 21
        @test target.core.xm ≈ markers.core.xm
        @test target.groups.metal.Xfem ≈ markers.groups.metal.Xfem
    end

    @testset "Configuration Variations in init_marker_arrays" begin
        coords = make_test_coordinates(33, 33)

        # Metal group without volatile metal partitioning
        cfg_fe = Erebus.override_config(
            default_config(),
            Dict(
                "thermodynamics.hr_fe" => true,
                "metal_partition.active" => false,
                "coreformation.percolation_active" => false,
                "coreformation.settling_active" => false,
            ),
        )
        markers_fe = init_marker_arrays(8, cfg_fe, coords)
        @test haskey(markers_fe.groups, :metal)
        @test length(markers_fe.groups.metal.Xfem) == 8
        @test all(markers_fe.groups.metal.Xfe_H_m .≈ 0.0)
        @test all(markers_fe.groups.metal.Xfe_C_m .≈ 0.0)

        # Volatiles group without degassing or volatiles active (transport only)
        cfg_trans = Erebus.override_config(
            default_config(),
            Dict(
                "melting.active" => true,
                "magma_transport.active" => true,
                "volatiles.active" => false,
                "magma_degassing.active" => false,
            ),
        )
        markers_trans = init_marker_arrays(8, cfg_trans, coords)
        @test haskey(markers_trans.groups, :volatiles)
        @test all(markers_trans.groups.volatiles.XH2Om .≈ 0.0)
        @test all(markers_trans.groups.volatiles.XCm .≈ 0.0)

        # Redox group without metal group (initial_xfe_bulk is nothing)
        cfg_rdx = Erebus.override_config(
            default_config(),
            Dict(
                "redox.active" => true,
                "coreformation.percolation_active" => false,
                "coreformation.settling_active" => false,
                "metal_partition.active" => false,
                "thermodynamics.hr_fe" => false,
            ),
        )
        markers_rdx = init_marker_arrays(8, cfg_rdx, coords)
        @test haskey(markers_rdx.groups, :redox)
        @test !haskey(markers_rdx.groups, :metal)
        @test length(markers_rdx.groups.redox.deltaIW_m) == 8
        @test length(markers_rdx.groups.redox.nFe0_m) == 8
    end

    @testset "Accretion Boundary Full Groups and Disk States" begin
        coords = make_test_coordinates(33, 33)
        cfg = Erebus.override_config(
            default_config(),
            Dict(
                "accretion.active" => true,
                "accretion.Xfe_bulk_accreted" => 0.15,
                "metal_partition.active" => true,
                "volatiles.active" => true,
                "redox.active" => true,
                "volatile_mixture.active" => true,
                "refractory.active" => true,
                "phase_tracking.active" => true,
            ),
        )
        markers = init_marker_arrays(10, cfg, coords)

        # Place markers 1 and 2 inside accretion envelope as sticky air (tm = 3)
        markers.tm[1] = 3
        markers.xm[1] = coords.xcenter + 500.0
        markers.ym[1] = coords.ycenter + 500.0
        markers.tm[2] = 3
        markers.xm[2] = coords.xcenter + 600.0
        markers.ym[2] = coords.ycenter + 600.0

        disk_state_wet = (;
            condensed_H2O=true,
            X_ice_H2O=0.12,
            X_ice_NH3=0.01,
            X_ice_CO2=0.02,
            X_ice_CO=0.005,
            X_ice_CH4=0.002,
            X_ice_N2=0.001,
            X_ice_H2S=0.003,
            X_ice_PH3=0.0001,
            f_refr_C=0.02,
            f_refr_S=0.015,
            f_refr_N=0.001,
            f_refr_P=0.0005,
            f_refr_H=0.0002,
        )

        n_conv1 = advance_accretion_boundary!(
            0.0, 1000.0, markers; cfg=cfg, disk_state=disk_state_wet, current_time=1.5e6
        )
        @test n_conv1 == 2
        @test markers.tm[1] == 2
        @test markers.tm[2] == 2
        @test markers.groups.metal.Xfe_bulk[1] ≈ 0.15
        @test markers.groups.volatiles.XH2Om[1] ≈ cfg.accretion.XH2O_wet_wtpct
        # deltaIW for phi_fe=0.15 via two-phase mass conversion (w_fe = 7/24)
        @test isapprox(markers.groups.redox.deltaIW_m[1], -1.088035858059; atol=1e-6)
        @test markers.groups.hcnspo.X_ice_H2O_m[1] ≈ 0.12
        @test markers.groups.phase.Xmin_troilite_m[1] ≈ 0.0
        @test markers.groups.accretion.t_accreted[1] ≈ 1.5e6

        # Dry disk state branch
        markers.tm[3] = 3
        markers.xm[3] = coords.xcenter + 700.0
        markers.ym[3] = coords.ycenter + 700.0
        disk_state_dry = (;
            condensed_H2O=false,
            X_ice_H2O=0.0,
            X_ice_NH3=0.0,
            X_ice_CO2=0.0,
            X_ice_CO=0.0,
            X_ice_CH4=0.0,
            X_ice_N2=0.0,
            X_ice_H2S=0.0,
            X_ice_PH3=0.0,
            f_refr_C=0.02,
            f_refr_S=0.015,
            f_refr_N=0.001,
            f_refr_P=0.0005,
            f_refr_H=0.0002,
        )
        n_conv2 = advance_accretion_boundary!(
            0.0, 1000.0, markers; cfg=cfg, disk_state=disk_state_dry, current_time=2.0e6
        )
        @test n_conv2 == 1
        @test markers.tm[3] == 2
        @test markers.groups.volatiles.XH2Om[3] ≈ cfg.accretion.XH2O_dry_wtpct

        # Advance without cfg (fallback defaults)
        markers.tm[4] = 3
        markers.xm[4] = coords.xcenter + 500.0
        markers.ym[4] = coords.ycenter + 500.0
        n_conv3 = advance_accretion_boundary!(
            0.0, 1000.0, markers; cfg=nothing, current_time=3.0e6
        )
        @test n_conv3 == 1
        @test markers.tm[4] == 2
        @test markers.tkm[4] ≈ 200.0
    end

    @testset "Checkpoint Dictionary Collection and State Reconstruction" begin
        coords = make_test_coordinates(33, 33)
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
        markers = init_marker_arrays(10, cfg, coords)

        # 1. _collect_checkpoint_marker_dict with MarkerArrays
        d_ma = Erebus._collect_checkpoint_marker_dict(
            markers;
            core_budgets=(; core_fe=1.0),
            M_atm_species=Dict("H2O" => 1.0),
            M_escaped_species=Dict("H2O" => 0.1),
            M_accreted_total=5.0e18,
            M_planet_val=6.0e19,
            regional_mineral_modes=(; troilite=0.01),
        )
        @test d_ma isa Dict{Symbol,Any}
        @test haskey(d_ma, :xm)
        @test haskey(d_ma, :deltaIW_m)
        @test d_ma[:M_accreted_total] ≈ 5.0e18
        @test d_ma[:M_planet_val] ≈ 6.0e19

        # 2. _collect_checkpoint_marker_dict with individual kwargs (fallback branch)
        N = 10
        d_raw = Erebus._collect_checkpoint_marker_dict(
            nothing;
            xm=zeros(N),
            ym=zeros(N),
            tm=ones(Int, N),
            tkm=fill(300.0, N),
            sxxm=zeros(N),
            sxym=zeros(N),
            etavpm=fill(1.0e20, N),
            phim=zeros(N),
            rhototalm=fill(3000.0, N),
            rhocptotalm=fill(3.0e6, N),
            etatotalm=fill(1.0e20, N),
            hrtotalm=zeros(N),
            ktotalm=fill(3.0, N),
            tkm_rhocptotalm=fill(9.0e8, N),
            etafluidcur_inv_kphim=zeros(N),
            inv_gggtotalm=fill(1.0e-10, N),
            fricttotalm=zeros(N),
            cohestotalm=zeros(N),
            tenstotalm=zeros(N),
            rhofluidcur=fill(1000.0, N),
            alphasolidcur=fill(3.0e-5, N),
            alphafluidcur=fill(1.0e-4, N),
            XWsolidm0=fill(0.1, N),
            F_extract_m=zeros(N),
            Xfem=zeros(N),
            Xfem0=zeros(N),
            Xfe_bulk=fill(0.2, N),
            XH2Om=fill(1.0, N),
            XCm=fill(100.0, N),
            XNm=fill(10.0, N),
            XSm=fill(1000.0, N),
            Xfe_H_m=zeros(N),
            Xfe_C_m=zeros(N),
            Xfe_N_m=zeros(N),
            Xfe_S_m=zeros(N),
            Xmin_troilite_m=zeros(N),
            Xmin_schreibersite_m=zeros(N),
            Xmin_cohenite_m=zeros(N),
            Xmin_nitride_m=zeros(N),
            Xmin_metal_matrix_m=zeros(N),
            Xmin_graphite_m=zeros(N),
            X_graphite_m=zeros(N),
            t_accreted=zeros(N),
            hcnspo_props=(; X_ice_H2O_m=zeros(N)),
            redox_props=(; deltaIW_m=fill(-2.0, N)),
            core_budgets=(; core_fe=1.0),
            M_atm_species=Dict("H2O" => 1.0),
            M_escaped_species=Dict("H2O" => 0.1),
            M_accreted_total=5.0e18,
            M_planet_val=6.0e19,
            regional_mineral_modes=(; troilite=0.01),
        )
        @test d_raw isa Dict{Symbol,Any}
        @test haskey(d_raw, :xm)
        @test haskey(d_raw, :Xmin_troilite_m)
        @test haskey(d_raw, :Xfe_H_m)
        @test d_raw[:Xfe_bulk][1] ≈ 0.2

        # 3. _reconstruct_checkpoint_marker_arrays: full, partial, minimal
        core = markers.core
        opt_full = (;
            Xfem=zeros(N),
            Xfem0=zeros(N),
            Xfe_bulk=fill(0.2, N),
            Xfe_H_m=zeros(N),
            Xfe_C_m=zeros(N),
            Xfe_N_m=zeros(N),
            Xfe_S_m=zeros(N),
            XH2Om=fill(1.0, N),
            XCm=fill(100.0, N),
            XNm=fill(10.0, N),
            XSm=fill(1000.0, N),
            X_graphite_m=zeros(N),
            F_extract_m=zeros(N),
            redox_props=(;
                nFe0_m=zeros(N),
                nFe2_m=zeros(N),
                nFe3_m=zeros(N),
                deltaIW_m=fill(-2.0, N),
                nC_graphite_m=zeros(N),
                nCO_m=zeros(N),
                nCO2_m=zeros(N),
                nCH4_m=zeros(N),
            ),
            hcnspo_props=(;
                X_ice_H2O_m=zeros(N),
                X_ice_NH3_m=zeros(N),
                X_ice_CO2_m=zeros(N),
                X_ice_CO_m=zeros(N),
                X_ice_CH4_m=zeros(N),
                X_ice_N2_m=zeros(N),
                X_ice_H2S_m=zeros(N),
                X_ice_PH3_m=zeros(N),
                X_refr_C_m=zeros(N),
                X_refr_S_m=zeros(N),
                X_refr_N_m=zeros(N),
                X_refr_P_m=zeros(N),
                X_refr_H_m=zeros(N),
            ),
            Xmin_troilite_m=zeros(N),
            Xmin_schreibersite_m=zeros(N),
            Xmin_cohenite_m=zeros(N),
            Xmin_graphite_m=zeros(N),
            Xmin_nitride_m=zeros(N),
            Xmin_metal_matrix_m=zeros(N),
            t_accreted=zeros(N),
        )
        rec_full = Erebus._reconstruct_checkpoint_marker_arrays(core, opt_full, N)
        @test haskey(rec_full.groups, :metal)
        @test haskey(rec_full.groups, :volatiles)
        @test haskey(rec_full.groups, :redox)
        @test haskey(rec_full.groups, :hcnspo)
        @test haskey(rec_full.groups, :phase)
        @test haskey(rec_full.groups, :accretion)

        # Partial (zeros fallback for Xfe_H_m etc.)
        opt_partial = (;
            Xfem=zeros(N),
            Xfem0=zeros(N),
            Xfe_bulk=fill(0.2, N),
            Xfe_H_m=nothing,
            Xfe_C_m=nothing,
            Xfe_N_m=nothing,
            Xfe_S_m=nothing,
            XH2Om=nothing,
            XCm=nothing,
            XNm=nothing,
            XSm=nothing,
            X_graphite_m=nothing,
            F_extract_m=zeros(N),
            redox_props=nothing,
            hcnspo_props=nothing,
            Xmin_troilite_m=nothing,
            Xmin_schreibersite_m=nothing,
            Xmin_cohenite_m=nothing,
            Xmin_graphite_m=nothing,
            Xmin_nitride_m=nothing,
            Xmin_metal_matrix_m=nothing,
            t_accreted=nothing,
        )
        rec_part = Erebus._reconstruct_checkpoint_marker_arrays(core, opt_partial, N)
        @test haskey(rec_part.groups, :metal)
        @test haskey(rec_part.groups, :volatiles)
        @test !haskey(rec_part.groups, :redox)
        @test length(rec_part.groups.metal.Xfe_H_m) == N

        # Minimal
        opt_min = (;
            Xfem=nothing,
            Xfem0=nothing,
            Xfe_bulk=nothing,
            Xfe_H_m=nothing,
            Xfe_C_m=nothing,
            Xfe_N_m=nothing,
            Xfe_S_m=nothing,
            XH2Om=nothing,
            XCm=nothing,
            XNm=nothing,
            XSm=nothing,
            X_graphite_m=nothing,
            F_extract_m=nothing,
            redox_props=nothing,
            hcnspo_props=nothing,
            Xmin_troilite_m=nothing,
            Xmin_schreibersite_m=nothing,
            Xmin_cohenite_m=nothing,
            Xmin_graphite_m=nothing,
            Xmin_nitride_m=nothing,
            Xmin_metal_matrix_m=nothing,
            t_accreted=nothing,
        )
        rec_min = Erebus._reconstruct_checkpoint_marker_arrays(core, opt_min, N)
        @test isempty(rec_min.groups)
    end
end
