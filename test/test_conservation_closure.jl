# test/test_conservation_closure.jl
# 20-step integration test on magma_ocean_cooling_turb_on_32.toml verifying
# non-triviality, 2D planar conservation, 3D reservoir balance, cross-frame
# relations, and oxygen accounting in redox-off and redox-on regimes.

using Test
using Random
using Erebus
using Erebus.Config
using Erebus.Geometry
using Erebus.Numerics
using Erebus.Particles
using Erebus.Physics
using Erebus.Simulation

@testset "Conservation Closure Integration" begin
    @testset "Redox Off Base Integration" begin
        cfg_base = load_config(
            normpath(
                joinpath(@__DIR__, "..", "configs", "magma_ocean_cooling_turb_on_32.toml")
            ),
        )
        cfg = Erebus.override_config(
            cfg_base,
            Dict(
                "grid.Nx" => 33,
                "grid.Ny" => 33,
                "time.n_steps" => 20,
                "time.dt_initial" => 0.05,
                "time.dt_longest" => 0.05,
                "solver.p2m_mode" => :tiled,
                "solver.hydromech_solver" => :direct,
                "solver.seed" => 42,
                "melting.active" => true,
                "volatiles.active" => true,
                "volatiles.initial_water_wtpct" => 1.0,
                "volatiles.initial_carbon_ppm" => 500.0,
                "volatiles.initial_nitrogen_ppm" => 50.0,
                "volatiles.initial_sulfur_ppm" => 1000.0,
                "venting.active" => true,
                "venting.species" => :H2,
                "retention.active" => true,
                "retention.venting_drainage_active" => true,
                "magma_degassing.active" => true,
                "magma_degassing.mode" => :dynamic_flux,
                "magma_degassing.degas_depth_fraction" => 0.90,
                "atmosphere.active" => true,
                "escape.active" => true,
                "escape.multi_species" => true,
                "metal_partition.active" => true,
                "reaction.active" => false,
                "accretion.active" => false,
                "telescoping.active" => false,
                "magma_transport.active" => false,
                "coreformation.percolation_active" => false,
                "coreformation.settling_active" => false,
                "hydrothermal.active" => false,
                "phase_tracking.active" => false,
                "disk.enabled" => false,
                "volatile_mixture.active" => false,
                "refractory.active" => false,
                "mpi.enable" => false,
                "redox.active" => false,
                "output.mode" => :telemetry,
                "output.save_final" => false,
                "output.output_dir" => tempname(),
            ),
        )

        @test validate_config(cfg) === nothing
        res = simulation_loop(cfg)
        @test length(res.transfers) >= 100

        # 1. Non-triviality assertions
        degas_records = [
            r for r in res.transfers if r.channel === :degassing && r.dM2 > 0.0
        ]
        vent_records = [
            r for r in res.transfers if
            (r.channel === :venting || r.channel === :pore_venting) && r.dM2 > 0.0
        ]
        @test length(degas_records) >= 10
        @test length(vent_records) >= 10
        @test count(>(0.0), (res.atm.escaped.H,)) >= 1
        @test count(>(0.40), res.markers.Fm) >= 1

        # 2. Cross-frame geometric assertion (dM3 == dM2 * 2 * rm)
        xc = cfg.geometry.xcenter
        yc = cfg.geometry.ycenter
        max_cross_rel_err = 0.0
        for r in res.transfers
            rm = hypot(r.x - xc, r.y - yc)
            expected_dM3 = r.dM2 * (2.0 * rm)
            if abs(expected_dM3) > 0.0
                err = abs(r.dM3 - expected_dM3) / abs(expected_dM3)
                max_cross_rel_err = max(max_cross_rel_err, err)
            end
        end
        @test max_cross_rel_err < 1.0e-12
        @test isapprox(max_cross_rel_err, 0.0; atol=1.0e-12)

        # 3. 3D reservoir balance: atm.elem + atm.escaped == sum(dM3)
        elem_dM3 = Dict(:H => 0.0, :C => 0.0, :N => 0.0, :S => 0.0)
        elem_dM2 = Dict(:H => 0.0, :C => 0.0, :N => 0.0, :S => 0.0)
        for r in res.transfers
            if haskey(elem_dM3, r.element)
                elem_dM3[r.element] += r.dM3
                elem_dM2[r.element] += r.dM2
            end
        end

        for E in [:H, :C, :N, :S]
            m_atm = getproperty(res.atm.elem, E)
            m_esc = getproperty(res.atm.escaped, E)
            tot_3d = m_atm + m_esc
            @test isapprox(tot_3d, elem_dM3[E]; rtol=1.0e-8, atol=1.0e-3)
        end

        # 4. 2D planar volatile conservation
        coords = GridCoordinates(
            cfg.grid.Nx, cfg.grid.Ny; xsize=cfg.grid.xsize, ysize=cfg.grid.ysize
        )
        Am = marker_area(coords)
        h_conv = 2.01588 / 18.01528
        rho_sil = cfg.materials.rhosolidm[1]
        rho_met = cfg.coreformation.rho_metal
        rho_fluid = cfg.materials.rhofluidm[2]
        tm = res.markers.tm
        marknum = length(tm)
        Xfe_bulk = res.markers.Xfe_bulk

        # Final 2D inventories
        H_fin = 0.0
        C_fin = 0.0
        N_fin = 0.0
        S_fin = 0.0
        for m in 1:marknum
            if tm[m] < 3
                phi_fe = Xfe_bulk !== nothing ? clamp(Float64(Xfe_bulk[m]), 0.0, 1.0) : 0.0
                phi_sil = max(0.0, 1.0 - phi_fe)
                m_sil_2d = phi_sil * rho_sil * Am
                m_met_2d = phi_fe * rho_met * Am
                m_pore_2d = res.markers.phim[m] * rho_fluid * Am

                H_fin +=
                    res.markers.XH2Om[m] * 0.01 * h_conv * m_sil_2d + m_pore_2d * h_conv
                C_fin += res.markers.XCm[m] * 1.0e-6 * m_sil_2d
                N_fin += res.markers.XNm[m] * 1.0e-6 * m_sil_2d
                S_fin += res.markers.XSm[m] * 1.0e-6 * m_sil_2d

                if haskey(res.markers, :Xfe_H_m)
                    H_fin += res.markers.Xfe_H_m[m] * 1.0e-6 * m_met_2d
                    C_fin += res.markers.Xfe_C_m[m] * 1.0e-6 * m_met_2d
                    N_fin += res.markers.Xfe_N_m[m] * 1.0e-6 * m_met_2d
                    S_fin += res.markers.Xfe_S_m[m] * 1.0e-6 * m_met_2d
                end
            end
        end

        # Initial 2D inventories recomputed from seed 42
        Random.seed!(Erebus.rgen, cfg.solver.seed)
        (xm_i, ym_i, tm_i, tkm_i, sxxm_i, sxym_i, etavpm_i, phim_i, phinewm_i, pfm0_i, XWsolidm_i, XWsolidm0_i, Fm_i) = setup_marker_properties(
            marknum, coords
        )
        (Xfem_i, Xfem0_i, Xfe_bulk_i) = setup_marker_metal_properties(marknum)
        (rhotot_i, rhocptot_i, etatot_i, hrtot_i, ktot_i, tkm_rhocptot_i, etafluid_inv_k_i, inv_ggg_i, frict_i, cohes_i, tens_i, rhofluid_i, alphasolid_i, alphafluid_i) = setup_marker_properties_helpers(
            marknum
        )
        define_markers!(
            xm_i,
            ym_i,
            tm_i,
            phim_i,
            etavpm_i,
            rhotot_i,
            rhocptot_i,
            etatot_i,
            hrtot_i,
            ktot_i,
            tkm_i,
            inv_ggg_i,
            frict_i,
            cohes_i,
            tens_i,
            rhofluid_i,
            alphasolid_i,
            alphafluid_i,
            XWsolidm0_i;
            coords=coords,
            xcenter_val=cfg.geometry.xcenter,
            ycenter_val=cfg.geometry.ycenter,
            rplanet_val=cfg.geometry.rplanet,
            rcrust_val=cfg.geometry.rcrust,
            XWsolidm_init_val=cfg.materials.XWsolidm_init,
            phim0_val=cfg.thermodynamics.phim0,
            Xfe_bulk=Xfe_bulk_i,
            Xfem=Xfem_i,
            Xfem0=Xfem0_i,
            Xfe_bulk_val=cfg.coreformation.Xfe_bulk,
            tkm0_val=cfg.materials.tkm0,
        )

        H_init = 0.0
        C_init = 0.0
        N_init = 0.0
        S_init = 0.0
        for m in 1:marknum
            if tm_i[m] < 3
                phi_fe = clamp(Float64(Xfe_bulk_i[m]), 0.0, 1.0)
                phi_sil = max(0.0, 1.0 - phi_fe)
                m_sil_2d = phi_sil * rho_sil * Am
                m_met_2d = phi_fe * rho_met * Am
                m_pore_2d_init = phim_i[m] * rho_fluid * Am

                H_init +=
                    (cfg.volatiles.initial_water_wtpct * 0.01 * h_conv) * m_sil_2d +
                    (cfg.metal_partition.initial_metal_h_ppm * 1.0e-6) * m_met_2d +
                    m_pore_2d_init * h_conv
                C_init +=
                    (cfg.volatiles.initial_carbon_ppm * 1.0e-6) * m_sil_2d +
                    (cfg.metal_partition.initial_metal_c_ppm * 1.0e-6) * m_met_2d
                N_init +=
                    (cfg.volatiles.initial_nitrogen_ppm * 1.0e-6) * m_sil_2d +
                    (cfg.metal_partition.initial_metal_n_ppm * 1.0e-6) * m_met_2d
                S_init +=
                    (cfg.volatiles.initial_sulfur_ppm * 1.0e-6) * m_sil_2d +
                    (cfg.metal_partition.initial_metal_s_ppm * 1.0e-6) * m_met_2d
            end
        end

        # Unvented compaction loss of pore fluid in interior Darcy flow
        H_compaction_unvented = 0.0
        for m in 1:marknum
            if tm[m] < 3
                dphi = phim_i[m] - res.markers.phim[m]
                H_compaction_unvented += dphi * rho_fluid * Am * h_conv
            end
        end
        # Subtract vented pore fluid component already counted in elem_dM2[:H]
        H_vented_pore = 0.0
        for r in res.transfers
            if r.channel === :pore_venting && r.element === :H
                H_vented_pore += r.dM2
            end
        end
        H_compaction_unvented = max(0.0, H_compaction_unvented - H_vented_pore)

        @test isapprox(C_fin + elem_dM2[:C], C_init; rtol=1.0e-8, atol=1.0e-6)
        @test isapprox(N_fin + elem_dM2[:N], N_init; rtol=1.0e-8, atol=1.0e-6)
        @test isapprox(S_fin + elem_dM2[:S], S_init; rtol=1.0e-8, atol=1.0e-6)
        @test isapprox(
            H_fin + elem_dM2[:H] + H_compaction_unvented, H_init; rtol=1.0e-8, atol=1.0e-6
        )

        # 5. Oxygen accounting in redox-off regime
        # Delivered water oxygen = degas H2O (delivered as H2O) + mineral H2O (delivered as H2O) + pore H2O (converted to H2)
        # In redox-off mode:
        # - Vented pore water converted to H2 leaves O in the rock buffer: dO_buffer = pore_H_3D * (15.9994 / 2.01588)
        # - Vented mineral water enters atmosphere as H2O (both H and O enter atmosphere)
        # - Degassed water enters atmosphere as H2O (both H and O enter atmosphere)
        total_atm_esc_O = res.atm.elem.O + res.atm.escaped.O
        total_transfer_O = sum(
            r.dM3 * (15.9994 / 2.01588) for r in res.transfers if r.element === :H
        )
        dO_buffer = sum(
            r.dM3 * (15.9994 / 2.01588) for
            r in res.transfers if r.channel === :pore_venting && r.element === :H
        )
        @test isapprox(
            total_atm_esc_O + dO_buffer, total_transfer_O; rtol=1.0e-8, atol=1.0e-3
        )
        @test count(>(1.0e10), (total_atm_esc_O,)) >= 1
    end

    @testset "Redox On Integration" begin
        cfg_base = load_config(
            normpath(
                joinpath(@__DIR__, "..", "configs", "magma_ocean_cooling_turb_on_32.toml")
            ),
        )
        cfg = Erebus.override_config(
            cfg_base,
            Dict(
                "grid.Nx" => 33,
                "grid.Ny" => 33,
                "time.n_steps" => 20,
                "time.dt_initial" => 0.05,
                "time.dt_longest" => 0.05,
                "solver.p2m_mode" => :tiled,
                "solver.hydromech_solver" => :direct,
                "solver.seed" => 42,
                "melting.active" => true,
                "volatiles.active" => true,
                "volatiles.initial_water_wtpct" => 1.0,
                "volatiles.initial_carbon_ppm" => 500.0,
                "volatiles.initial_nitrogen_ppm" => 50.0,
                "volatiles.initial_sulfur_ppm" => 1000.0,
                "venting.active" => true,
                "venting.species" => :H2,
                "retention.active" => true,
                "retention.venting_drainage_active" => true,
                "magma_degassing.active" => true,
                "magma_degassing.mode" => :dynamic_flux,
                "magma_degassing.degas_depth_fraction" => 0.90,
                "atmosphere.active" => true,
                "escape.active" => true,
                "escape.multi_species" => true,
                "metal_partition.active" => true,
                "reaction.active" => false,
                "accretion.active" => false,
                "telescoping.active" => false,
                "magma_transport.active" => false,
                "coreformation.percolation_active" => false,
                "coreformation.settling_active" => false,
                "hydrothermal.active" => false,
                "phase_tracking.active" => false,
                "disk.enabled" => false,
                "volatile_mixture.active" => false,
                "refractory.active" => false,
                "mpi.enable" => false,
                "redox.active" => true,
                "output.mode" => :telemetry,
                "output.save_final" => false,
                "output.output_dir" => tempname(),
            ),
        )

        @test validate_config(cfg) === nothing
        res = simulation_loop(cfg)
        @test length(res.transfers) >= 100

        # Oxygen accounting in redox-on regime
        total_atm_esc_O = res.atm.elem.O + res.atm.escaped.O
        total_transfer_O = sum(
            r.dM3 * (15.9994 / 2.01588) for r in res.transfers if r.element === :H
        )
        dO_buffer = sum(
            r.dM3 * (15.9994 / 2.01588) for
            r in res.transfers if r.channel === :pore_venting && r.element === :H
        )
        @test isapprox(
            total_atm_esc_O + dO_buffer, total_transfer_O; rtol=1.0e-8, atol=1.0e-3
        )
        @test count(>(1.0e10), (total_atm_esc_O,)) >= 1
    end
end
