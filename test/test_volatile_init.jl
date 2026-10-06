using Test
using Erebus
using Random
using StaticArrays

make_test_coordinates(Nx=33, Ny=33) = default_grid_coordinates()

@testset "Volatile Budgets and Initial Inventories (PR 3: tl/volatile-init)" begin
    @testset "Finding F05: Accretion Water & Ice Separation" begin
        coords = make_test_coordinates(33, 33)
        cfg = Erebus.override_config(
            default_config(),
            Dict(
                "accretion.active" => true,
                "accretion.snowline_coupling" => true,
                "accretion.XH2O_wet_wtpct" => 10.0,
                "accretion.XH2O_dry_wtpct" => 0.1,
                "accretion.XWsolid_wet" => 0.40,
                "accretion.XWsolid_dry" => 0.0,
                "accretion.Xfe_bulk_accreted" => 0.15,
                "volatile_mixture.active" => true,
                "volatile_mixture.X_ice_H2O" => 0.85,
                "volatiles.active" => true,
                "refractory.active" => true,
                "redox.active" => true,
            ),
        )

        markers = init_marker_arrays(10, cfg, coords)
        # Set marker 1 and 2 as air markers inside the accretion shell
        markers.tm[1] = 3
        markers.xm[1] = coords.xcenter + 500.0
        markers.ym[1] = coords.ycenter + 500.0
        markers.tm[2] = 3
        markers.xm[2] = coords.xcenter + 600.0
        markers.ym[2] = coords.ycenter + 600.0

        disk_state_wet = (;
            condensed_H2O=true,
            condensed_NH3=true,
            condensed_CO2=true,
            condensed_CO=false,
            condensed_CH4=false,
            condensed_N2=false,
            condensed_H2S=true,
            condensed_PH3=false,
            X_ice_H2O=0.85,
            X_ice_NH3=0.03,
            X_ice_CO2=0.05,
            X_ice_CO=0.0,
            X_ice_CH4=0.0,
            X_ice_N2=0.0,
            X_ice_H2S=0.01,
            X_ice_PH3=0.0,
            f_refr_C=0.02,
            f_refr_S=0.015,
            f_refr_N=0.001,
            f_refr_P=0.0005,
            f_refr_H=0.0002,
        )

        # Accrete shell below snowline
        n_conv = advance_accretion_boundary!(
            0.0, 1000.0, markers; cfg=cfg, disk_state=disk_state_wet, current_time=1.0e5
        )
        @test n_conv == 2
        @test markers.tm[1] == 2
        @test markers.tm[2] == 2

        # Rock mineral water must follow accretion wet water settings, NOT cryogenic ice fraction * 100
        @test markers.groups.volatiles.XH2Om[1] ≈ 10.0
        @test markers.groups.volatiles.XH2Om[1] != 85.0
        @test markers.core.XWsolidm[1] ≈ 0.40
        @test markers.core.XWsolidm[1] != 0.85
        @test markers.core.XWsolidm0[1] ≈ 0.40

        # Dedicated HCN-S-P-O ice arrays receive cryogenic ice mass fraction
        @test markers.groups.hcnspo.X_ice_H2O_m[1] ≈ 0.85
        @test markers.groups.hcnspo.X_ice_CO2_m[1] ≈ 0.05
        @test markers.groups.hcnspo.X_refr_C_m[1] ≈ 0.02

        # Dry accretion shell test
        markers.tm[3] = 3
        markers.xm[3] = coords.xcenter + 700.0
        markers.ym[3] = coords.ycenter + 700.0

        disk_state_dry = (;
            condensed_H2O=false,
            condensed_NH3=false,
            condensed_CO2=false,
            condensed_CO=false,
            condensed_CH4=false,
            condensed_N2=false,
            condensed_H2S=false,
            condensed_PH3=false,
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
            0.0, 1000.0, markers; cfg=cfg, disk_state=disk_state_dry, current_time=2.0e5
        )
        @test n_conv2 == 1
        @test markers.tm[3] == 2
        @test markers.groups.volatiles.XH2Om[3] ≈ 0.1
        @test markers.core.XWsolidm[3] ≈ 0.0
        @test markers.groups.hcnspo.X_ice_H2O_m[3] ≈ 0.0

        # Direct advance_accretion_boundary_hcnspo! method verification
        marknum_h = 4
        xm_h = [
            coords.xcenter + 200.0,
            coords.xcenter + 400.0,
            coords.xcenter + 800.0,
            coords.xcenter + 1200.0,
        ]
        ym_h = fill(coords.ycenter, marknum_h)
        tm_h = fill(3, marknum_h)
        tkm_h = fill(150.0, marknum_h)
        phim_h = fill(0.35, marknum_h)
        XWsolidm0_h = fill(0.0, marknum_h)
        hcn_h = HcnspoGroup([zeros(marknum_h) for _ in 1:13]...)

        n_hcn = advance_accretion_boundary_hcnspo!(
            0.0,
            500.0,
            xm_h,
            ym_h,
            tm_h,
            tkm_h,
            phim_h,
            XWsolidm0_h,
            hcn_h,
            disk_state_wet;
            xcenter=coords.xcenter,
            ycenter=coords.ycenter,
            XWsolid_accreted=0.40,
        )
        @test n_hcn == 2
        @test tm_h[1] == 2
        @test tm_h[2] == 2
        @test XWsolidm0_h[1] ≈ 0.40
        @test XWsolidm0_h[1] != 0.85
        @test hcn_h.X_ice_H2O_m[1] ≈ 0.85
    end

    @testset "Finding F27: Initial HCNSPO Inventories on All Planet Markers" begin
        coords = make_test_coordinates(33, 33)
        cfg = Erebus.override_config(
            default_config(),
            Dict(
                "geometry.rplanet" => 50000.0,
                "geometry.rcrust" => 50000.0,
                "volatile_mixture.active" => true,
                "volatile_mixture.X_ice_H2O" => 0.85,
                "volatile_mixture.X_ice_CO2" => 0.05,
                "refractory.active" => true,
                "refractory.f_refr_C" => 0.03,
                "refractory.f_refr_S" => 0.02,
            ),
        )

        marknum = coords.start_marknum
        markers = init_marker_arrays(marknum, cfg, coords)
        Erebus.define_markers!(
            markers; cfg=cfg, coords=coords, rplanet_val=50000.0, rcrust_val=50000.0
        )

        n_planet = 0
        n_air = 0
        for m in 1:marknum
            if markers.tm[m] < 3
                n_planet += 1
                # All planet rock markers must have initial volatile mixture and refractory inventories
                @test markers.groups.hcnspo.X_ice_H2O_m[m] ≈ 0.85
                @test markers.groups.hcnspo.X_ice_CO2_m[m] ≈ 0.05
                @test markers.groups.hcnspo.X_refr_C_m[m] ≈ 0.03
                @test markers.groups.hcnspo.X_refr_S_m[m] ≈ 0.02
            else
                n_air += 1
                # Sticky air markers must remain vacuum (iszero)
                @test iszero(markers.groups.hcnspo.X_ice_H2O_m[m])
                @test iszero(markers.groups.hcnspo.X_ice_CO2_m[m])
                @test iszero(markers.groups.hcnspo.X_refr_C_m[m])
                @test iszero(markers.groups.hcnspo.X_refr_S_m[m])
            end
        end

        @test n_planet + n_air == marknum
        @test n_planet == count(m -> markers.tm[m] < 3, 1:marknum)
        @test n_air == count(m -> markers.tm[m] == 3, 1:marknum)
    end

    @testset "Element Conservation and 3D Budget Closure" begin
        coords = make_test_coordinates(33, 33)
        r_planet = 40000.0
        cfg = Erebus.override_config(
            default_config(),
            Dict(
                "geometry.rplanet" => r_planet,
                "geometry.rcrust" => r_planet,
                "volatile_mixture.active" => true,
                "volatile_mixture.X_ice_H2O" => 0.10,
                "refractory.active" => true,
                "refractory.f_refr_C" => 0.02,
                "thermodynamics.phim0" => 0.05,
            ),
        )

        rho_rock = cfg.materials.rhosolidm[1]
        marknum = coords.start_marknum
        markers = init_marker_arrays(marknum, cfg, coords)
        Erebus.define_markers!(
            markers; cfg=cfg, coords=coords, rplanet_val=r_planet, rcrust_val=r_planet
        )

        # 3D spherical analytical volume and solid rock mass with porosity
        phim0_val = Erebus.phim0
        v_sphere = (4.0 / 3.0) * pi * r_planet^3
        m_sphere_expected = v_sphere * rho_rock * (1.0 - phim0_val)
        m_carbon_expected = m_sphere_expected * 0.02

        # 2D Cartesian cylindrical proxy to 3D sphere: weight is 2 * r_m
        # Marker area: coords.dxm * coords.dym
        m_area = Erebus.marker_area(coords)

        m_sphere_integrated = 0.0
        m_carbon_integrated = 0.0
        for m in 1:marknum
            if markers.tm[m] < 3
                dx_m = markers.xm[m] - coords.xcenter
                dy_m = markers.ym[m] - coords.ycenter
                r_m = sqrt(dx_m^2 + dy_m^2)
                w3d = markers.w3d_m[m]
                dm_3d = rho_rock * (1.0 - markers.phim[m]) * m_area * w3d
                m_sphere_integrated += dm_3d
                m_carbon_integrated += dm_3d * markers.groups.hcnspo.X_refr_C_m[m]
            end
        end

        # Discretization tolerance for 1000 markers on Cartesian grid
        @test isapprox(m_sphere_integrated, m_sphere_expected, rtol=0.10)
        @test isapprox(m_carbon_integrated, m_carbon_expected, rtol=0.10)
    end
end
