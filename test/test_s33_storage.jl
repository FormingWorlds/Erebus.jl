using Test
using JLD2
using Erebus

# Small simulation state; `extra` overrides config fields.
function _s33_state(extra=Dict{String,Any}())
    cfg = Erebus.override_config(
        default_config(), merge(Dict{String,Any}("grid.Nx" => 17, "grid.Ny" => 17), extra)
    )
    state, coords, ws, _, _ = init_simulation(cfg)
    return cfg, state, coords, ws
end

@testset "Out-of-plane stress S33 storage" begin
    @testset "GridArrays and CoreGroup carry the S33 fields" begin
        cfg, state, coords, ws = _s33_state()
        @test all(fn -> hasfield(Erebus.GridArrays, fn), (:S33, :S330, :DS33))
        @test size(state.grids.S33) == size(state.grids.SXX)
        @test size(state.grids.S330) == size(state.grids.SXX0)
        @test size(state.grids.DS33) == size(state.grids.DSXX)
        @test fieldnames(Erebus.CoreGroup)[findfirst(==(:sxym), fieldnames(Erebus.CoreGroup)) + 1] ===
            :s33m
        @test all(iszero, state.markers.core.s33m)
        @test length(state.markers.core.s33m) == length(state.markers)
    end

    @testset "marker to grid interpolation gives S330 = c" begin
        c = 3.7e5
        for mode in (:buffered, :tiled)
            cfg, state, coords, ws = _s33_state(Dict{String,Any}("solver.p2m_mode" => mode))
            state.markers.core.s33m .= c
            state.grids.S330 .= -1.0
            interpolate_markers_to_grid!(
                state,
                coords,
                cfg;
                p2m_workspace=ws.p2m,
                thread_buffers=ws.thread_buffers,
                interp_arrays=ws.interp_arrays,
            )
            WTPSUM = ws.interp_arrays[33]
            filled = WTPSUM .> 0.0
            @test any(filled)
            @test all(isapprox.(state.grids.S330[filled], c; rtol=1.0e-12))
            @test all(iszero, state.grids.S330[.!filled])
        end
    end

    @testset "thread buffers reset and reduce the S33 sum" begin
        coords = GridCoordinates(default_config().grid)
        bufs = Erebus.allocate_thread_interpolation_buffers(2, coords)
        for b in bufs
            b.S33SUM .= 2.5
        end
        master = zeros(size(bufs[1].S33SUM))
        reduce_dummy = Erebus.setup_interpolated_properties(coords)
        Erebus.reduce_thread_buffers!(reduce_dummy[1:30]..., bufs; S33SUM=master)
        @test all(isapprox.(master, 5.0; rtol=1.0e-12))
        @test all(iszero, reduce_dummy[34])
        Erebus.reset_thread_buffers!(bufs)
        @test all(all(iszero, b.S33SUM) for b in bufs)
        @test length(reduce_dummy) == 34
    end

    @testset "DS33 updates the marker stress" begin
        c, d = 2.0e5, -7.5e3
        cfg, state, coords, ws = _s33_state()
        core = state.markers.core
        core.s33m .= c
        state.grids.DS33 .= d
        Erebus.update_marker_stress!(
            core.xm,
            core.ym,
            core.sxxm,
            core.sxym,
            state.grids.DSXX,
            state.grids.DSXY,
            length(state.markers);
            coords=coords,
            s33m=core.s33m,
            DS33=state.grids.DS33,
        )
        @test all(isapprox.(core.s33m, c + d; rtol=1.0e-12))
        @test all(iszero, core.sxxm)
    end

    @testset "subgrid relaxation moves s33m toward S330 and corrects DS33" begin
        c = 1.0e6
        cfg, state, coords, ws = _s33_state()
        core = state.markers.core
        g = state.grids
        core.s33m .= c
        g.S330 .= 0.0
        g.DS33 .= 0.0
        SXXSUM, SXYSUM = ws.interp_arrays[27], ws.interp_arrays[4]
        WTPSUM, WTSUM, S33SUM = ws.interp_arrays[33],
        ws.interp_arrays[8],
        ws.interp_arrays[34]
        Erebus.apply_subgrid_stress_diffusion!(
            core.xm,
            core.ym,
            core.tm,
            core.inv_gggtotalm,
            core.sxxm,
            core.sxym,
            g.SXX0,
            g.SXY0,
            g.DSXX,
            g.DSXY,
            SXXSUM,
            SXYSUM,
            WTPSUM,
            WTSUM,
            1.0e10,
            length(state.markers);
            coords=coords,
            dsubgrids=1.0,
            etavpm=core.etavpm,
            etasolidm=cfg.materials.etasolidm,
            s33m=core.s33m,
            S330=g.S330,
            DS33=g.DS33,
            S33SUM=S33SUM,
        )
        @test all(0.0 - 1.0e-9 .<= core.s33m .<= c)
        @test any(<(c), core.s33m)
        # The grid increment is the opposite of the marker relaxation, so their sum is conserved.
        @test all(0.0 .<= g.DS33 .<= c * (1 + 1.0e-12))
        @test !iszero(maximum(g.DS33))
    end

    @testset "replenish_markers! copies s33m of the source marker" begin
        cfg, state, coords, ws = _s33_state()
        markers = state.markers
        core = markers.core
        n0 = length(markers)
        core.xm .*= 0.5
        for m in 1:n0
            core.s33m[m] = 1.0e3 + m
            core.pfm0[m] = 2.0 * core.s33m[m]
        end
        mdis, mnum = Erebus.setup_marker_geometry_helpers(coords)
        n1 = Erebus.replenish_markers!(markers, mdis, mnum; randomized=false, coords=coords)
        @test n1 > n0
        @test length(core.s33m) == n1
        new = (n0 + 1):n1
        @test all(isapprox.(core.pfm0[new], 2.0 .* core.s33m[new]; rtol=1.0e-12))
        @test all(core.s33m[new] .> 1.0e3)
    end

    @testset "telescoping and accretion give zero s33m for new markers" begin
        cfg, state, coords, ws = _s33_state()
        markers = state.markers
        n0 = length(markers)
        markers.core.s33m .= 4.0e4
        new_coords = Erebus.compute_telescoped_coordinates(coords)
        n1 = Erebus.telescope_marker_arrays!(
            markers; old_coords=coords, new_coords=new_coords, cfg=cfg
        )
        @test n1 > n0
        @test length(markers.core.s33m) == n1
        @test all(isapprox.(markers.core.s33m[1:n0], 4.0e4; rtol=1.0e-12))
        @test all(iszero, markers.core.s33m[(n0 + 1):n1])

        core = Erebus.init_marker_arrays(3, cfg, coords).core
        core.s33m .= 9.0
        params = (phi_accreted=0.1, T_accreted=200.0, XW_acc=0.0)
        Erebus._init_accreted_marker_core!(core, 2, params, cfg)
        @test iszero(core.s33m[2])
        @test all(isapprox.(core.s33m[[1, 3]], 9.0; rtol=1.0e-12))
    end

    @testset "telescope_domain! remaps S33, S330 and DS33 like SXX" begin
        cfg, state, coords, ws = _s33_state(Dict{String,Any}("telescoping.active" => true))
        g = state.grids
        g.S33 .= 1.0
        g.S330 .= 2.0
        g.DS33 .= 3.0
        state.accumulators.rplanet = 0.99 * coords.xsize / 2
        coords_ref = Ref(coords)
        @test Erebus.telescope_domain!(state, coords_ref, cfg, ws)
        g2 = state.grids
        @test size(g2.S33) == size(g2.SXX)
        @test size(g2.S330) == size(g2.SXX0)
        @test size(g2.DS33) == size(g2.DSXX)
        @test all(isapprox.(extrema(g2.S33), (0.0, 1.0); atol=1.0e-12))
        @test maximum(g2.S330) ≈ 2.0 atol = 1.0e-12
        @test maximum(g2.DS33) ≈ 3.0 atol = 1.0e-12
    end

    @testset "checkpoint keeps S33, S330, DS33 and s33m" begin
        mktempdir() do dir
            cfg, state, coords, ws = _s33_state()
            g = state.grids
            g.S33 .= 1.5e3
            g.S330 .= -2.5e3
            g.DS33 .= 0.75e3
            state.markers.core.s33m .= 6.25e4
            path = Erebus.save_state(dir, state, coords, cfg)
            loaded, _, _ = Erebus.load_simulation_state(path)
            @test loaded.grids.S33 == g.S33
            @test loaded.grids.S330 == g.S330
            @test loaded.grids.DS33 == g.DS33
            @test loaded.markers.core.s33m == state.markers.core.s33m
            @test all(isapprox.(loaded.markers.core.s33m, 6.25e4; rtol=1.0e-12))
        end
    end

    @testset "checkpoint without S33 fields: zeros, or an error with spherical_flow" begin
        mktempdir() do dir
            for (flow, extra) in (
                (false, Dict{String,Any}()),
                (true, Dict{String,Any}("geometry.spherical_flow" => true)),
            )
                cfg, state, coords, ws = _s33_state(extra)
                path = Erebus.save_state(
                    joinpath(dir, "flow_$(flow).jld2"), state, coords, cfg
                )
                data = JLD2.load(path)
                for k in ("S33", "S330", "DS33", "s33m")
                    delete!(data, k)
                end
                older = joinpath(dir, "older_$(flow).jld2")
                JLD2.save(older, data)
                if flow
                    err = try
                        Erebus.load_simulation_state(older)
                        nothing
                    catch e
                        e
                    end
                    @test err isa Erebus.CheckpointError
                    @test occursin("has no S33 stress history", err.msg)
                else
                    loaded, _, _ = Erebus.load_simulation_state(older)
                    @test all(iszero, loaded.grids.S33)
                    @test all(iszero, loaded.markers.core.s33m)
                    @test size(loaded.grids.DS33) == size(loaded.grids.DSXX)
                end
            end
        end
    end
end
