using Test
using Random
using Erebus
using Erebus.Numerics: assemble_thermal_lse!, compute_fluid_advection_heat!

"""Build a P grid temperature field that is linear in x and y."""
function _linear_temperature(coords, T0, gx, gy)
    return [T0 + gx * coords.xp[j] + gy * coords.yp[i] for i in 1:coords.Ny1, j in 1:coords.Nx1]
end

"""Max error of the steady advection-conduction profile, march in time to steady state."""
function _steady_profile_error(n, Pe)
    coords = Erebus.GridCoordinates(n, 3; xsize=1.0, ysize=3.0 / n)
    h = coords.dx
    L = coords.xp[n + 1] - coords.xp[2]
    u = 1.0
    kappa = u * L / Pe
    x = (coords.xp .- coords.xp[2]) ./ L
    T = zeros(coords.Ny1, coords.Nx1)
    T[:, (n + 1):end] .= 1.0
    qx = fill(u, coords.Ny1, coords.Nx1)
    qy = zeros(coords.Ny1, coords.Nx1)
    Q = zeros(coords.Ny1, coords.Nx1)
    dt = 0.4 * min(h / u, h^2 / (2.0 * kappa))
    for _ in 1:ceil(Int, 20.0 * L / (u * dt))
        compute_fluid_advection_heat!(Q, T, qx, qy, _ -> 1.0; coords=coords)
        Tnew = copy(T)
        for j in 3:n, i in axes(T, 1)
            i_q = clamp(i, 2, coords.Ny)
            Tnew[i, j] = T[i, j] + dt * (Q[i_q, j] + kappa * (T[i, j - 1] - 2.0 * T[i, j] + T[i, j + 1]) / h^2)
        end
        T = Tnew
    end
    exact = @. (exp(Pe * x[2:(n + 1)]) - 1.0) / (exp(Pe) - 1.0)
    return maximum(abs, T[2, 2:(n + 1)] .- exact)
end

@testset "Fluid advection heat on the P grid" begin
    coords = Erebus.GridCoordinates(20, 16; xsize=1.0e4, ysize=8.0e3)
    f = 2.5e6
    Q = zeros(coords.Ny1, coords.Nx1)

    @testset "uniform flow and gradient: exact value on all 4 inflow faces" begin
        for (qxv, qyv, gx, gy) in (
            (3.0e-8, 0.0, 0.05, 0.0),
            (-3.0e-8, 0.0, 0.05, 0.0),
            (0.0, 3.0e-8, 0.0, 0.05),
            (0.0, -3.0e-8, 0.0, 0.05),
            (2.0e-8, -1.0e-8, 0.03, 0.07),
        )
            T = _linear_temperature(coords, 500.0, gx, gy)
            qx = fill(qxv, coords.Ny1, coords.Nx1)
            qy = fill(qyv, coords.Ny1, coords.Nx1)
            compute_fluid_advection_heat!(Q, T, qx, qy, _ -> f; coords=coords)
            @test isapprox(Q[8, 10], -f * (qxv * gx + qyv * gy); rtol=1e-10)
            @test isapprox(Q[3, 17], Q[8, 10]; rtol=1e-10)
        end
    end

    @testset "zero flux and uniform temperature give zero heat" begin
        T = _linear_temperature(coords, 500.0, 0.05, 0.02)
        compute_fluid_advection_heat!(Q, T, zeros(size(Q)), zeros(size(Q)), _ -> f; coords=coords)
        @test all(iszero, Q)
        qx = fill(1.0e-8, size(Q))
        compute_fluid_advection_heat!(Q, fill(500.0, size(Q)), qx, qx, _ -> f; coords=coords)
        @test all(iszero, Q)
    end

    @testset "a face takes the upstream temperature, ghost nodes stay zero" begin
        T = fill(300.0, size(Q))
        i0, j0 = 8, 10
        T[i0, j0 + 1] = 1300.0
        qx = fill(1.0e-8, size(Q))
        # Flow runs from j toward j+1: the hot downstream node changes nothing upstream.
        compute_fluid_advection_heat!(Q, T, qx, zeros(size(Q)), _ -> f; coords=coords)
        @test isapprox(Q[i0, j0], 0.0; atol=1e-20)
        # The hot node itself cools from inflow of 300 K fluid.
        @test isapprox(Q[i0, j0 + 1], f * 1.0e-8 * (300.0 - 1300.0) / coords.dx; rtol=1e-12)
        @test isapprox(Q[i0, j0 + 2], f * 1.0e-8 * (1300.0 - 300.0) / coords.dx; rtol=1e-12)
        @test all(iszero, (Q[1, :], Q[end, :], Q[:, 1], Q[:, end]))
    end

    @testset "heat capacity is read at the upstream temperature" begin
        T = fill(300.0, size(Q))
        T[8, 9] = 400.0
        qx = fill(1.0e-8, size(Q))
        compute_fluid_advection_heat!(Q, T, qx, zeros(size(Q)), t -> t * 1.0e4; coords=coords)
        @test isapprox(Q[8, 10], 400.0e4 * 1.0e-8 * 100.0 / coords.dx; rtol=1e-12)
        # The next node downstream takes 300 K fluid at the capacity of 300 K.
        @test isapprox(Q[8, 9], 300.0e4 * 1.0e-8 * (300.0 - 400.0) / coords.dx; rtol=1e-12)
    end

    @testset "rough flow with div q: no new extremum over 200 steps" begin
        rng = MersenneTwister(7)
        T = 300.0 .+ 1200.0 .* rand(rng, size(Q)...)
        T0min, T0max = extrema(T[2:(end - 1), 2:(end - 1)])
        qx = 1.0e-8 .* randn(rng, size(Q)...)
        qy = 1.0e-8 .* randn(rng, size(Q)...)
        C = 3.0e6
        # Courant sum over at most 4 inflow faces stays below 1 with this step.
        dt = 0.2 * C * min(coords.dx, coords.dy) / (f * maximum(abs, vcat(vec(qx), vec(qy))))
        T[1, :] .= 800.0; T[end, :] .= 800.0; T[:, 1] .= 800.0; T[:, end] .= 800.0
        for _ in 1:200
            compute_fluid_advection_heat!(Q, T, qx, qy, _ -> f; coords=coords)
            T .+= dt .* Q ./ C
        end
        Tmin, Tmax = extrema(T[2:(end - 1), 2:(end - 1)])
        @test isapprox(max(Tmax - T0max, 0.0), 0.0; atol=1e-9)
        @test isapprox(max(T0min - Tmin, 0.0), 0.0; atol=1e-9)
        # A rough field with div q must change: guards a term that does nothing.
        @test !isapprox(Tmax - Tmin, T0max - T0min; rtol=1e-3)
    end

    @testset "steady advection-conduction converges to the analytic profile at first order" begin
        errs = [_steady_profile_error(n, 5.0) for n in (25, 50, 100)]
        @test isapprox(log2(errs[1] / errs[2]), 1.0; atol=0.15)
        @test isapprox(log2(errs[2] / errs[3]), 1.0; atol=0.1)
        @test isapprox(errs[3], 0.0; atol=0.012)
    end

    @testset "assembled thermal system adds the term once and closes its heat sum" begin
        T = _linear_temperature(coords, 500.0, 0.05, 0.02)
        qx = fill(2.0e-8, size(Q))
        qy = fill(-1.0e-8, size(Q))
        compute_fluid_advection_heat!(Q, T, qx, qy, _ -> f; coords=coords)
        Z = zeros(size(Q))
        RHOCP = fill(3.3e6, size(Q))
        dt = 1.0e9
        RT = zeros(coords.Ny1 * coords.Nx1)
        LT = assemble_thermal_lse!(T, RHOCP, Z, Z, Z, Z, Z, Z, RT, dt; coords=coords, Q_fluid=Q)
        T2 = reshape(LT.cscmatrix \ RT, size(Q))
        I2 = 2:(coords.Ny1 - 1)
        J2 = 2:(coords.Nx1 - 1)
        @test isapprox(
            sum(RHOCP[I2, J2] .* (T2[I2, J2] .- T[I2, J2])), sum(Q[I2, J2]) * dt; rtol=1e-8
        )
        # The same system without the term leaves temperature unchanged.
        RT0 = zeros(length(RT))
        LT0 = assemble_thermal_lse!(T, RHOCP, Z, Z, Z, Z, Z, Z, RT0, dt; coords=coords)
        @test isapprox(reshape(LT0.cscmatrix \ RT0, size(Q))[I2, J2], T[I2, J2]; rtol=1e-10)
    end

    @testset "input contracts" begin
        T = fill(300.0, size(Q))
        qx = zeros(size(Q))
        @test_throws ArgumentError compute_fluid_advection_heat!(
            zeros(3, 3), T, qx, qx, _ -> f; coords=coords
        )
        @test_throws ArgumentError compute_fluid_advection_heat!(
            Q, T, zeros(3, 3), qx, _ -> f; coords=coords
        )
        @test_throws DomainError compute_fluid_advection_heat!(
            Q, T, fill(1.0e-8, size(Q)), qx, _ -> -1.0; coords=coords
        )
    end
end

@testset "Thermal solve holds the fluid term built from the Darcy fluxes" begin
    path = mktempdir()
    cfg = override_config(
        load_config(joinpath(@__DIR__, "..", "configs", "lunar_growth_tutorial.toml")),
        Dict{String,Any}("time.n_steps" => 3),
    )
    st = Erebus.simulation_loop(cfg; output_path=path)
    g = st.grids
    # The domain can telescope to more nodes at the configured spacing.
    base = Erebus.GridCoordinates(cfg.grid)
    Ny, Nx = size(g.tk1) .- 1
    coords = Erebus.GridCoordinates(Nx, Ny; xsize=base.dx * (Nx - 1), ysize=base.dy * (Ny - 1))
    Q = zeros(size(g.Q_fluid_grid))
    compute_fluid_advection_heat!(
        Q,
        g.tk1,
        g.qxD,
        g.qyD,
        t -> Erebus.compute_rhocpfluidm(t, Erebus.marker_property_mode; rhocpfluidm=cfg.materials.rhocpfluidm);
        coords=coords,
    )
    @test !all(iszero, Q)
    @test isapprox(g.Q_fluid_grid, Q; rtol=1e-12)
end

@testset "Thermal solve right side holds the fluid term once" begin
    cfg = load_config(joinpath(@__DIR__, "..", "configs", "lunar_growth_tutorial.toml"))
    st, coords, ws, _, _ = init_simulation(cfg; output_path=mktempdir())
    Erebus.interpolate_markers_to_grid!(
        st,
        coords,
        cfg;
        p2m_workspace=ws.p2m,
        thread_buffers=ws.thread_buffers,
        interp_arrays=ws.interp_arrays,
    )
    Erebus.recompute_bulk_viscosity!(
        st.grids.ETA, st.grids.ETAP, st.grids.ETAPHI, st.grids.PHI, cfg.solver.etaphikoef
    )
    Erebus.snapshot_step_start_inventories!(ws, st, cfg)
    g = st.grids
    g.tk1 .= [300.0 + 0.002 * coords.xp[j] for i in 1:(coords.Ny1), j in 1:(coords.Nx1)]
    g.qxD .= 1.0e-10
    fill!(g.qyD, 0.0)
    dt_solve = st.dt
    Erebus.solve_thermal_energy!(st, coords, cfg, ws; titer=1, DHP_pyro=nothing)
    @test !all(iszero, g.Q_fluid_grid)
    RT_without = zeros(length(ws.RT))
    assemble_thermal_lse!(
        g.tk1,
        g.RHOCP,
        g.KX,
        g.KY,
        g.HR,
        g.HA,
        g.HS,
        g.DHP,
        RT_without,
        dt_solve;
        coords=coords,
        Q_metric=g.Q_metric,
        Q_lat=g.Q_lat_grid,
        Q_seg=g.Q_seg_grid,
    )
    @test isapprox(ws.RT .- RT_without, vec(g.Q_fluid_grid); rtol=1e-4)
end
