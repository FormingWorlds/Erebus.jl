using Erebus
using LinearAlgebra
using Random
using SparseArrays
using Test

isdefined(@__MODULE__, :PorousShellValidation) ||
    include(joinpath(@__DIR__, "..", "benchmarks", "validation_porous_shell.jl"))
FS_PSV = PorousShellValidation

# Rock conductivity [W/(m K)], heat capacity [J/(m^3 K)], heating [W/m^3], surface T [K]
const FS_K = 2.0
const FS_RHOCP = 3.0e6
const FS_H = 1.0e-6
const FS_TS = 250.0
const FS_R = 20.0e3
const FS_L = 50.0e3

"""
Grid of `N` x `N` basic nodes in a box of `FS_L`, sphere of radius `FS_R` centred at the
box centre plus `shift` cells in x and y, and its fixed-surface mask.
"""
function fs_grid(N::Int, shift::Float64)
    coords = GridCoordinates(N, N; xsize=FS_L, ysize=FS_L)
    xc = FS_L / 2 + shift * coords.dx
    yc = FS_L / 2 + shift * coords.dy
    fixed = zeros(Bool, coords.Ny1, coords.Nx1)
    nfix = Erebus.compute_fixed_surface_mask!(
        fixed, coords; xcenter=xc, ycenter=yc, rplanet=FS_R
    )
    return coords, xc, yc, fixed, nfix
end

"""
Take `nsteps` backward-Euler steps of `dt` [s] of the spherical-metric heat equation with
uniform rock properties and heating `H` inside `FS_R` and `FS_TS` at the fixed nodes.
`cut_cell = false` keeps the node staircase. `rhocp_fault` multiplies the heat capacity of
the free nodes next to a fixed node. Returns the fields of the last step.
"""
function fs_run!(
    tk, coords, xc, yc, fixed; dt, nsteps, H=FS_H, cut_cell=true, rhocp_fault=1.0
)
    Ny1, Nx1 = coords.Ny1, coords.Nx1
    KX = fill(FS_K, Ny1, Nx1)
    KY = fill(FS_K, Ny1, Nx1)
    if cut_cell
        Erebus.apply_cut_cell_conductance!(
            KX, KY, fixed, coords; xcenter=xc, ycenter=yc, rplanet=FS_R
        )
    end
    RHOCP = fill(FS_RHOCP, Ny1, Nx1)
    for j in 2:(Nx1 - 1), i in 2:(Ny1 - 1)
        fixed[i, j] && continue
        if fixed[i - 1, j] || fixed[i + 1, j] || fixed[i, j - 1] || fixed[i, j + 1]
            RHOCP[i, j] *= rhocp_fault
        end
    end
    HR = fill(H, Ny1, Nx1)
    Z = zeros(Ny1, Nx1)
    Q = zeros(Ny1, Nx1)
    RT = zeros(Ny1 * Nx1)
    tk[fixed] .= FS_TS
    tk_start = copy(tk)
    for _ in 1:nsteps
        tk_start .= tk
        Erebus.compute_spherical_metric_heat_source!(
            Q, tk, KX, KY, coords; xcenter=xc, ycenter=yc, rplanet=FS_R, reg_cells=0.5
        )
        LT = Erebus.assemble_thermal_lse!(
            tk, RHOCP, KX, KY, HR, Z, Z, Z, RT, dt; coords=coords, Q_metric=Q, fixed=fixed
        )
        tk .= reshape(LT.cscmatrix \ RT, Ny1, Nx1)
    end
    return (;
        tk1=tk_start,
        tk2=copy(tk),
        RHOCP,
        KX,
        KY,
        HR,
        HA=Z,
        HS=Z,
        DHP=Z,
        Q_metric=Q,
        Q_lat_grid=Z,
        Q_seg_grid=Z,
        Q_fluid_grid=Z,
    )
end

"""
Effective radius of a steady sphere from the free nodes with `2 dx < r < R/2`:
`R_eff^2 = 6 k (T_n - T_s) / H + r_n^2`, averaged over the nodes.
"""
function fs_reff(tk, coords, xc, yc, fixed)
    vals = Float64[]
    for j in 2:(coords.Nx1 - 1), i in 2:(coords.Ny1 - 1)
        fixed[i, j] && continue
        r = hypot(coords.xp[j] - xc, coords.yp[i] - yc)
        (2 * coords.dx < r < FS_R / 2) || continue
        push!(vals, sqrt(6 * FS_K * (tk[i, j] - FS_TS) / FS_H + r^2))
    end
    return sum(vals) / length(vals)
end

function fs_steady_error(N, shift; cut_cell=true)
    coords, xc, yc, fixed, _ = fs_grid(N, shift)
    tk = fill(FS_TS, coords.Ny1, coords.Nx1)
    fs_run!(tk, coords, xc, yc, fixed; dt=1.0e30, nsteps=80, cut_cell=cut_cell)
    return abs(fs_reff(tk, coords, xc, yc, fixed) - FS_R), coords.dx
end

# Configuration of the porous shell case with the thermodynamics field replaced.
function fs_config(cfg::SimulationConfig; surface_temperature=:fixed, kwargs...)
    th = cfg.thermodynamics
    th_new = ThermalConfig(;
        (f => getfield(th, f) for f in fieldnames(ThermalConfig))...,
        surface_temperature=surface_temperature,
    )
    return SimulationConfig(;
        (f => getfield(cfg, f) for f in fieldnames(SimulationConfig))...,
        thermodynamics=th_new,
        kwargs...,
    )
end

function fs_case(; Nx=17, n_steps=2)
    return FS_PSV.make_porous_shell_config(;
        Nx=Nx, phi_init=0.25, t_form_yr=2.5e6, n_steps=n_steps, dt_yr=5.0e4
    )
end

@testset "Fixed surface temperature (thermodynamics.surface_temperature = :fixed)" begin
    @testset "S7 config: the option and its refusals" begin
        base = fs_config(
            default_config();
            surface_temperature=:sticky_air,
            disk=DiskConfig(; enabled=true),
        )
        @test default_config().thermodynamics.surface_temperature === :sticky_air
        valid = try
            Erebus.validate_config(fs_config(base))
            true
        catch
            false
        end
        @test valid
        @test_throws ArgumentError Erebus.validate_config(
            fs_config(base; surface_temperature=:radiative)
        )
        no_disk = fs_config(base; disk=DiskConfig(; enabled=false))
        @test_throws ArgumentError Erebus.validate_config(no_disk)
        for (name, kw) in (
            ("accretion", (; accretion=AccretionConfig(; active=true))),
            ("telescoping", (; telescoping=TelescopingConfig(; active=true))),
            ("venting", (; venting=VentingConfig(; active=true))),
            ("atmosphere", (; atmosphere=AtmosphereConfig(; active=true))),
        )
            err = try
                Erebus.validate_config(fs_config(base; kw...))
                nothing
            catch e
                e
            end
            @test err isa ArgumentError
            @test occursin("surface_temperature", sprint(showerror, err))
        end
        th_rad = ThermalConfig(;
            (f => getfield(base.thermodynamics, f) for f in fieldnames(ThermalConfig))...,
            surface_radiation=true,
            surface_temperature=:fixed,
        )
        cfg_rad = SimulationConfig(;
            (f => getfield(base, f) for f in fieldnames(SimulationConfig))...,
            thermodynamics=th_rad,
        )
        @test_throws ArgumentError Erebus.validate_config(cfg_rad)
        # The default stays valid and the option round-trips through the TOML dict.
        d = Erebus.config_to_dict(fs_config(base))
        @test d["thermodynamics"]["surface_temperature"] in (:fixed, "fixed")
    end

    @testset "Mask: nodes outside R_p are fixed, box edges are not" begin
        coords, xc, yc, fixed, nfix = fs_grid(33, 0.13)
        @test nfix == count(fixed)
        @test !any(fixed[1, :]) && !any(fixed[end, :])
        @test !any(fixed[:, 1]) && !any(fixed[:, end])
        n_out = 0
        for j in 2:(coords.Nx1 - 1), i in 2:(coords.Ny1 - 1)
            n_out += hypot(coords.xp[j] - xc, coords.yp[i] - yc) > FS_R
        end
        @test nfix == n_out
        @test 0 < nfix < (coords.Nx1 - 2) * (coords.Ny1 - 2)
        @test_throws DomainError Erebus.compute_fixed_surface_mask!(
            fixed, coords; xcenter=xc, ycenter=yc, rplanet=-1.0
        )
        @test_throws ArgumentError Erebus.compute_fixed_surface_mask!(
            zeros(Bool, 3, 3), coords; xcenter=xc, ycenter=yc, rplanet=FS_R
        )
    end

    @testset "Cut-cell face: k dx / d with d the distance to R_p" begin
        coords, xc, yc, fixed, _ = fs_grid(33, 0.37)
        KX = fill(FS_K, coords.Ny1, coords.Nx1)
        KY = fill(FS_K, coords.Ny1, coords.Nx1)
        Erebus.apply_cut_cell_conductance!(
            KX, KY, fixed, coords; xcenter=xc, ycenter=yc, rplanet=FS_R
        )
        # The face east of the last inside node on the row through the centre
        i = argmin(abs.(coords.yp .- yc))
        j = findlast(jj -> !fixed[i, jj] && fixed[i, jj + 1], 1:(coords.Nx1 - 1))
        dyc = coords.yp[i] - yc
        d = xc + sqrt(FS_R^2 - dyc^2) - coords.xp[j]
        @test 0.0 < d <= coords.dx
        @test isapprox(KX[i, j], FS_K * coords.dx / max(d, 1.0e-3 * coords.dx); rtol=1e-12)
        # Faces between 2 free nodes and between 2 fixed nodes keep their value
        @test isapprox(KX[i, j - 1], FS_K; rtol=0, atol=0)
        @test isapprox(KX[i, j + 1], FS_K; rtol=0, atol=0)
        @test all(>=(FS_K * (1 - 1e-12)), KX) && all(>=(FS_K * (1 - 1e-12)), KY)
        @test maximum(KX) <= 1000 * FS_K * (1 + 1e-12)
    end

    @testset "S1 steady sphere: R_eff converges at second order" begin
        for shift in (0.0, 0.13, 0.37)
            errs = Float64[]
            dxs = Float64[]
            for N in (33, 65, 129)
                e, dx = fs_steady_error(N, shift)
                push!(errs, e)
                push!(dxs, dx)
                println("S1 shift=$shift N=$N |R_eff-R_p|=$(e) m = $(e / dx) dx")
            end
            orders = [log2(errs[k] / errs[k + 1]) for k in 1:2]
            println("S1 shift=$shift orders=$orders")
            @test errs[3] < errs[1]
            @test errs[3] < 0.05 * dxs[3]
            # second order, not first (1) or a super-convergent fluke (> 3)
            @test 1.5 < sum(orders) / 2 < 3.0
        end
        # Without the cut-cell face the error is of the order of the staircase offset
        e_stair, dx65 = fs_steady_error(65, 0.37; cut_cell=false)
        e_cut, _ = fs_steady_error(65, 0.37)
        println("S1 staircase N=65 shift=0.37 |R_eff-R_p|=$(e_stair / dx65) dx")
        @test e_stair > 5 * e_cut
    end

    @testset "S2 fixed nodes hold and free rows equal the main assembly" begin
        coords, xc, yc, fixed, _ = fs_grid(33, 0.13)
        Ny1, Nx1 = coords.Ny1, coords.Nx1
        rng = MersenneTwister(7)
        tk1 = 300.0 .+ 500.0 .* rand(rng, Ny1, Nx1)
        tk1[fixed] .= FS_TS
        KX = 1.0 .+ rand(rng, Ny1, Nx1)
        KY = 1.0 .+ rand(rng, Ny1, Nx1)
        RHOCP = FS_RHOCP .* (1.0 .+ rand(rng, Ny1, Nx1))
        HR = fill(FS_H, Ny1, Nx1)
        Z = zeros(Ny1, Nx1)
        dt = 1.0e12
        RT_m = zeros(Ny1 * Nx1)
        RT_f = zeros(Ny1 * Nx1)
        L_main = Erebus.assemble_thermal_lse!(
            tk1, RHOCP, KX, KY, HR, Z, Z, Z, RT_m, dt; coords=coords
        ).cscmatrix
        L_fix = Erebus.assemble_thermal_lse!(
            tk1, RHOCP, KX, KY, HR, Z, Z, Z, RT_f, dt; coords=coords, fixed=fixed
        ).cscmatrix
        n_same = 0
        n_bad = 0
        for j in 2:(Nx1 - 1), i in 2:(Ny1 - 1)
            gk = (j - 1) * Ny1 + i
            if fixed[i, j]
                n_bad += !isapprox(L_fix[gk, gk], 1.0; atol=0)
                n_bad += !isapprox(sum(abs, L_fix[gk, :]), 1.0; atol=0)
                n_bad += !isapprox(RT_f[gk], FS_TS; atol=0)
            elseif !(
                fixed[i - 1, j] || fixed[i + 1, j] || fixed[i, j - 1] || fixed[i, j + 1]
            )
                n_same += 1
                n_bad += !isapprox(Vector(L_fix[gk, :]), Vector(L_main[gk, :]); atol=0)
                n_bad += !isapprox(RT_f[gk], RT_m[gk]; atol=0)
            end
        end
        @test 100 < n_same < (Nx1 - 2) * (Ny1 - 2)
        @test n_bad == 0
        tk2 = reshape(L_fix \ RT_f, Ny1, Nx1)
        @test all(isapprox.(tk2[fixed], FS_TS; atol=1e-9))
        @test all(isfinite, tk2)
        # Matrix-free operator: same product and unit diagonal at the fixed nodes
        op = MatrixFreeThermalOperator(RHOCP, KX, KY, dt; coords=coords, fixed=fixed)
        x = rand(rng, Ny1 * Nx1)
        y = zeros(Ny1 * Nx1)
        Erebus.apply_thermal_A!(y, op, x)
        @test isapprox(y, L_fix * x; rtol=1e-12)
        dgl = Erebus.compute_operator_diagonal(op)
        @test isapprox(dgl, Vector(diag(L_fix)); rtol=1e-12)
        ST = zeros(Ny1 * Nx1)
        Erebus.solve_thermal_iterative!(op, RT_f, ST; rtol=1e-12, atol=1e-12, maxiter=2000)
        @test isapprox(ST, vec(tk2); rtol=1e-6)
        @test all(isapprox.(reshape(ST, Ny1, Nx1)[fixed], FS_TS; atol=1e-6))
        # A node 1e-5 dx inside R_p: the face takes the floor d = 1e-3 dx (1000 k)
        jc = (Nx1 + 1) ÷ 2
        xc = coords.xp[jc]
        yc = coords.yp[jc]
        R = coords.xp[jc + 10] - xc + 1.0e-5 * coords.dx
        fixed2 = zeros(Bool, Ny1, Nx1)
        Erebus.compute_fixed_surface_mask!(
            fixed2, coords; xcenter=xc, ycenter=yc, rplanet=R
        )
        K2x = fill(FS_K, Ny1, Nx1)
        K2y = fill(FS_K, Ny1, Nx1)
        Erebus.apply_cut_cell_conductance!(
            K2x, K2y, fixed2, coords; xcenter=xc, ycenter=yc, rplanet=R
        )
        @test isapprox(K2x[jc, jc + 10], 1000 * FS_K; rtol=1e-12)
        tk1[fixed2] .= FS_TS
        tk1[.!fixed2] .= 600.0
        RT2 = zeros(Ny1 * Nx1)
        L2 = Erebus.assemble_thermal_lse!(
            tk1, RHOCP, K2x, K2y, HR, Z, Z, Z, RT2, dt; coords=coords, fixed=fixed2
        ).cscmatrix
        x_direct = L2 \ RT2
        op2 = MatrixFreeThermalOperator(RHOCP, K2x, K2y, dt; coords=coords, fixed=fixed2)
        ST2 = zeros(Ny1 * Nx1)
        _, st2 = Erebus.solve_thermal_iterative!(
            op2, RT2, ST2; rtol=1e-12, atol=1e-12, maxiter=2000
        )
        @test st2.solved
        @test isapprox(ST2, x_direct; rtol=1e-6)
    end

    @testset "S8 transient heating against the eigenmode series" begin
        coords, xc, yc, fixed, _ = fs_grid(65, 0.13)
        kappa = FS_K / FS_RHOCP
        tau = FS_R^2 / kappa
        dt = 0.005 * tau
        n = 40
        function max_rel_err(rhocp_fault)
            tk = fill(FS_TS, coords.Ny1, coords.Nx1)
            fs_run!(tk, coords, xc, yc, fixed; dt=dt, nsteps=n, rhocp_fault=rhocp_fault)
            err = 0.0
            rise_max = 0.0
            for j in 2:(coords.Nx1 - 1), i in 2:(coords.Ny1 - 1)
                fixed[i, j] && continue
                r = hypot(coords.xp[j] - xc, coords.yp[i] - yc)
                r < 0.8 * FS_R || continue
                T_ref = FS_PSV.eigenmode_T(
                    r;
                    R=FS_R,
                    T_surf=FS_TS,
                    kappa=kappa,
                    rhocp=FS_RHOCP,
                    H=fill(FS_H, n),
                    dt=dt,
                )
                err = max(err, abs(tk[i, j] - T_ref))
                rise_max = max(rise_max, T_ref - FS_TS)
            end
            return err / rise_max, rise_max
        end
        e_ok, rise = max_rel_err(1.0)
        e_fault, _ = max_rel_err(1000.0)
        println("S8 N=65 max |T - T_series| / rise = $e_ok, with air capacity = $e_fault")
        @test 1.0 < rise < 1.0e4
        @test e_ok < 0.02
        @test e_fault > 5 * e_ok
    end

    @testset "S5 energy bookkeeping with the boundary flux" begin
        coords, xc, yc, fixed, _ = fs_grid(33, 0.37)
        tk = fill(FS_TS, coords.Ny1, coords.Nx1)
        dt = 0.01 * FS_R^2 * FS_RHOCP / FS_K
        g = fs_run!(tk, coords, xc, yc, fixed; dt=dt, nsteps=3)
        res = FS_PSV.thermal_step_energy_residual(g, dt, coords.dx, coords.dy; fixed=fixed)
        @test res < 1e-10
        @test isfinite(res)
        # Without the mask the fixed nodes and the boundary flux break the identity
        res_plain = FS_PSV.thermal_step_energy_residual(g, dt, coords.dx, coords.dy)
        @test 1e-3 < res_plain <= 1.0
    end

    @testset "S4 subgrid diffusion and limit at fixed nodes, marker heat identity (D3)" begin
        coords, xc, yc, fixed, _ = fs_grid(17, 0.13)
        Ny1, Nx1 = coords.Ny1, coords.Nx1
        rng = MersenneTwister(11)
        marknum = 4000
        xm = coords.dx .+ rand(rng, marknum) .* (FS_L - 2 * coords.dx)
        ym = coords.dy .+ rand(rng, marknum) .* (FS_L - 2 * coords.dy)
        tm = [hypot(xm[m] - xc, ym[m] - yc) <= FS_R ? 1 : 3 for m in 1:marknum]
        tk1 = 400.0 .+ 200.0 .* rand(rng, Ny1, Nx1)
        tk1[fixed] .= FS_TS
        tkm = 400.0 .+ 200.0 .* rand(rng, marknum)
        rhocpm = FS_RHOCP .* (1.0 .+ rand(rng, marknum))
        km = fill(FS_K, marknum)
        phim = fill(0.1, marknum)
        DT_solve = 5.0 .* randn(rng, Ny1, Nx1)
        DT_solve[fixed] .= 0.0
        DT = copy(DT_solve)
        TKSUM = zeros(Ny1, Nx1)
        RHOCPSUM = zeros(Ny1, Nx1)
        tkm0 = copy(tkm)
        Erebus.apply_subgrid_temperature_diffusion!(
            xm,
            ym,
            tm,
            tkm,
            phim,
            tk1,
            DT,
            TKSUM,
            RHOCPSUM,
            1.0e12,
            marknum,
            2;
            coords=coords,
            dsubgridt=1.0,
            rhocptotalm=rhocpm,
            ktotalm=km,
            fixed=fixed,
        )
        @test all(isapprox.(DT[fixed], 0.0; atol=0))
        @test any(>(0.0), RHOCPSUM[fixed])
        # Marker increment from the grid
        for m in 1:marknum
            tm[m] < 3 || continue
            i, j, w = Erebus.fix_weights(
                xm[m],
                ym[m],
                coords.xp,
                coords.yp,
                coords.dx,
                coords.dy,
                1,
                Nx1 - 1,
                1,
                Ny1 - 1,
            )
            tkm[m] += Erebus.dot4(Erebus.grid_vector(i, j, DT), w)
        end
        marker_heat = sum(rhocpm[m] * (tkm[m] - tkm0[m]) for m in 1:marknum if tm[m] < 3)
        grid_heat = sum(DT_solve[.!fixed] .* RHOCPSUM[.!fixed])
        boundary_share = sum(TKSUM[fixed])
        scale = sum(abs, DT_solve .* RHOCPSUM) + sum(abs, TKSUM)
        @test isapprox(marker_heat, grid_heat + boundary_share; atol=1e-12 * scale)
        @test abs(boundary_share) > 1e-6 * scale
        # The limit gives no capacity at fixed nodes; the heat it moves stays bounded
        tkm_l = copy(tkm0)
        DEF = zeros(Ny1, Nx1)
        CAP = zeros(Ny1, Nx1)
        tk2 = tk1 .+ DT_solve
        st = Erebus.limit_marker_temperature_increment!(
            xm,
            ym,
            tm,
            tkm_l,
            rhocpm,
            DT_solve,
            tk2,
            DEF,
            CAP,
            marknum;
            coords=coords,
            fixed=fixed,
        )
        @test all(isapprox.(CAP[fixed], 0.0; atol=0))
        moved = sum(rhocpm[m] * (tkm_l[m] - tkm0[m]) for m in 1:marknum if tm[m] < 3)
        @test abs(moved) <= st.heat_lost * (1 + 1e-9) + 1e-9 * st.heat_rejected
        @test st.n_limited >= 0 && all(isfinite, tkm_l)
    end

    @testset "S3 rock-only node properties on the scatter paths" begin
        cfg_fixed = fs_config(fs_case())
        results = Dict{Symbol,Any}()
        for mode in (:tiled, :buffered)
            cfg = fs_config(
                cfg_fixed;
                surface_temperature=:fixed,
                solver=SolverConfig(;
                    (
                        f => getfield(cfg_fixed.solver, f) for f in fieldnames(SolverConfig)
                    )...,
                    p2m_mode=mode,
                ),
            )
            state, coords, ws, _, _ = init_simulation(cfg; output_path=mktempdir())
            Erebus.prepare_step_ambient!(state, coords, cfg, ws)
            interpolate_markers_to_grid!(
                state,
                coords,
                cfg;
                p2m_workspace=ws.p2m,
                thread_buffers=ws.thread_buffers,
                interp_arrays=ws.interp_arrays,
            )
            results[mode] = (; g=deepcopy(state.grids), state, coords, cfg)
        end
        g = results[:tiled].g
        @test g.fixed_tk isa Matrix{Bool}
        @test isapprox(g.RHOCP, results[:buffered].g.RHOCP; rtol=1e-12)
        @test isapprox(g.tk1, results[:buffered].g.tk1; rtol=1e-12)
        @test isapprox(g.KX, results[:buffered].g.KX; rtol=1e-12)
        # Brute force at the inside nodes with air markers in their stencil
        (; state, coords) = results[:tiled]
        core = state.markers.core
        acc = state.accumulators
        Ny1, Nx1 = coords.Ny1, coords.Nx1
        W = zeros(Ny1, Nx1)
        Wr = zeros(Ny1, Nx1)
        C = zeros(Ny1, Nx1)
        E = zeros(Ny1, Nx1)
        for m in 1:length(state.markers)
            i, j, w = Erebus.fix_weights(
                core.xm[m],
                core.ym[m],
                coords.xp,
                coords.yp,
                coords.dx,
                coords.dy,
                1,
                Nx1 - 1,
                1,
                Ny1 - 1,
            )
            Erebus.interpolate_add_to_grid!(i, j, w, 1.0, W)
            if core.tm[m] < 3
                Erebus.interpolate_add_to_grid!(i, j, w, 1.0, Wr)
                Erebus.interpolate_add_to_grid!(i, j, w, core.rhocptotalm[m], C)
                Erebus.interpolate_add_to_grid!(
                    i, j, w, core.rhocptotalm[m] * core.tkm[m], E
                )
            end
        end
        n_mixed = 0
        n_bad = 0
        for j in 2:(Nx1 - 1), i in 2:(Ny1 - 1)
            g.fixed_tk[i, j] && continue
            (Wr[i, j] < W[i, j] * (1 - 1e-9)) || continue
            n_mixed += 1
            n_bad += !isapprox(g.RHOCP[i, j], C[i, j] / W[i, j]; rtol=1e-10)
            n_bad += !isapprox(g.tk1[i, j], E[i, j] / C[i, j]; rtol=1e-10)
        end
        @test 0 < n_mixed < count(.!g.fixed_tk)
        @test n_bad == 0
        T_space = FS_PSV.T_SPACE
        @test all(isapprox.(g.tk1[g.fixed_tk], T_space; atol=1e-9))
    end

    @testset "S2/S6 short runs: fixed nodes hold, default unchanged" begin
        base = fs_case(; Nx=17, n_steps=2)
        st_fixed = mktempdir() do d
            simulation_loop(fs_config(base); output_path=d)
        end
        g = st_fixed.grids
        @test all(isapprox.(g.tk2[g.fixed_tk], FS_PSV.T_SPACE; atol=1e-9))
        @test all(isapprox.(g.DT0[g.fixed_tk], 0.0; atol=1e-9))
        @test all(isfinite, g.tk2)
        st_def = mktempdir() do d
            simulation_loop(base; output_path=d)
        end
        st_explicit = mktempdir() do d
            simulation_loop(fs_config(base; surface_temperature=:sticky_air); output_path=d)
        end
        @test st_def.grids.fixed_tk === nothing
        @test isequal(st_def.grids.tk2, st_explicit.grids.tk2)
        @test isequal(st_def.markers.core.tkm, st_explicit.markers.core.tkm)
    end
end
