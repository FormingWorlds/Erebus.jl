using LinearAlgebra: Tridiagonal
using SHA
using Test
using TOML
using Erebus

isdefined(@__MODULE__, :SphericalPressureDiffusionValidation) || include(
    joinpath(@__DIR__, "..", "benchmarks", "validation_spherical_pressure_diffusion.jl")
)

include(joinpath(@__DIR__, "..", "tools", "check_validation.jl"))
isdefined(@__MODULE__, :SphericalPressureDiffusionExport) || include(
    joinpath(@__DIR__, "..", "tools", "export_validation_spherical_pressure_diffusion.jl"),
)

const SPDV = SphericalPressureDiffusionValidation
const SPDE = SphericalPressureDiffusionExport

"""
    DiffusionInterruptingBench

Benchmark stand-in whose every property access throws `InterruptException`, as a Ctrl-C
during the benchmark does.
"""
struct DiffusionInterruptingBench end
Base.getproperty(::DiffusionInterruptingBench, ::Symbol) = throw(InterruptException())

"""
    DiffusionCrashingBench

Benchmark stand-in whose every property access throws an `ErrorException`, as a solver
failure during the benchmark does.
"""
struct DiffusionCrashingBench end
function Base.getproperty(::DiffusionCrashingBench, ::Symbol)
    return error("solver failure inside the benchmark")
end

"""
    write_diffusion_toml(path, dict)

Write `dict` as TOML to `path` and return `path`.

Parameters
----------
path : String
    Output file path.
dict : AbstractDict
    Table to write.

Returns
-------
String
    `path`.
"""
function write_diffusion_toml(path::String, dict::AbstractDict)
    open(path, "w") do io
        return TOML.print(io, dict)
    end
    return path
end

"""
    capture_diffusion_stdio(f)

Call `f()` with standard output and standard error redirected to files and return the
result with the captured text.

Parameters
----------
f : Function
    Function of no arguments.

Returns
-------
NamedTuple
    `value`: the return value of `f`; `out`, `err`: the captured standard output and
    standard error.
"""
function capture_diffusion_stdio(f)
    return mktempdir() do dir
        out_path, err_path = joinpath(dir, "out.txt"), joinpath(dir, "err.txt")
        value = open(out_path, "w") do out
            open(err_path, "w") do err
                return redirect_stdio(f; stdout=out, stderr=err)
            end
        end
        return (; value, out=read(out_path, String), err=read(err_path, String))
    end
end

"""
    diffusion_coords(N; shift=0.0)

Return a `NamedTuple` with the staggered P node coordinates that the error measures read,
laid out as `GridCoordinates` does for an `N` x `N` grid of the box `SPDV.L_BOX`, with every
coordinate moved by `shift` [m] in x and y.

Parameters
----------
N : Int
    Basic grid resolution in x and y.
shift : Float64, optional
    Offset of all coordinates [m] (default `0.0`).

Returns
-------
NamedTuple
    `Nx`, `Ny`, `Nx1`, `Ny1`, `dx`, `dy`, `xp` and `yp`.
"""
function diffusion_coords(N::Int; shift::Float64=0.0)
    d = SPDV.L_BOX / (N - 1)
    xp = collect(range(-d / 2, SPDV.L_BOX + d / 2; length=N + 1)) .+ shift
    return (; Nx=N, Ny=N, Nx1=N + 1, Ny1=N + 1, dx=d, dy=d, xp, yp=copy(xp))
end

"""
    radial_be_fd(d, m; n_cells=4000, dtau=SPDV.DTAU)

Return the cell centres and the nondimensional pressure after `m` backward-Euler steps of
the radial diffusion equation `dp/dtau = x^(1 - d) d/dx(x^(d - 1) dp/dx)` on `0 <= x <= 1`,
with `p = 1` at `tau = 0` and `p = 0` at `x = 1`, by a cell-centred finite-volume solve.

Parameters
----------
d : Int
    Flow dimension, 2 for the disk and 3 for the sphere.
m : Int
    Number of time steps.
n_cells : Int, optional
    Number of cells on `[0, 1]` (default 4000).
dtau : Float64, optional
    Nondimensional time step (default `SPDV.DTAU`).

Returns
-------
NamedTuple
    `x`: cell centres; `p`: pressure at the cell centres.
"""
function radial_be_fd(d::Int, m::Int; n_cells::Int=4000, dtau::Float64=SPDV.DTAU)
    h = 1.0 / n_cells
    faces = collect(range(0.0, 1.0; length=n_cells + 1))
    x = faces[1:(end - 1)] .+ h / 2
    vol = [(faces[i + 1]^d - faces[i]^d) / d for i in 1:n_cells]
    area = faces .^ (d - 1)
    # Conductance of each interior face and of the Dirichlet face at x = 1.
    cond = [area[i + 1] / h for i in 1:(n_cells - 1)]
    c_out = area[end] / (h / 2)
    p = ones(n_cells)
    for _ in 1:m
        lo = zeros(n_cells - 1)
        di = vol ./ dtau
        up = zeros(n_cells - 1)
        for i in 1:(n_cells - 1)
            di[i] += cond[i]
            di[i + 1] += cond[i]
            up[i] = -cond[i]
            lo[i] = -cond[i]
        end
        di[end] += c_out
        p = Tridiagonal(lo, di, up) \ (vol ./ dtau .* p)
    end
    return (; x, p)
end

@testset "Spherical Pressure Diffusion Validation" begin
    @testset "Bessel functions J0 and J1 against tabulated values" begin
        j0, j1 = SPDV._j0j1(10.0)
        @test isapprox(j0, -0.2459357645; atol=1.0e-10)
        @test isapprox(j1, 0.0434727462; atol=1.0e-10)
        j0_0, j1_0 = SPDV._j0j1(0.0)
        @test isapprox(j0_0, 1.0; rtol=1.0e-15)
        @test isapprox(j1_0, 0.0; atol=1.0e-15)
        # J0 is even and J1 odd in x.
        @test isapprox(SPDV._j0j1(-3.0)[1], SPDV._j0j1(3.0)[1]; rtol=1.0e-14)
        @test isapprox(SPDV._j0j1(-3.0)[2], -SPDV._j0j1(3.0)[2]; rtol=1.0e-14)
        # J0' = -J1: a wrong phase in the integrand breaks this.
        for x in (0.5, 2.0, 7.5, 40.0)
            h = 1.0e-5
            fd = (SPDV._j0j1(x + h)[1] - SPDV._j0j1(x - h)[1]) / (2h)
            @test isapprox(fd, -SPDV._j0j1(x)[2]; atol=1.0e-9)
        end
        @test_throws DomainError SPDV._j0j1(NaN)
        @test_throws DomainError SPDV._j0j1(Inf)
    end

    @testset "zeros of J0 by Newton iteration" begin
        z = SPDV._j0_zeros(5)
        @test isapprox(z[1], 2.404825557695773; rtol=1.0e-14)
        @test isapprox(z[2], 5.520078110286311; rtol=1.0e-14)
        @test isapprox(SPDV._j0j1(z[1])[2], 0.5191474973; atol=1.0e-10)
        @test all(abs(SPDV._j0j1(zk)[1]) < 1.0e-13 for zk in z)
        # Consecutive zeros are about pi apart and in increasing order.
        @test all(2.9 < dz < pi + 0.2 for dz in diff(z))
        @test length(SPDV.J0_ZEROS) == SPDV.DISK_MODES
        @test isapprox(SPDV.J0_ZEROS[3], z[3]; rtol=1.0e-15)
        @test isapprox(SPDV.J1_AT_ZEROS[1], 0.5191474973; atol=1.0e-10)
        @test_throws ArgumentError SPDV._j0_zeros(0)
    end

    @testset "exact pressure series: limits, time stepping form and input contracts" begin
        for d in (2, 3)
            # The surface holds p = 0, and the pressure falls from the centre outwards.
            @test isapprox(SPDV.exact_p(1.0, d, 20), 0.0; atol=1.0e-12)
            xs = range(0.0, 1.0; length=41)
            @test all(diff([SPDV.exact_p(x, d, 20) for x in xs]) .< 0)
            for m in (5, 20)
                @test all(0.0 < SPDV.exact_p(x, d, m) < 1.0 for x in xs[1:(end - 1)])
            end
        end
        # Continuum centre values at tau = 0.1 from the first 4 modes summed by hand:
        # 0.707101 for the sphere and 0.84836 for the disk.
        cont3 = SPDV.exact_p(0.0, 3, 20; be=false)
        @test isapprox(cont3, 0.707101; atol=2.0e-6)
        @test isapprox(SPDV.exact_p(0.0, 2, 20; be=false), 0.84836; atol=1.0e-5)
        # Backward Euler lags the continuum: the centre values differ by 0.4 %.
        be3 = SPDV.exact_p(0.0, 3, 20)
        @test isapprox(be3, 0.70998; atol=1.0e-5)
        @test isapprox(SPDV.exact_p(0.0, 2, 20), 0.84486; atol=1.0e-5)
        @test 0.003 < be3 / cont3 - 1 < 0.005
        # The sphere drains faster than the disk.
        @test SPDV.exact_p(0.0, 3, 20) < SPDV.exact_p(0.0, 2, 20)
        @test_throws DomainError SPDV.exact_p(-0.1, 3, 20)
        @test_throws DomainError SPDV.exact_p(1.1, 2, 20)
        @test_throws DomainError SPDV.exact_p(NaN, 2, 20)
        @test_throws DomainError SPDV.exact_p(0.5, 1, 20)
        @test_throws DomainError SPDV.exact_p(0.5, 3, 0)
        @test_throws DomainError SPDV.exact_p(0.5, 3, 2.5)
        @test_throws DomainError SPDV.exact_p(0.5, 3, 20; dtau=0.0)
        @test_throws DomainError SPDV.exact_p(0.5, 3, 20; dtau=Inf)
        # One step needs more modes than the tables hold: an error, not a truncated sum.
        @test_throws ErrorException SPDV.exact_p(0.5, 3, 1)
        @test_throws ErrorException SPDV.exact_p(0.5, 2, 1)
    end

    @testset "series against a 1-D radial backward-Euler finite-volume solve (d = $d)" for d in
                                                                                           (
        2, 3
    )
        for m in (5, 20)
            fd = radial_be_fd(d, m)
            idx = 1:40:length(fd.x)
            series = SPDV.exact_p_cells(fd.x[idx], d, m)
            @test maximum(abs, series .- fd.p[idx]) < 1.0e-4
            @test length(series) == length(idx)
        end
        # The continuum form is a different reference: the solve departs from it.
        fd20 = radial_be_fd(d, 20)
        @test !isapprox(SPDV.exact_p(fd20.x[1], d, 20; be=false), fd20.p[1]; atol=1.0e-3)
    end

    @testset "series modes stop at the tail bound" begin
        # The terms beyond the stop change the sum by less than 1e-12.
        for d in (2, 3), m in (5, 20), x in (0.0, 0.37, 0.9)
            full = SPDV._series_sum(x, d, m, SPDV.DTAU, true; tail=0.0)
            @test isapprox(SPDV.exact_p(x, d, m), full; atol=1.0e-12)
        end
        @test SPDV.n_modes(3, 20, SPDV.DTAU, true) < SPDV.n_modes(3, 5, SPDV.DTAU, true)
        @test SPDV.n_modes(2, 5, SPDV.DTAU, true) < SPDV.DISK_MODES
        @test_throws ErrorException SPDV.n_modes(3, 1, SPDV.DTAU, true)
        @test_throws DomainError SPDV.n_modes(4, 5, SPDV.DTAU, true)
    end

    @testset "storage coefficient of the Pt and Pf rows in the rigid limit" begin
        bphi, bs, bf, phi = SPDV.BETA_PHI, SPDV.BETA_S, SPDV.BETA_F, SPDV.PHI
        s = SPDV.storage_coefficient(bphi, bs)
        kbw = SPDV.biot_willis_closed_form(bphi, bs)
        @test isapprox(s, 1.3e-10; rtol=1.0e-14)
        @test isapprox(kbw, 0.6; rtol=1.0e-14)
        # S = bd kbw (1/ksk - kbw) with bd = (bphi + bs) / (1 - phi), kbw = 1 - bs / bd and
        # 1/ksk = 1 + phi (bf - bs) / (bd - bs), the coefficients of the 2 rows.
        bd = (bphi + bs) / (1 - phi)
        inv_ksk = 1 + phi * (bf - bs) / (bd - bs)
        @test isapprox(1 - bs / bd, kbw; rtol=1.0e-14)
        @test isapprox(bd * kbw * (inv_ksk - kbw), s; rtol=1.0e-14)
        # The solver helpers give the same row coefficients for the case constants.
        bd_src = Erebus.compute_drained_compressibility(bphi, phi, bs)
        kbw_src = Erebus.compute_biot_willis_coefficient(bd_src, bs)
        ksk_src = Erebus.compute_skempton_coefficient(bd_src, phi, bs, bf)
        @test isapprox(bd_src * kbw_src * (1 / ksk_src - kbw_src), s; rtol=1.0e-12)
        # Limits: an incompressible solid stores phi bf only and kbw = 1; a rigid pore space
        # (bphi = 0) gives kbw = phi.
        @test isapprox(SPDV.storage_coefficient(bphi, 0.0), phi * bf; rtol=1.0e-14)
        @test isapprox(SPDV.biot_willis_closed_form(bphi, 0.0), 1.0; rtol=1.0e-14)
        @test isapprox(SPDV.biot_willis_closed_form(0.0, bs), phi; rtol=1.0e-14)
        # The second solid compressibility changes S by a quarter, so the diffusion time does.
        @test isapprox(
            SPDV.storage_coefficient(bphi, SPDV.BETA_S_ALT), 9.7273e-11; rtol=1.0e-4
        )
        @test_throws DomainError SPDV.storage_coefficient(-1.0e-10, bs)
        @test_throws DomainError SPDV.storage_coefficient(0.0, 0.0)
        @test_throws DomainError SPDV.storage_coefficient(bphi, NaN)
        @test_throws DomainError SPDV.biot_willis_closed_form(bphi, -bs)
    end

    @testset "time step: tau = 0.1 in 20 steps, from the year length of the solver" begin
        dt = SPDV.time_step_s(SPDV.BETA_S)
        @test isapprox(dt, 1.625e7; rtol=1.0e-12)
        @test isapprox(
            SPDV.K_PERM * dt / (SPDV.MU_F * 1.3e-10 * SPDV.R_PLANET^2),
            SPDV.DTAU;
            rtol=1.0e-12,
        )
        cfg = SPDV.make_config(33; spherical_flow=true)
        @test isapprox(cfg.time.dt_initial * Erebus.SECONDS_PER_YEAR, dt; rtol=1.0e-14)
        @test isapprox(cfg.time.dt_initial * cfg.time.yearlength, dt; rtol=1.0e-14)
        alt = SPDV.make_config(33; spherical_flow=true, beta_s=SPDV.BETA_S_ALT)
        @test isapprox(
            alt.time.dt_initial / cfg.time.dt_initial, 9.7273e-11 / 1.3e-10; rtol=1.0e-4
        )
    end

    @testset "case set-up: configuration, domain of the storage formula, input contracts" begin
        cfg = SPDV.make_config(SPDV.NX_REFERENCE; spherical_flow=false)
        pe = cfg.poroelasticity
        # C5: bf >= bs, bphi + bs > 0 and phimin < phi < phimax, so no clamp acts.
        @test pe.betafluid >= pe.betasolid
        @test 0.0 < SPDV.BETA_PHI + pe.betasolid < 1.0e-9
        @test pe.phimin < cfg.thermodynamics.phim0 < pe.phimax
        @test !cfg.geometry.spherical_flow
        @test cfg.solver.darcy_elimination
        @test isapprox(cfg.venting.k_vent, SPDV.K_VENT; rtol=1.0e-15)
        p_vent = Erebus.compute_venting_pressure(
            SPDV.T_SURF, SPDV.P_AMB; species=:H2O, L_sub=cfg.venting.L_sublimation
        )
        @test isapprox(p_vent, SPDV.P_AMB; rtol=1.0e-12)
        @test_throws ArgumentError SPDV.make_config(SPDV.NX_MIN - 1; spherical_flow=true)
        @test_throws ArgumentError SPDV.make_config(33; spherical_flow=true, beta_s=0.0)
        @test_throws ArgumentError SPDV.make_config(33; spherical_flow=true, beta_s=NaN)
        @test_throws ArgumentError SPDV.make_config(
            33; spherical_flow=true, beta_s=2 * SPDV.BETA_F
        )
        @test_throws ArgumentError SPDV.run_case(SPDV.NX_MIN - 1; spherical_flow=true)
        @test_throws ArgumentError SPDV.evaluate_spherical_pressure_diffusion_benchmark(;
            nx_list=(17,)
        )
        @test_throws ArgumentError SPDV.evaluate_spherical_pressure_diffusion_benchmark(;
            nx_list=(33.0,)
        )
    end

    @testset "C2: the vent drain cap stays inactive on every vent cell (margin 10)" begin
        phi_avail = SPDV.PHI - 1.0e-4
        for n in SPDV.NX_LIST, sf in (false, true), bs in (SPDV.BETA_S, SPDV.BETA_S_ALT)
            cap = SPDV.drain_cap_ratio(n; spherical_flow=sf, beta_s=bs)
            @test cap.ratio <= 0.1
            # The staircase of the circle has 2 (Nx - 1) vent cells: 256 at Nx = 129.
            @test cap.n_vent == 2 * (n - 1)
        end
        # Without the face weights a vent cell drains k_vent / (mu dx^2) per open face; the
        # cell with the most open faces sets the ratio.
        n = 129
        cap = SPDV.drain_cap_ratio(n; spherical_flow=false)
        dx = SPDV.L_BOX / (n - 1)
        per_face =
            SPDV.K_VENT / (SPDV.MU_F * dx^2) * SPDV.P0 * SPDV.time_step_s(SPDV.BETA_S)
        @test isapprox(cap.ratio, cap.max_faces * per_face / phi_avail; rtol=1.0e-12)
        @test cap.max_faces in (1, 2, 3)
        # A vent conductance 100 times larger makes the cap act (ratio above 1).
        wide = SPDV.drain_cap_ratio(n; spherical_flow=false, k_vent=1.0e-9).ratio
        @test isapprox(wide, 100 * cap.ratio; rtol=1.0e-12)
        @test 1.0 < wide < 1.5
    end

    @testset "set_fields!: initial pressures, face resistances and rigid matrix" begin
        cfg = SPDV.make_config(33; spherical_flow=true)
        mktempdir() do out
            state, coords, _, tio, _ = SPDV._init_case(cfg, out)
            tio === nothing || close(tio)
            g = state.grids
            inside(i, j) =
                hypot(coords.xp[j] - SPDV.X_CENTRE, coords.yp[i] - SPDV.Y_CENTRE) <=
                SPDV.R_PLANET
            cells = [(i, j) for j in 1:(coords.Nx1), i in 1:(coords.Ny1)]
            @test all(
                isapprox(
                    g.pf[i, j], SPDV.P_AMB + (inside(i, j) ? SPDV.P0 : 0.0); rtol=1.0e-15
                ) for (i, j) in cells
            )
            @test isapprox(g.pr, g.pf; rtol=1.0e-15)
            @test isapprox(g.pr0, g.pr; rtol=1.0e-15)
            @test isapprox(g.pf0, g.pf; rtol=1.0e-15)
            n_rock_faces = count(
                j < coords.Nx1 && inside(i, j) && inside(i, j + 1) for (i, j) in cells
            )
            @test count(x -> isapprox(x, SPDV.MU_F / SPDV.K_PERM; rtol=1.0e-15), g.RX) ==
                n_rock_faces
            @test all(x -> isapprox(x, SPDV.BETA_PHI; rtol=1.0e-15), g.BETAPHI)
            # C4 before the solve, on the basic nodes: ETA >= 1e27 and G dt >= 1e6 ETA.
            dt = SPDV.time_step_s(SPDV.BETA_S)
            @test maximum(inv, g.ETA) <= 1.0e-27
            @test maximum(g.ETA ./ (g.GGG .* dt)) <= 1.0e-6
            @test isapprox(state.dt, dt; rtol=1.0e-14)
            @test state.timestep == 2
        end
    end

    # Runs at Nx = 33 and 65 only. The Nx = 129 and order metrics take their recorded values,
    # which `validate_record_dict` checks against their reference and tolerance.
    bench_runs = SPDV.evaluate_spherical_pressure_diffusion_benchmark(; nx_list=(33, 65))
    runs = bench_runs.runs
    measured = bench_runs.metrics

    @testset "solver converges, fields are finite, dt stays fixed (spherical_flow = $sf)" for sf in
                                                                                              (
        false, true
    )
        for key in
            ((33, sf, :base), (65, sf, :base), (33, sf, :elim_off), (33, sf, :beta_alt))
            r = runs[key]
            @test r.converged
            for field in (:pf, :pr, :qxD, :qyD)
                @test all(isfinite, getfield(r, field))
            end
            # C1: timestep 2 keeps the compressibilities; C3: the solve keeps dt.
            @test r.timestep == 2
            @test r.dt_dev <= 1.0e-12
        end
    end

    @testset "C4: rigid matrix, no compaction or elastic softening (spherical_flow = $sf)" for sf in
                                                                                               (
        false, true
    )
        for n in (33, 65)
            r = runs[(n, sf, :base)]
            # 1/(ETAPHI (1 - phi)) <= 1e-6 bd/dt on the interior P cells after the first solve.
            @test r.c4_compaction <= 1.0e-6
            # G dt >= 1e6 eta on the basic nodes (GGG, ETA) and the P cells (GGGP, ETAP).
            @test r.c4_elastic <= 1.0e-6
            @test r.c4_elastic_p <= 1.0e-6
            @test r.eta_inv_max <= 1.0e-27
            # The solid does not move: measured below 1e-18 m/s.
            @test r.v_max <= 1.0e-15
        end
    end

    @testset "centre pressure below the backward-Euler value, first order (spherical_flow = $sf)" for sf in
                                                                                                      (
        false, true
    )
        e33 = runs[(33, sf, :base)].centre_err_m20
        e65 = runs[(65, sf, :base)].centre_err_m20
        # The vent cells sit inside R, so the fluid drains through a smaller sphere and the
        # centre pressure falls below the exact value.
        @test -0.1 < e33 < 0.0
        @test -0.04 < e65 < 0.0
        # Measured order 1.36 (disk) and 1.31 (sphere) from Nx = 33 to 65; the record
        # checks the order from Nx = 65 to 129.
        @test 0.7 < log2(e33 / e65) < 1.6
        @test isapprox(
            measured["centre_err_m20_$(sf ? "on" : "off")_nx65"], abs(e65); rtol=1.0e-14
        )
        e5 = runs[(65, sf, :base)].centre_err_m5
        @test -0.01 < e5 < 0.0
        # The L2 error over x <= 0.9 falls at first order.
        l33, l65 = runs[(33, sf, :base)].l2_m20, runs[(65, sf, :base)].l2_m20
        @test isapprox(log2(l33 / l65), 1.0; atol=0.35)
        @test 0.0 < l65 < 0.08
    end

    @testset "centre fluid pressure at Nx = 65, m = 20 (pinned)" begin
        p_sph = runs[(65, true, :base)].p_centre_m20
        p_disk = runs[(65, false, :base)].p_centre_m20
        @test isapprox(p_sph, 0.684822155440; rtol=1.0e-6)
        @test isapprox(p_disk, 0.827932945493; rtol=1.0e-6)
        # Scale guard: (Pf - P_amb) / P0 lies between the vent value 0 and the start value 1.
        @test 0.5 < p_sph < p_disk < 1.0
    end

    @testset "total pressure follows kbw (Pf - pf0) with the closed-form kbw (spherical_flow = $sf)" for sf in
                                                                                                         (
        false, true
    )
        for n in (33, 65)
            r = runs[(n, sf, :base)]
            @test r.pt_kbw <= 1.0e-9
            # The centre total pressure: (Pt - P_amb) / P0 = 1 - kbw (1 - p).
            kbw = SPDV.biot_willis_closed_form(SPDV.BETA_PHI, SPDV.BETA_S)
            @test isapprox(r.pt_centre_m20, 1 - kbw * (1 - r.p_centre_m20); atol=1.0e-9)
            # kbw = 1 (the Pt row without the coupling) would put Pt on Pf.
            @test !isapprox(r.pt_centre_m20, r.p_centre_m20; atol=0.05)
        end
    end

    @testset "vent cells hold P_amb and vent the stored fluid (spherical_flow = $sf)" for sf in
                                                                                          (
        false, true
    )
        for n in (33, 65)
            r = runs[(n, sf, :base)]
            # Measured 1.2e-3 to 1.9e-3 at Nx = 33 and 6.4e-4 to 1.0e-3 at Nx = 65.
            @test r.vent_dev < 3.0e-3
            @test r.vent_closure <= 1.0e-8
            @test isapprox(r.vent_flux, -r.storage_flux; rtol=1.0e-8)
        end
        @test runs[(65, sf, :base)].vent_dev < runs[(33, sf, :base)].vent_dev
    end

    @testset "sphere to disk centre ratio against the backward-Euler series" begin
        ratio = measured["centre_ratio_m20_nx65"]
        x_c = runs[(65, true, :base)].x_centre
        @test isapprox(x_c, SPDV.centre_radius(65); rtol=1.0e-12)
        @test isapprox(x_c, runs[(65, false, :base)].x_centre; rtol=1.0e-14)
        ref = SPDV.expected_metric_references()["centre_ratio_m20_nx65"]
        @test isapprox(
            ref, SPDV.exact_p(x_c, 3, 20) / SPDV.exact_p(x_c, 2, 20); rtol=1.0e-14
        )
        @test isapprox(ratio, ref; rtol=0.022)
        # Face weights of 1 in the Pf row give the disk value: a ratio near 1.
        @test !isapprox(ratio, 1.0; rtol=0.1)
        @test isapprox(
            ratio,
            runs[(65, true, :base)].p_centre_m20 / runs[(65, false, :base)].p_centre_m20;
            rtol=1.0e-14,
        )
    end

    @testset "nondimensional pressure is the same for 2 solid compressibilities (spherical_flow = $sf)" for sf in
                                                                                                            (
        false, true
    )
        # bs = 1e-11 changes S by a quarter and dt with it; tau and p must not change.
        s_inv = measured["s_invariance_$(sf ? "on" : "off")_nx33"]
        @test s_inv <= 1.0e-9
        alt, base = runs[(33, sf, :beta_alt)], runs[(33, sf, :base)]
        @test isapprox(alt.dt / base.dt, 9.7273e-11 / 1.3e-10; rtol=1.0e-4)
        @test isapprox(alt.p_centre_m20, base.p_centre_m20; atol=1.0e-9)
    end

    @testset "6-variable path equals the 4-variable path (spherical_flow = $sf)" for sf in (
        false, true
    )
        six, four = runs[(33, sf, :elim_off)], runs[(33, sf, :base)]
        @test measured["elim_$(sf ? "on" : "off")_nx33"] <= 1.0e-10
        @test maximum(abs, six.pf .- four.pf) / SPDV.P0 <= 1.0e-10
        @test maximum(abs, six.qxD .- four.qxD) / maximum(abs, four.qxD) <= 1.0e-10
    end

    record_path = joinpath(
        @__DIR__, "..", "benchmarks", "records", "spherical_pressure_diffusion.toml"
    )
    raw_record_all = TOML.parsefile(record_path)
    recorded = Dict(
        String(m["id"]) => Float64(m["value"]) for m in raw_record_all["metrics"]
    )
    recorded_only = filter(kv -> !haskey(measured, kv.first), recorded)
    bench = (; runs, metrics=merge(measured, recorded_only))

    @testset "metrics: 27 from the runs, 13 recorded values for Nx = 129 and the orders" begin
        @test length(measured) == 27
        @test length(recorded_only) == 13
        @test all(endswith(k, "129") for k in keys(recorded_only))
        @test Set(keys(bench.metrics)) == Set(SPDV.EXPECTED_METRIC_IDS)
    end

    @testset "vent deviation covers every surface rock cell and is infinite for a closed vent" begin
        coords = diffusion_coords(41)
        n = coords.Ny1
        p_in = SPDV.P_AMB + SPDV.P0
        cells = SPDV._rock_cells(coords)
        surface = [(i, j) for (i, j, _) in cells if SPDV._surface_rock(coords, i, j)]
        @test length(surface) == 2 * (41 - 1)
        pf = fill(SPDV.P_AMB, n, n)
        s_vent = zeros(n, n)
        for (i, j, _) in cells
            pf[i, j] = p_in
        end
        for (i, j) in surface
            pf[i, j] = SPDV.P_AMB
            s_vent[i, j] = 1.0e-12
        end
        kbw = SPDV.biot_willis_closed_form(SPDV.BETA_PHI, SPDV.BETA_S)
        grids(pf, pr, s) = (; pf, pr, pf0=copy(pf), S_vent_grid=s)
        pr = SPDV.P_AMB .+ SPDV.P0 .+ kbw .* (pf .- p_in)
        dt = SPDV.time_step_s(SPDV.BETA_S)
        m = SPDV._final_metrics(grids(pf, pr, s_vent), coords, false, dt, SPDV.BETA_S)
        @test isapprox(m.vent_dev, 0.0; atol=1.0e-15)
        @test isapprox(m.pt_kbw, 0.0; atol=1.0e-12)
        @test m.n_rock == length(cells)
        # A surface cell 0.3 P0 above P_amb sets the deviation; a closed one makes it Inf.
        i0, j0 = first(surface)
        pf_high = copy(pf)
        pf_high[i0, j0] = SPDV.P_AMB + 0.3 * SPDV.P0
        m_high = SPDV._final_metrics(
            grids(pf_high, pr, s_vent), coords, false, dt, SPDV.BETA_S
        )
        @test isapprox(m_high.vent_dev, 0.3; rtol=1.0e-12)
        # Pf raised by 0.3 P0 with Pt kept gives the total pressure check 0.3 kbw.
        @test isapprox(m_high.pt_kbw, 0.3 * kbw; rtol=1.0e-12)
        s_closed = copy(s_vent)
        s_closed[i0, j0] = 0.0
        m_closed = SPDV._final_metrics(
            grids(pf, pr, s_closed), coords, false, dt, SPDV.BETA_S
        )
        @test isequal(m_closed.vent_dev, Inf)
        # Cells away from the sphere: no rock, no error measure.
        far = diffusion_coords(41; shift=10 * SPDV.R_PLANET)
        @test isempty(SPDV._rock_cells(far))
        @test_throws ArgumentError SPDV._l2_inner((; pf), far, 3, 20)
        @test_throws ArgumentError SPDV.centre_radius(SPDV.NX_MIN - 1)
    end

    @testset "metrics_from_runs with every grid of NX_LIST" begin
        fake(e) = (;
            centre_err_m5=(-e / 10),
            centre_err_m20=(-e),
            l2_m20=2e,
            pt_kbw=1.0e-11,
            p_centre_m20=1.0 - e,
            vent_dev=1.0e-4,
            vent_closure=1.0e-12,
            v_max=1.0e-20,
            pf=fill(1.0, 2, 2),
            pr=fill(1.0, 2, 2),
            qxD=fill(1.0, 2, 2),
            qyD=fill(1.0, 2, 2),
            xp=[SPDV.X_CENTRE, SPDV.X_CENTRE + 1.0],
            yp=[SPDV.Y_CENTRE, SPDV.Y_CENTRE + 1.0],
        )
        fake_runs = Dict{Tuple,Any}()
        for (n, e) in ((33, 0.04), (65, 0.02), (129, 0.01)), sf in (false, true)
            fake_runs[(n, sf, :base)] = fake(sf ? e : 1.5e)
        end
        for sf in (false, true), v in (:elim_off, :beta_alt)
            fake_runs[(33, sf, v)] = fake(sf ? 0.04 : 0.06)
        end
        all_metrics = SPDV.metrics_from_runs(fake_runs)
        @test Set(keys(all_metrics)) == Set(SPDV.EXPECTED_METRIC_IDS)
        # The error halves from Nx = 65 to 129, so every order is 1.
        for tag in ("on", "off"), name in ("centre", "l2")
            @test isapprox(
                all_metrics["$(name)_order_m20_$(tag)_65_129"], 1.0; rtol=1.0e-12
            )
        end
        @test isapprox(all_metrics["centre_ratio_m20_nx129"], 0.99 / 0.985; rtol=1.0e-14)
        @test isapprox(all_metrics["centre_err_m20_on_nx65"], 0.02; rtol=1.0e-14)
        @test isapprox(all_metrics["elim_on_nx33"], 0.0; atol=1.0e-15)
        delete!(fake_runs, (129, true, :base))
        delete!(fake_runs, (33, true, :beta_alt))
        partial = SPDV.metrics_from_runs(fake_runs)
        @test !haskey(partial, "centre_order_m20_on_65_129")
        @test !haskey(partial, "centre_ratio_m20_nx129")
        @test !haskey(partial, "s_invariance_on_nx33")
        @test haskey(partial, "centre_order_m20_off_65_129")
    end

    @testset "Validation record schema and drift guards" begin
        @test isfile(record_path)
        raw_record = TOML.parsefile(record_path)
        valid_ok, valid_errors = SPDV.validate_record_dict(
            raw_record; require_doc_pages=true
        )
        @test valid_ok
        @test isempty(valid_errors)
        rec_ids = [String(m["id"]) for m in raw_record["metrics"]]
        @test rec_ids == collect(SPDV.EXPECTED_METRIC_IDS)
        @test Set(keys(SPDV.PHYSICS_REGRESSION_BANDS)) == Set(SPDV.EXPECTED_METRIC_IDS)
        missing_key = deepcopy(raw_record)
        delete!(missing_key["config"], "storage_per_Pa")
        ok_missing, errs_missing = SPDV.validate_record_dict(
            missing_key; require_doc_pages=false, check_doc_sync=false
        )
        @test !ok_missing
        @test any(occursin("storage_per_Pa", e) for e in errs_missing)

        for (mid, value) in bench.metrics
            band = SPDV.PHYSICS_REGRESSION_BANDS[mid]
            @test isapprox(value, recorded[mid]; atol=band.atol, rtol=band.rtol)
        end
        ok_full, errs_full = SPDV.check_physics_regression_bands(bench.metrics, record_path)
        @test ok_full
        @test isempty(errs_full)
        # A storage coefficient 2 % too large moves the centre error at Nx = 65 by about
        # 0.016 for the sphere (measured with RX, RY scaled by 1.02): outside the band.
        target = "centre_err_m20_on_nx65"
        moved = copy(bench.metrics)
        moved[target] += 2 * SPDV.PHYSICS_REGRESSION_BANDS[target].atol
        ok_moved, errs_moved = SPDV.check_physics_regression_bands(moved, record_path)
        @test !ok_moved
        @test occursin(target, only(errs_moved))
        missing_metric = copy(bench.metrics)
        delete!(missing_metric, target)
        ok_gone, errs_gone = SPDV.check_physics_regression_bands(
            missing_metric, record_path
        )
        @test !ok_gone
        @test any(occursin("Measured metrics missing '$target'", e) for e in errs_gone)
        no_metric = deepcopy(raw_record)
        filter!(m -> m["id"] != target, no_metric["metrics"])
        mktempdir() do tmpdir
            path = write_diffusion_toml(joinpath(tmpdir, "short.toml"), no_metric)
            ok_short, errs_short = SPDV.check_physics_regression_bands(bench.metrics, path)
            @test !ok_short
            @test any(occursin("record missing metric '$target'", e) for e in errs_short)
        end
        @test_throws DomainError SPDV.check_physics_regression_bands(bench.metrics, "")
        @test !first(SPDV.check_physics_regression_bands(bench.metrics, record_path * ".x"))
    end

    @testset "Record matches the configuration and the reference values" begin
        raw_record = TOML.parsefile(record_path)
        cfg = SPDV.make_config(SPDV.NX_REFERENCE; spherical_flow=true)
        sha = SPDV.canonical_config_sha256(cfg)
        @test sha == raw_record["metadata"]["config_sha256"]
        @test length(sha) == 64
        @test SPDV.canonical_config_sha256(
            SPDV.make_config(SPDV.NX_REFERENCE; spherical_flow=false)
        ) != sha
        @test SPDV.canonical_config_sha256(
            SPDV.make_config(SPDV.NX_REFERENCE; spherical_flow=true, beta_s=SPDV.BETA_S_ALT)
        ) != sha
        text = SPDV.canonical_config_text(cfg)
        for name in SPDV.HASHED_CASE_CONSTANTS
            @test occursin("case.$(name)=", text)
        end
        @test length(SPDV.HASHED_CASE_CONSTANTS) == 32
        @test bytes2hex(SHA.sha256(text)) == sha
        summary = SPDV.expected_config_summary(cfg)
        rec_cfg = raw_record["config"]
        @test Set(keys(summary)) == Set(keys(rec_cfg))
        for (key, expected) in summary
            if expected isa AbstractFloat
                @test isapprox(Float64(rec_cfg[key]), expected; rtol=1.0e-12, atol=0.0)
            else
                @test rec_cfg[key] == expected
            end
        end
        @test isapprox(summary["storage_per_Pa"], 1.3e-10; rtol=1.0e-14)
        @test isapprox(summary["kbw"], 0.6; rtol=1.0e-14)
        refs = SPDV.expected_metric_references()
        @test Set(keys(refs)) == Set(SPDV.EXPECTED_METRIC_IDS)
        # The reference centre ratio depends on the centre radius of the grid.
        @test isapprox(refs["centre_ratio_m20_nx65"], 0.8401398; rtol=1.0e-6)
        @test isapprox(refs["centre_ratio_m20_nx129"], 0.8403029; rtol=1.0e-6)
        for m in raw_record["metrics"]
            @test isapprox(
                Float64(m["reference"]), refs[String(m["id"])]; rtol=1.0e-12, atol=1.0e-12
            )
        end
    end

    @testset "Reproducibility check, CLI table and export round trip" begin
        raw_record = TOML.parsefile(record_path)
        ok_repro, errs_repro = SPDV.check_validation_reproducibility(bench, record_path)
        @test ok_repro
        @test isempty(errs_repro)
        @test run_diffusion_validation_check(record_path; bench=bench)
        @test run_validation_cli(["--check", record_path]; diffusion_bench=bench) == 0
        @test !run_diffusion_validation_check(record_path * ".missing"; bench=bench)
        @test_throws DomainError run_diffusion_validation_check(""; bench=bench)
        @test validation_case("spherical_pressure_diffusion").bench_keyword ==
            :diffusion_bench

        mktempdir() do tmpdir
            empty_bench = (; runs=Dict{Tuple,Any}(), metrics=Dict{String,Float64}())
            no_arg = capture_diffusion_stdio(
                () -> run_validation_cli(
                    ["--check"];
                    bench=empty_bench,
                    darcy_bench=empty_bench,
                    compaction_bench=empty_bench,
                    diffusion_bench=bench,
                ),
            )
            @test no_arg.value == 1
            @test occursin("spherical_pressure_diffusion.toml", no_arg.out)
            @test occursin("[reproducibility]", no_arg.err)
            @test last(VALIDATION_CASES).id == "spherical_pressure_diffusion"

            exported = joinpath(tmpdir, "exported.toml")
            commit = String(raw_record["metadata"]["src_commit"])
            @test SPDE.run_export_cli(
                [exported]; src_commit=commit, bench=bench, check_clean_tree=false
            ) == 0
            raw_exported = TOML.parsefile(exported)
            ok_valid, errs_valid = SPDV.validate_record_dict(
                raw_exported; require_doc_pages=false, check_doc_sync=false
            )
            @test ok_valid
            @test isempty(errs_valid)
            @test run_diffusion_validation_check(exported; bench=bench)
            @test [String(m["id"]) for m in raw_exported["metrics"]] == collect(SPDV.EXPECTED_METRIC_IDS)

            # A moved value is refused without the flag and written with it.
            target = "centre_err_m20_on_nx65"
            moved_metrics = copy(bench.metrics)
            moved_metrics[target] += 2 * SPDV.PHYSICS_REGRESSION_BANDS[target].atol
            moved_bench = (; runs=bench.runs, metrics=moved_metrics)
            @test SPDE.run_export_cli(
                [exported]; src_commit=commit, bench=moved_bench, check_clean_tree=false
            ) == 1
            @test SPDE.run_export_cli(
                ["--accept-physics-change", exported];
                src_commit=commit,
                bench=moved_bench,
                check_clean_tree=false,
            ) == 0
            @test SPDE.run_export_cli(
                ["--bad-flag"]; src_commit=commit, bench=bench, check_clean_tree=false
            ) == 2
            @test SPDE.run_export_cli(
                [exported];
                src_commit=commit,
                bench=bench,
                check_clean_tree=true,
                dirty_status=" M src/config.jl",
            ) == 1
            @test_throws DomainError SPDE.run_export_cli(
                [exported]; src_commit="abc", bench=bench, check_clean_tree=false
            )
            # A value that fails its own tolerance is never written.
            failing = copy(bench.metrics)
            failing["vent_closure_on_nx65"] = 1.0
            @test SPDE.run_export_cli(
                ["--accept-physics-change", exported];
                src_commit=commit,
                bench=(; runs=bench.runs, metrics=failing),
                check_clean_tree=false,
            ) == 1
            cases = (
                (0.01, 0.0, "le_atol"), (0.01, 1.0, "le_atol"), (1.0, 1.001, "isapprox")
            )
            @test isequal(
                [SPDE.metric_status(v, ref, 0.05, 0.0, c) for (v, ref, c) in cases],
                ["PASS", "FAIL", "PASS"],
            )
            @test isequal(SPDE.metric_status(1.0, 1.0, 0.0, 0.0, "other"), "FAIL")
            @test occursin(r"^[0-9a-f]{40}$", SPDE.last_src_commit())
            @test first(SPDE.check_clean_source_tree()) isa Bool

            # A changed configuration value and hash are rejected.
            changed_cfg = deepcopy(raw_record)
            changed_cfg["config"]["k_m2"] *= 2
            changed_cfg["metadata"]["config_sha256"] = repeat("0", 64)
            cfg_path = write_diffusion_toml(joinpath(tmpdir, "changed.toml"), changed_cfg)
            ok_cfg, errs_cfg = SPDV.check_validation_reproducibility(bench, cfg_path)
            @test !ok_cfg
            @test any(occursin("config.k_m2", e) for e in errs_cfg)
            @test any(occursin("config_sha256", e) for e in errs_cfg)
            # Other metric ids and a moved reference; a missing config key on its own.
            no_kbw = deepcopy(raw_record)
            delete!(no_kbw["config"], "kbw")
            @test SPDV._config_mismatches(
                no_kbw, SPDV.make_config(65; spherical_flow=true)
            ) == ["Recorded config missing key 'kbw'"]
            odd = deepcopy(raw_record)
            odd["metrics"][23]["reference"] = 0.5
            odd["metrics"][23]["value"] = 0.5
            odd["metrics"][end]["id"] = "not_a_metric"
            odd_path = write_diffusion_toml(joinpath(tmpdir, "odd.toml"), odd)
            ok_odd, errs_odd = SPDV.check_validation_reproducibility(bench, odd_path)
            @test !ok_odd
            @test any(occursin("Recorded metric IDs", e) for e in errs_odd)
            @test any(occursin("does not match recomputed reference", e) for e in errs_odd)
            @test any(occursin("Benchmark output missing metric", e) for e in errs_odd)
            bad_comp = deepcopy(raw_record)
            bad_comp["metrics"][1]["comparison"] = "other"
            bad_comp_path = write_diffusion_toml(joinpath(tmpdir, "comp.toml"), bad_comp)
            ok_comp, errs_comp = SPDV.check_validation_reproducibility(bench, bad_comp_path)
            @test !ok_comp
            @test any(occursin("unsupported comparison", e) for e in errs_comp)
            failing_bench = (; runs=bench.runs, metrics=failing)
            ok_fail, errs_fail = SPDV.check_validation_reproducibility(
                failing_bench, record_path
            )
            @test !ok_fail
            @test any(occursin("fails tolerance", e) for e in errs_fail)
            @test any(occursin("differs from recorded value", e) for e in errs_fail)
            # A record of another Julia minor version skips the replay of recorded values.
            other_julia = deepcopy(raw_record)
            other_julia["metadata"]["julia_version"] = "1.0.0"
            other_path = write_diffusion_toml(joinpath(tmpdir, "julia.toml"), other_julia)
            ok_other, _ = SPDV.check_validation_reproducibility(bench, other_path)
            @test ok_other
            ok_norep, _ = SPDV.check_validation_reproducibility(
                bench, record_path; check_recorded_values=false
            )
            @test ok_norep
        end
    end

    @testset "checks stop on an interrupt and report other errors with a stack trace" begin
        commit = String(raw_record_all["metadata"]["src_commit"])
        @test_throws InterruptException run_diffusion_validation_check(
            record_path; bench=DiffusionInterruptingBench()
        )
        crash = capture_diffusion_stdio(
            () ->
                run_diffusion_validation_check(record_path; bench=DiffusionCrashingBench()),
        )
        @test crash.value === false
        @test occursin("solver failure inside the benchmark", crash.err)
        @test occursin("Stacktrace", crash.err)
        mktempdir() do tmpdir
            target = joinpath(tmpdir, "record.toml")
            @test_throws InterruptException SPDE.run_export_cli(
                [target];
                src_commit=commit,
                bench=DiffusionInterruptingBench(),
                check_clean_tree=false,
            )
            failed = capture_diffusion_stdio(
                () -> SPDE.run_export_cli(
                    [target];
                    src_commit=commit,
                    bench=DiffusionCrashingBench(),
                    check_clean_tree=false,
                ),
            )
            @test failed.value == 1
            @test occursin("Stacktrace", failed.err)
            @test !isfile(target)
            # A clean-tree check whose `git status` fails is reported and gives exit code 1.
            not_repo = mkpath(joinpath(tmpdir, "not_a_repo"))
            withenv("GIT_DIR" => joinpath(not_repo, "absent.git")) do
                for cli in (SPDE.run_export_cli, SPDE.export_main)
                    no_git = capture_diffusion_stdio(
                        () -> cli(
                            [target];
                            src_commit=commit,
                            bench=bench,
                            check_clean_tree=true,
                            repo_root=not_repo,
                        ),
                    )
                    @test no_git.value == 1
                    @test occursin("status --porcelain", no_git.err)
                    @test !isfile(target)
                end
            end
        end
    end

    @testset "export command line: one record path, unknown flags and an empty path" begin
        commit = String(raw_record_all["metadata"]["src_commit"])
        mktempdir() do tmpdir
            cd(tmpdir) do
                code = SPDE.run_export_cli(
                    ["-h"]; src_commit=commit, bench=bench, check_clean_tree=false
                )
                @test code == 2
                @test !isfile(joinpath(tmpdir, "-h"))
            end
            p1, p2 = joinpath(tmpdir, "first.toml"), joinpath(tmpdir, "second.toml")
            @test SPDE.run_export_cli(
                [p1, p2]; src_commit=commit, bench=bench, check_clean_tree=false
            ) == 2
            @test !isfile(p1) && !isfile(p2)
            empty_run = capture_diffusion_stdio(
                () -> SPDE.export_main(
                    [""]; src_commit=commit, bench=bench, check_clean_tree=false
                ),
            )
            @test empty_run.value == 1
            @test occursin("Record path must be non-empty", empty_run.err)
            @test SPDE.export_main(
                [p2]; src_commit=commit, bench=bench, check_clean_tree=false
            ) == 0
            @test isfile(p2)
            # The export tool reads the validation module of the module that includes it.
            @test SPDE.diffusion_validation_module() === SPDV
            @test_throws ErrorException SPDE.diffusion_validation_module((Module(:Empty),))
        end
    end

    @testset "malformed records are reported, not thrown" begin
        mktempdir() do tmpdir
            broken = joinpath(tmpdir, "broken.toml")
            write(broken, "[metadata\nid = ")
            parsed = capture_diffusion_stdio(
                () -> run_diffusion_validation_check(broken; bench=bench)
            )
            @test parsed.value === false
            @test occursin("[schema] Cannot parse", parsed.err)
            for check in (
                () -> SPDV.check_validation_reproducibility(bench, broken),
                () -> SPDV.check_physics_regression_bands(bench.metrics, broken),
            )
                ok_parse, errs_parse = check()
                @test !ok_parse
                @test occursin("Cannot parse", only(errs_parse))
            end
            @test_throws DomainError SPDV.check_validation_reproducibility(bench, "")
            @test !first(SPDV.check_validation_reproducibility(bench, broken * ".x"))
            schema_bad = deepcopy(raw_record_all)
            delete!(schema_bad, "config")
            schema_path = write_diffusion_toml(joinpath(tmpdir, "schema.toml"), schema_bad)
            ok_schema, errs_schema = SPDV.check_validation_reproducibility(
                bench, schema_path
            )
            @test !ok_schema
            @test any(occursin("config", e) for e in errs_schema)
        end
    end

    @testset "figure data export: CSV columns and the exported files" begin
        script_host = Module(:DiffusionFigureData)
        Core.eval(script_host, :(include(path) = Base.include($script_host, path)))
        Base.include(
            script_host,
            joinpath(
                @__DIR__, "..", "scripts", "export_validation_pressure_diffusion_data.jl"
            ),
        )
        write_columns = getfield(script_host, :write_columns)
        export_all = getfield(script_host, :export_all)
        mktempdir() do tmpdir
            path = joinpath(tmpdir, "columns.csv")
            write(path, "kept\n")
            for (header, columns) in (
                (["a"], Vector{Float64}[]),
                (["a", "b", "c"], [[1.0, 2.0], [3.0, 4.0]]),
                (["a", "b"], [[1.0, 2.0], [3.0, 4.0, 5.0]]),
            )
                @test_throws ArgumentError Base.invokelatest(
                    write_columns, path, header, columns
                )
            end
            @test read(path, String) == "kept\n"
            outdir = joinpath(tmpdir, "data")
            @test Base.invokelatest(export_all, outdir; runs=runs) == outdir
            files = sort(readdir(outdir))
            @test files == sort([
                "convergence.csv",
                "profile_off_65.csv",
                "profile_on_65.csv",
                "reference_curves.csv",
            ])
            conv = readlines(joinpath(outdir, "convergence.csv"))
            @test conv[1] == "nx,l2_off,l2_on,centre_off,centre_on"
            @test length(conv) == 1 + length(SPDV.NX_LIST)
            row65 = parse.(Float64, split(conv[3], ","))
            @test isapprox(row65[2], recorded["l2_m20_off_nx65"]; rtol=1.0e-15)
            @test isapprox(row65[5], recorded["centre_err_m20_on_nx65"]; rtol=1.0e-15)
            curves = readlines(joinpath(outdir, "reference_curves.csv"))
            @test curves[1] == "x,p_be_m5_off,p_be_m5_on,p_be_m20_off,p_be_m20_on"
            # The exact curves start at the centre value and end at 0 on the surface.
            first_curve = parse.(Float64, split(curves[2], ","))
            last_curve = parse.(Float64, split(curves[end], ","))
            @test isapprox(first_curve[5], SPDV.exact_p(0.0, 3, 20); rtol=1.0e-15)
            @test isapprox(last_curve[1], 1.0; rtol=1.0e-15)
            @test all(isapprox(v, 0.0; atol=1.0e-12) for v in last_curve[2:end])
            prof = readlines(joinpath(outdir, "profile_on_65.csv"))
            @test prof[1] == "x,p_m5,p_m20"
            vals = [parse.(Float64, split(l, ",")) for l in prof[2:end]]
            @test length(vals) == runs[(65, true, :base)].n_rock
            # The cell values stay within 0.1 of the series at m = 20.
            @test all(abs(v[3] - SPDV.exact_p(v[1], 3, 20)) < 0.1 for v in vals)
        end
    end
end
