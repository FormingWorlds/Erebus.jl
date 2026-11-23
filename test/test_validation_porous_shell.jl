using SHA
using Test
using TOML
using Erebus

isdefined(@__MODULE__, :PorousShellValidation) ||
    include(joinpath(@__DIR__, "..", "benchmarks", "validation_porous_shell.jl"))
include(joinpath(@__DIR__, "..", "tools", "check_validation.jl"))
isdefined(@__MODULE__, :PorousShellExport) ||
    include(joinpath(@__DIR__, "..", "tools", "export_validation_porous_shell.jl"))
isdefined(@__MODULE__, :PorousShellSuiteExport) ||
    include(joinpath(@__DIR__, "..", "tools", "export_validation_porous_shell_suite.jl"))
const PSV = PorousShellValidation

"""
    PorousShellInterruptingBench

Benchmark stand-in whose every property access throws `InterruptException`, as a Ctrl-C
during the benchmark does.
"""
struct PorousShellInterruptingBench end
Base.getproperty(::PorousShellInterruptingBench, ::Symbol) = throw(InterruptException())

"""
    PorousShellCrashingBench

Benchmark stand-in whose every property access throws an `ErrorException`, as a solver
failure during the benchmark does.
"""
struct PorousShellCrashingBench end
function Base.getproperty(::PorousShellCrashingBench, ::Symbol)
    return error("solver failure inside the benchmark")
end

"""
    write_porous_shell_toml(path, dict)

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
function write_porous_shell_toml(path::String, dict::AbstractDict)
    open(path, "w") do io
        return TOML.print(io, dict)
    end
    return path
end

"""
    capture_porous_shell_stdio(f)

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
function capture_porous_shell_stdio(f)
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
    write_suite_folder(dir, raw; tamper=Dict())

Write a folder of suite result files, as `scripts/run_validation_porous_shell_suite.jl`
writes them, with the depths and temperature errors of the suite record `raw`.

Parameters
----------
dir : String
    Output folder.
raw : AbstractDict
    Parsed suite record.
tamper : AbstractDict, optional
    Map from a tag to a `(column, text)` pair that replaces one field of its row.

Returns
-------
String
    `dir`.
"""
function write_suite_folder(dir::String, raw::AbstractDict; tamper=Dict())
    PSVm = PorousShellValidation
    value(id) = only(m["value"] for m in raw["metrics"] if m["id"] == id)
    reference(id) = only(m["reference"] for m in raw["metrics"] if m["id"] == id)
    write(joinpath(dir, "header.csv"), join(PSVm.SUITE_RESULT_HEADER, ",") * "\n")
    for (tag, nx, tf, phi, n, k_air) in PSVm.SUITE_CASES
        d, dref = value("shell_depth_$tag"), reference("shell_depth_$tag")
        row = Any[
            tag,
            nx,
            tf,
            phi,
            n,
            something(k_air, -1.0),
            48_000.0 / (nx - 1),
            d,
            d,
            d,
            dref,
            0.0,
            value("tmax_rel_l2_$tag"),
            0.0,
            800.0,
            800.0,
            1.0e-10,
            1.0,
        ]
        if haskey(tamper, tag)
            col, text = tamper[tag]
            row[findfirst(==(col), PSVm.SUITE_RESULT_HEADER)] = text
        end
        write(joinpath(dir, "result_$tag.csv"), join(row, ",") * "\n")
    end
    return dir
end

"""
    energy_test_grids(rhocp, dT, q; q_metric=true)

Return 4 x 4 grid fields whose 2 x 2 interior has the heat capacity `rhocp` [J/(m^3 K)],
the temperature change `dT` [K] from 300 K and the heat source `q` [W/m^3], split over
the 8 source terms of the thermal right-hand side. The boundary nodes hold values that
`thermal_step_energy_residual` must ignore.

Parameters
----------
rhocp, dT, q : Matrix{Float64}
    2 x 2 interior values.
q_metric : Bool, optional
    Whether `Q_metric` is an array (`true`) or `nothing` (default `true`).

Returns
-------
NamedTuple
    Fields `tk1`, `tk2`, `RHOCP`, `HR`, `HA`, `HS`, `DHP`, `Q_metric`, `Q_lat_grid`,
    `Q_seg_grid`, `Q_fluid_grid`.
"""
function energy_test_grids(rhocp, dT, q; q_metric::Bool=true)
    field(inner, edge) = (a=fill(edge, 4, 4); a[2:3, 2:3].=inner; a)
    tk1 = fill(300.0, 4, 4)
    names = q_metric ? 8 : 7
    part = q ./ names
    return (;
        tk1,
        tk2=field(300.0 .+ dT, 1.0e4),
        RHOCP=field(rhocp, 3.0e10),
        HR=field(part, 1.0e9),
        HA=field(part, 1.0e9),
        HS=field(part, 1.0e9),
        DHP=field(part, 1.0e9),
        Q_metric=q_metric ? field(part, 1.0e9) : nothing,
        Q_lat_grid=field(part, 1.0e9),
        Q_seg_grid=field(part, 1.0e9),
        Q_fluid_grid=field(part, 1.0e9),
    )
end

@testset "Porous Shell Validation" begin
    @testset "Reference laws agree with the Erebus porosity code" begin
        for P in (0.0, 0.05, 0.13, 1.3, 40.0)
            @test isapprox(
                PSV.phi_cold(P), Erebus.porosity_phi0(P, 1.0; phi_cap=0.99); rtol=1e-14
            )
        end
        for phi in (0.0, 0.1, prevfloat(0.2), 0.2, 0.3, 0.4, nextfloat(0.4), 0.6, 0.85)
            @test isapprox(
                PSV.k_porous_ref(3.0, phi), Erebus.k_porous(3.0, phi); rtol=1e-14
            )
        end
        # States over both Eq. 43 branches, near onset and at the contact-overlap limit
        for (phi, phi0, P, T) in (
            (0.40, 0.42, 5.0, 900.0),
            (0.30, 0.43, 0.5, 1100.0),
            (0.08, 0.42, 20.0, 800.0),
            (0.4199, 0.42, 2.0, 1000.0),
            (0.05, 0.75, 3.0, 950.0),
            (0.6, 0.75, 1.0, 1200.0),
        )
            r_ref = PSV.sinter_rate_ref(phi, phi0, P, T)
            r_src = Erebus.porosity_decrease_rate(phi, phi0, P, T)
            @test 0.0 < r_ref < Inf
            @test isapprox(r_ref, r_src; rtol=1e-9)
        end
        # Below the sintering temperature and at P = 0 the rate is 0
        @test iszero(PSV.sinter_rate_ref(0.40, 0.42, 5.0, 699.0))
        @test iszero(PSV.sinter_rate_ref(0.40, 0.42, 0.0, 900.0))
    end

    @testset "Sintering step converges to a fine RK4 integration" begin
        phi, phi0, P, T = 0.40, 0.42, 5.0, 900.0
        dt = 2.0e4 * PSV.YEARLENGTH_S
        rk = phi
        nrk = 20_000
        h = dt / nrk
        f(x) = -PSV.sinter_rate_ref(x, phi0, P, T)
        for _ in 1:nrk
            k1 = f(rk)
            k2 = f(rk + 0.5h * k1)
            k3 = f(rk + 0.5h * k2)
            k4 = f(rk + h * k3)
            rk += h / 6 * (k1 + 2k2 + 2k3 + k4)
        end
        step = PSV.sinter_step_ref(phi, phi0, P, T, dt)
        @test 0.0 < rk < phi
        @test isapprox(step, rk; atol=2.0e-6)
        # A step at 650 K changes nothing
        @test PSV.sinter_step_ref(phi, phi0, P, 650.0, dt) === phi
    end

    @testset "Hydrostatic pressure of a uniform sphere" begin
        R, Nr = 20_000.0, 2000
        P = PSV.hydrostatic_pressure(fill(0.1, Nr), R)
        rho = PSV.RHO_SOL * 0.9
        r = [(i - 0.5) * R / Nr for i in 1:Nr]
        exact = (2pi / 3) * PSV.G_GRAV * rho^2 .* (R^2 .- r .^ 2)
        @test isapprox(P[1], exact[1]; rtol=1e-6)
        @test maximum(abs.(P .- exact)) / exact[1] < 1e-6
    end

    @testset "Cold-pressed start and an adiabatic porous interior" begin
        R, dt_yr, nsteps = 20_000.0, 1.0e4, 20
        res = PSV.solve_porous_sphere(;
            R=R,
            phi_init=0.75,
            t_form_yr=1.5e6,
            t_end_yr=nsteps * dt_yr,
            dt_yr=dt_yr,
            Nr=400,
        )
        # Below 700 K the porosity stays on the cold-pressed profile of the start
        P = PSV.hydrostatic_pressure(res.phi, R)
        expect = [min(0.75, PSV.phi_cold(p * 1.0e-5), PSV.PHI_CAP) for p in P]
        @test maximum(abs.(res.phi .- expect)) < 1e-12
        @test 0.42 < res.phi[1] < 0.43 && isapprox(res.phi[end], 0.75; atol=1e-12)
        # Heat and heat capacity both scale with 1 - phi: the centre heats as solid rock
        dt = dt_yr * PSV.YEARLENGTH_S
        t0 = 1.5e6 * PSV.YEARLENGTH_S
        rise =
            sum(
                PSV.RHO_BULK * PSV.q_al_step(t0 + (k - 1) * dt, dt, PSV.RATIO_AL_L16) * dt
                for k in 1:nsteps
            ) / PSV.RHOCP_SOL
        @test 50.0 < rise < 600.0
        @test isapprox(res.T[1] - PSV.T_SPACE, rise; rtol=1e-9)
    end

    @testset "Shell depth of the 1-D profile" begin
        mask = [true, true, true, true, true, false, true, false, false, false]
        @test isapprox(PSV._shell_depth_1d(mask, 10.0, 1.0), 3.0; atol=1e-12)
        @test isapprox(PSV._shell_depth_1d(falses(10), 10.0, 1.0), 10.0; atol=1e-12)
        @test isapprox(PSV._shell_depth_1d(trues(10), 10.0, 1.0), 0.0; atol=1e-12)
    end

    @testset "Without porosity the run follows the eigenmode series" begin
        R, dt_yr, nsteps = 20_000.0, 1.0e4, 100
        res = PSV.solve_porous_sphere(;
            R=R, phi_init=0.0, t_form_yr=2.0e6, t_end_yr=nsteps * dt_yr, dt_yr=dt_yr, Nr=400
        )
        dt = dt_yr * PSV.YEARLENGTH_S
        t0 = 2.0e6 * PSV.YEARLENGTH_S
        H = [
            PSV.RHO_BULK * PSV.q_al_step(t0 + (k - 1) * dt, dt, PSV.RATIO_AL_L16) for
            k in 1:nsteps
        ]
        kw = (;
            R=R,
            T_surf=PSV.T_SPACE,
            kappa=PSV.K_SOL / PSV.RHOCP_SOL,
            rhocp=PSV.RHOCP_SOL,
            H=H,
            dt=dt,
        )
        rise = PSV.eigenmode_T(0.0; kw...) - PSV.T_SPACE
        @test res.n_steps == nsteps
        @test 100.0 < rise < 2000.0
        for i in (1, 100, 200, 300, 390)
            @test isapprox(res.T[i], PSV.eigenmode_T(res.r[i]; kw...); atol=1e-4 * rise)
        end
        @test all(iszero, res.phi) && all(iszero, res.dphi_sinter)
    end

    @testset "Case configuration" begin
        kw = (; Nx=33, phi_init=0.25, t_form_yr=2.5e6, n_steps=3, dt_yr=1.0e4)
        cfg = PSV.make_porous_shell_config(; kw...)
        @test cfg.porosity.model === :henke2012
        @test isapprox(cfg.porosity.phi_init, 0.25; rtol=1e-14)
        @test isapprox(cfg.thermodynamics.ratio_al, 5.85e-5; rtol=1e-14)
        @test isapprox(cfg.time.start_time, 2.5e6; rtol=1e-14)
        @test isapprox(cfg.time.dt_initial, cfg.time.dt_longest; rtol=1e-14)
        # Matched air: the conductivity of the cold-pressed surface rock (phi 0.25)
        @test isapprox(cfg.materials.ksolidm[3], Erebus.k_porous(3.0, 0.25); rtol=1e-12)
        @test isapprox(
            PSV.make_porous_shell_config(; kw..., k_air=3000.0).materials.ksolidm[3],
            3000.0;
            rtol=1e-14,
        )
        # phi_init above the cold-pressed surface value 0.85: air as rock at 0.85
        hi = PSV.make_porous_shell_config(; kw..., phi_init=0.9)
        @test isapprox(hi.materials.ksolidm[3], Erebus.k_porous(3.0, 0.85); rtol=1e-12)
        # The reference constants are those of the configuration and of the code
        @test isapprox(PSV.T_HALF_AL_S, cfg.thermodynamics.t_half_al; rtol=1e-14)
        @test isapprox(PSV.F_AL, cfg.thermodynamics.f_al; rtol=1e-14)
        @test isapprox(PSV.E_AL, cfg.thermodynamics.E_al; rtol=1e-14)
        @test isapprox(PSV.RHO_SOL, cfg.materials.rhosolidm[1]; rtol=1e-14)
        @test isapprox(PSV.RHO_MET, cfg.coreformation.rho_metal; rtol=1e-14)
        @test isapprox(PSV.X_FE, Erebus.X_FE_REF_CHONDRITE; rtol=1e-14)
        @test isapprox(PSV.G_GRAV, Erebus.G_GRAV; rtol=1e-14)
        @test isapprox(PSV.T_SINTER, cfg.porosity.T_sinter; rtol=1e-14)
        @test isapprox(PSV.PHI_CAP, cfg.porosity.phi_cap; rtol=1e-14)
        @test_throws DomainError PSV.make_porous_shell_config(; kw..., Nx=4)
        @test_throws DomainError PSV.make_porous_shell_config(; kw..., n_steps=0)
        @test_throws DomainError PSV.make_porous_shell_config(; kw..., dt_yr=0.0)
        @test_throws DomainError PSV.make_porous_shell_config(; kw..., dt_yr=NaN)
        @test_throws DomainError PSV.make_porous_shell_config(; kw..., t_form_yr=-1.0)
        @test_throws DomainError PSV.make_porous_shell_config(; kw..., phi_init=1.5)
        @test_throws DomainError PSV.make_porous_shell_config(; kw..., k_air=0.0)
    end

    @testset "Radial bin medians" begin
        r = [0.5, 0.6, 0.7, 2.5, 2.6, 9.9, 10.5]
        v = [1.0, 5.0, 2.0, 4.0, 8.0, 7.0, 99.0]
        c, m = PSV.radial_bin_medians(r, v, 10.0, 1.0)
        @test isapprox(c, [0.5, 2.5, 9.5]; atol=1e-12)
        @test isapprox(m, [2.0, 6.0, 7.0]; atol=1e-12)
        # R / dr = 10 / 3: the partial bin [9, 10) keeps the marker at 9.5, within dr of R,
        # with its centre at 9.5; markers at R, above R and below 0 are left out
        rp = [1.0, 9.5, 10.0, 10.5, 11.9, -0.5]
        vp = [3.0, 4.0, 50.0, 60.0, 70.0, 80.0]
        cp, mp = PSV.radial_bin_medians(rp, vp, 10.0, 3.0)
        @test isapprox(cp, [1.5, 9.5]; atol=1e-12)
        @test isapprox(mp, [3.0, 4.0]; atol=1e-12)
        @test length(cp) == 2
        # One value per radius, and a positive, finite radius and bin width
        @test_throws ArgumentError PSV.radial_bin_medians(r, v[1:(end - 1)], 10.0, 1.0)
        @test_throws ArgumentError PSV.radial_bin_medians(r[1:(end - 1)], v, 10.0, 1.0)
        for bad in (0.0, -1.0, NaN, Inf)
            @test_throws DomainError PSV.radial_bin_medians(r, v, bad, 1.0)
            @test_throws DomainError PSV.radial_bin_medians(r, v, 10.0, bad)
        end
    end

    @testset "Run clock: the last step and the final time of the case" begin
        cfg = PSV.make_porous_shell_config(;
            Nx=33, phi_init=0.25, t_form_yr=2.5e6, n_steps=80, dt_yr=5.0e4
        )
        dt_s = 5.0e4 * 31_557_600.0
        t_end = 2.5e6 * 31_557_600.0 + 80 * dt_s
        # 4.0 Myr after a start at 2.5 Myr: 6.5 Myr, about 2.05e14 s
        @test 2.0e14 < t_end < 2.1e14
        @test PSV.check_run_clock((; dt=dt_s, timesum=t_end), cfg) === nothing
        @test PSV.check_run_clock((; dt=dt_s, timesum=t_end * (1 + 1e-13)), cfg) === nothing
        # One step missing, a run in years and a last step cut by DTmax fail
        @test_throws ErrorException PSV.check_run_clock(
            (; dt=dt_s, timesum=t_end - dt_s), cfg
        )
        @test_throws ErrorException PSV.check_run_clock((; dt=dt_s, timesum=6.5e6), cfg)
        @test_throws ErrorException PSV.check_run_clock(
            (; dt=0.5 * dt_s, timesum=t_end), cfg
        )
        @test_throws ErrorException PSV.check_run_clock(
            (; dt=dt_s, timesum=t_end * (1 + 1e-11)), cfg
        )
    end

    @testset "Erebus and the reference heat alike before sintering" begin
        res = PSV.run_porous_shell_case(;
            Nx=33, phi_init=0.25, t_form_yr=2.5e6, n_steps=3, dt_yr=1.0e4, Nr=400
        )
        # The adiabatic centre: the same heat and heat capacity in both codes
        @test 300.0 < res.T_peak_ref < 310.0
        @test isapprox(res.T_peak, res.T_peak_ref; rtol=1e-8)
        @test isapprox(res.d_shell, PSV.R_PLANET; atol=1e-9)
        @test isapprox(res.d_shell_ref, PSV.R_PLANET; atol=1e-9)
        @test isapprox(res.dx, 1500.0; rtol=1e-14)
        @test res.phi_rms < 1e-12
    end

    @testset "Porous run: energy per step, time and radial convergence" begin
        run(Nr, dt_yr) = PSV.solve_porous_sphere(;
            R=20_000.0,
            phi_init=0.5,
            t_form_yr=1.5e6,
            t_end_yr=1.0e6,
            dt_yr=dt_yr,
            Nr=Nr,
        )
        a = run(400, 2.0e4)
        @test a.energy_residual < 1e-10
        @test 0.0 < maximum(a.dphi_sinter) <= 0.5
        # Order 1 in dt at the centre; near the sintering front k changes by up to 100x
        # and the largest difference only falls (not yet in the asymptotic range)
        b, c = run(400, 1.0e4), run(400, 5.0e3)
        @test 1.8 < abs(a.T[1] - b.T[1]) / abs(b.T[1] - c.T[1]) < 2.3
        e1 = maximum(abs.(a.T .- b.T))
        e2 = maximum(abs.(b.T .- c.T))
        @test 1.5 < e1 / e2 < 4.0
        # Radial refinement: the temperature at the centre converges
        r1, r2, r3 = run(100, 2.0e4), run(200, 2.0e4), a
        @test abs(r2.T[1] - r3.T[1]) < abs(r1.T[1] - r2.T[1])
        @test abs(r2.T[1] - r3.T[1]) < 1e-3 * (a.T[1] - PSV.T_SPACE)
    end

    @testset "Energy residual of a thermal step: balance, imbalance and scale guard" begin
        dt = 1.0e10
        rhocp = [1.0e6 2.0e6; 3.0e6 4.0e6]
        q = [1.0e-3 2.0e-3; 3.0e-3 4.0e-3]
        balanced = energy_test_grids(rhocp, dt .* q ./ rhocp, q)
        res = PSV.thermal_step_energy_residual(balanced, dt, 750.0, 750.0)
        @test isapprox(res, 0.0; atol=1.0e-15)
        # The boundary nodes (large RHOCP, tk2 and sources) do not enter.
        @test isapprox(balanced.RHOCP[1, 1], 3.0e10; rtol=1e-14)
        # Stored heat twice the source heat: |2S - S| / (2S + S) = 1/3
        doubled = energy_test_grids(rhocp, 2 .* dt .* q ./ rhocp, q)
        r2 = PSV.thermal_step_energy_residual(doubled, dt, 750.0, 750.0)
        @test isapprox(r2, 1 / 3; rtol=1e-12)
        # dx dy cancel; a missing Q_metric counts as 0
        @test isapprox(
            PSV.thermal_step_energy_residual(doubled, dt, 1.0, 2.0), r2; rtol=1e-12
        )
        no_metric = energy_test_grids(rhocp, dt .* q ./ rhocp, q; q_metric=false)
        @test isapprox(
            PSV.thermal_step_energy_residual(no_metric, dt, 750.0, 750.0), 0.0; atol=1e-15
        )
        # Opposite changes cancel in the sum, not in the scale: residual 1, not NaN
        swing = energy_test_grids(rhocp, [1.0 -0.5; 0.0 0.0], zeros(2, 2))
        @test isapprox(
            PSV.thermal_step_energy_residual(swing, dt, 750.0, 750.0), 0.0; atol=1e-15
        )
        lone = energy_test_grids(rhocp, [1.0 0.0; 0.0 0.0], zeros(2, 2))
        @test isapprox(
            PSV.thermal_step_energy_residual(lone, dt, 750.0, 750.0), 1.0; rtol=1e-14
        )
        still = energy_test_grids(rhocp, zeros(2, 2), zeros(2, 2))
        @test_throws ArgumentError PSV.thermal_step_energy_residual(still, dt, 1.0, 1.0)
    end

    @testset "Shell depth tolerance, metric references and metrics from runs" begin
        # One cell, or 15 % of the reference depth when that is larger
        @test isapprox(PSV.shell_depth_atol(33, 4000.0), 1500.0; rtol=1e-14)
        @test isapprox(PSV.shell_depth_atol(65, 4000.0), 750.0; rtol=1e-14)
        @test isapprox(PSV.shell_depth_atol(65, 10_000.0), 1500.0; rtol=1e-14)
        @test_throws DomainError PSV.shell_depth_atol(1, 4000.0)
        @test_throws DomainError PSV.shell_depth_atol(33, -1.0)
        @test_throws DomainError PSV.shell_depth_atol(33, NaN)
        # The CI case runs at Nx = 65 only: 3 metrics
        @test PSV.NX_LIST == (65,)
        refs = PSV.expected_metric_references(; d_shell_ref=4321.0)
        @test Set(keys(refs)) == Set(PSV.EXPECTED_METRIC_IDS)
        @test isapprox(refs["shell_depth_nx65"], 4321.0; rtol=1e-14)
        @test iszero(refs["tmax_rel_l2_nx65"]) && iszero(refs["energy_lse_residual_nx65"])
        lim = PSV.expected_metric_limits(; d_shell_ref=6000.0)
        @test Set(keys(lim)) == Set(PSV.EXPECTED_METRIC_IDS)
        @test isapprox(lim["shell_depth_nx65"].atol, 900.0; rtol=1e-14)
        @test lim["shell_depth_nx65"].comparison == "isapprox"
        @test isapprox(lim["tmax_rel_l2_nx65"].atol, 0.10; rtol=1e-14)
        @test isapprox(lim["energy_lse_residual_nx65"].atol, 1.0e-9; rtol=1e-14)
        @test all(l -> iszero(l.rtol), values(lim))
        @test_throws DomainError PSV.expected_metric_limits(; d_shell_ref=-1.0)
        fake(d, ref) =
            (; d_shell=d, d_shell_ref=ref, tmax_rel_l2=0.05, energy_lse_residual=1e-12)
        runs = Dict{Int,Any}(65 => fake(5500.0, 4000.0))
        m = PSV.metrics_from_runs(runs)
        @test Set(keys(m)) == Set(PSV.EXPECTED_METRIC_IDS)
        @test isapprox(m["shell_depth_nx65"], 5500.0; rtol=1e-14)
        @test isapprox(PSV.bench_shell_reference((; runs, metrics=m)), 4000.0; rtol=1e-14)
        bad = Dict{Int,Any}(65 => fake(5500.0, 4000.0), 129 => fake(5000.0, 4100.0))
        @test_throws ArgumentError PSV.bench_shell_reference((; runs=bad))
        @test isempty(PSV.metrics_from_runs(Dict{Int,Any}()))
        @test_throws ArgumentError PSV.evaluate_porous_shell_benchmark(; nx_list=(33,))
        @test_throws ArgumentError PSV.evaluate_porous_shell_benchmark(; nx_list=(65.0,))
    end

    @testset "Reference literals match the configuration and the code" begin
        ok_lit, errs_lit = PSV.check_reference_literals_match_config()
        @test ok_lit
        @test isempty(errs_lit)
        base = PSV.make_case_config()
        changed = SimulationConfig(;
            thermodynamics=ThermalConfig(; ratio_al=6.0e-5),
            porosity=PorosityConfig(; model=:henke2012, A=5.0e-5),
            disk=base.disk,
            materials=base.materials,
        )
        ok_bad, errs_bad = PSV.check_reference_literals_match_config(changed; x_fe_ref=0.19)
        @test !ok_bad
        @test length(errs_bad) == 3
        @test any(occursin("RATIO_AL_L16", e) for e in errs_bad)
        @test any(occursin("CREEP_A", e) for e in errs_bad)
        @test any(occursin("X_FE", e) for e in errs_bad)
    end

    @testset "Suite record of the 6 full runs: static check and export" begin
        suite_path = joinpath(
            @__DIR__, "..", "benchmarks", "records", "porous_shell_suite.toml"
        )
        raw_suite = TOML.parsefile(suite_path)
        ok_suite, errs_suite = PSV.check_suite_record_static(suite_path)
        @test ok_suite
        @test isempty(errs_suite)
        @test raw_suite["metadata"]["config_sha256"] == PSV.suite_config_sha256()
        @test Set(keys(raw_suite["config"])) ==
            Set(PSV.REQUIRED_CONFIG_KEYS["porous_shell_suite"])
        # 6 runs, a shell depth each and a temperature error at Nx >= 65 each: 12 metrics
        @test raw_suite["config"]["tags"] == [String(c[1]) for c in PSV.SUITE_CASES]
        @test length(PSV.SUITE_CASES) == 6
        @test length(raw_suite["metrics"]) == 12
        for m in raw_suite["metrics"]
            nx, ref = m["nx"], Float64(m["reference"])
            if startswith(m["id"], "shell_depth_")
                # One cell of 48 km / (Nx - 1) or 15 % of the 1-D depth
                @test isapprox(m["atol"], max(48_000.0 / (nx - 1), 0.15 * ref); rtol=1e-14)
                @test 0.0 < ref < 20_000.0
                @test abs(m["value"] - ref) <= m["atol"]
            else
                @test isapprox(m["atol"], 0.10; rtol=1e-14) && m["comparison"] == "le_atol"
                @test 0.0 <= m["value"] <= 0.10
            end
        end
        # The suite has no run below Nx = 65; a run at Nx = 33 would have no tmax metric
        @test all(c -> c[2] >= 65, PSV.SUITE_CASES)
        @test PSV.suite_metric_ids()[1:2] ==
            ["shell_depth_s65_25_025", "tmax_rel_l2_s65_25_025"]
        @test_throws ArgumentError PSV.suite_metric_limit("other_s65", 65, 4000.0)
        mktempdir() do tmpdir
            # Tampered records fail the static check
            function static_errors(f)
                raw = deepcopy(raw_suite)
                f(raw)
                path = write_porous_shell_toml(joinpath(tmpdir, "suite.toml"), raw)
                return PSV.check_suite_record_static(path; check_doc_sync=false)
            end
            ok_same, errs_same = static_errors(identity)
            @test ok_same && isempty(errs_same)
            ok_w, errs_w = static_errors(r -> (r["metrics"][1]["atol"] *= 2.0))
            @test !ok_w
            @test any(occursin("differs from the rule", e) for e in errs_w)
            ok_t, errs_t = static_errors(r -> (r["metrics"][2]["atol"] = 0.5))
            @test !ok_t
            @test any(occursin("tmax_rel_l2_s65_25_025", e) for e in errs_t)
            ok_c, errs_c = static_errors(r -> deleteat!(r["metrics"], 11:12))
            @test !ok_c
            @test any(occursin("case list", e) for e in errs_c)
            @test any(occursin("SUITE_CASES", e) for e in errs_c)
            ok_h, errs_h = static_errors(
                r -> (r["metadata"]["config_sha256"] = repeat("0", 64))
            )
            @test !ok_h
            @test any(occursin("config_sha256", e) for e in errs_h)
            ok_n, errs_n = static_errors(r -> (r["config"]["n_steps"][2] = 740))
            @test !ok_n
            @test any(occursin("config.n_steps", e) for e in errs_n)
            ok_x, errs_x = static_errors(r -> (r["metrics"][3]["nx"] = 65))
            @test !ok_x
            @test any(occursin("has nx 65", e) for e in errs_x)
            ok_v, errs_v = static_errors(r -> (r["metrics"][3]["value"] += 1.0e4))
            @test !ok_v
            @test any(occursin("fails comparison", e) for e in errs_v)
            # A 1-D depth other than the pinned one fails, also when the value passes it
            ok_r, errs_r = static_errors(r -> (r["metrics"][1]["reference"] = 4660.0))
            @test !ok_r
            @test any(
                occursin("shell_depth_s65_25_025", e) && occursin("pinned 1-D depth", e) for
                e in errs_r
            )
            # A shell depth outside [0, R_PLANET] fails, also within its own tolerance
            function out_of_body(value, ref)
                return r -> begin
                    m = r["metrics"][1]
                    m["reference"], m["value"] = ref, value
                    m["atol"] = PSV.shell_depth_atol(65, ref)
                end
            end
            for (value, ref) in ((20_500.0, 19_900.0), (-300.0, 100.0))
                ok_o, errs_o = static_errors(out_of_body(value, ref))
                @test !ok_o
                @test any(occursin("outside [0, 20000.0]", e) for e in errs_o)
            end
            @test !first(static_errors(out_of_body(NaN, 4560.0)))
            # The pinned depths cover the cases, lie inside the body and enter the hash
            pinned_tags = [t for (t, _) in PSV.SUITE_REFERENCE_DEPTHS]
            @test pinned_tags == [String(c[1]) for c in PSV.SUITE_CASES]
            @test all(0.0 < d < PSV.R_PLANET for (_, d) in PSV.SUITE_REFERENCE_DEPTHS)
            @test occursin("case.SUITE_REFERENCE_DEPTHS=", PSV.suite_config_text())
            for m in raw_suite["metrics"]
                startswith(m["id"], "shell_depth_") || continue
                @test isapprox(
                    m["reference"], PSV.suite_reference_depth(m["case"]); atol=1e-9
                )
            end
            @test_throws ArgumentError PSV.suite_reference_depth("s33_other")
            @test !first(PSV.check_suite_record_static(suite_path * ".missing"))
            @test_throws DomainError PSV.check_suite_record_static("")

            # Static branch of the command line
            cli = capture_porous_shell_stdio(
                () -> run_validation_cli(["--check", suite_path])
            )
            @test cli.value == 0
            @test occursin("static check", cli.out)
            @test static_validation_case("porous_shell_suite").record ==
                normpath(suite_path)
            wide = deepcopy(raw_suite)
            wide["metrics"][1]["atol"] *= 2.0
            bad_path = write_porous_shell_toml(joinpath(tmpdir, "bad_suite.toml"), wide)
            cli_bad = capture_porous_shell_stdio(() -> run_validation_cli([bad_path]))
            @test cli_bad.value == 1
            @test occursin("[static]", cli_bad.err)

            # Export from a folder of result files reproduces the record
            folder = write_suite_folder(mkpath(joinpath(tmpdir, "suite")), raw_suite)
            out = joinpath(tmpdir, "exported_suite.toml")
            commit = String(raw_suite["metadata"]["src_commit"])
            SE = PorousShellSuiteExport
            @test SE.run_suite_export_cli(
                [folder, out]; src_commit=commit, check_clean_tree=false
            ) == 0
            raw_out = TOML.parsefile(out)
            @test first(PSV.check_suite_record_static(out; check_doc_sync=false))
            @test PSV.render_validation_markdown_block(raw_out) ==
                PSV.render_validation_markdown_block(raw_suite)
            @test raw_out["config"] == raw_suite["config"]
            # A row of another case, a missing file and a bad command line are refused
            bad = write_suite_folder(
                mkpath(joinpath(tmpdir, "bad")),
                raw_suite;
                tamper=Dict("s129_25_050" => ("Nx", 65)),
            )
            @test_throws ArgumentError SE.read_suite_results(bad)
            err_run = capture_porous_shell_stdio(
                () -> SE.run_suite_export_cli(
                    [bad, joinpath(tmpdir, "x.toml")];
                    src_commit=commit,
                    check_clean_tree=false,
                ),
            )
            @test err_run.value == 1
            @test !isfile(joinpath(tmpdir, "x.toml"))
            # Uncommitted changes in src/, benchmarks/ or tools/ stop the export
            dirty = capture_porous_shell_stdio(
                () -> SE.run_suite_export_cli(
                    [folder, joinpath(tmpdir, "dirty.toml")];
                    src_commit=commit,
                    dirty_status=" M benchmarks/validation_porous_shell.jl",
                ),
            )
            @test dirty.value == 1
            @test occursin("uncommitted changes", dirty.err)
            @test !isfile(joinpath(tmpdir, "dirty.toml"))
            @test first(
                PorousShellExport.check_clean_source_tree(;
                    record="benchmarks/records/porous_shell_suite.toml"
                ),
            ) isa Bool
            # An energy residual above ENERGY_RESIDUAL_LIMIT (1e-9) is refused
            hot = write_suite_folder(
                mkpath(joinpath(tmpdir, "hot")),
                raw_suite;
                tamper=Dict("s129_25_050" => ("energy_res", 2.0e-9)),
            )
            e_hot = try
                SE.read_suite_results(hot)
            catch err
                err
            end
            @test e_hot isa ArgumentError
            @test occursin("energy_res", e_hot.msg)
            nan_e = write_suite_folder(
                mkpath(joinpath(tmpdir, "nan_e")),
                raw_suite;
                tamper=Dict("s65_25_025" => ("energy_res", NaN)),
            )
            @test_throws ArgumentError SE.read_suite_results(nan_e)
            # A 1-D depth other than the pinned depth of the case is refused
            moved_ref = write_suite_folder(
                mkpath(joinpath(tmpdir, "moved_ref")),
                raw_suite;
                tamper=Dict("s129_225_025" => ("d_ref", 2400.0)),
            )
            e_ref = try
                SE.read_suite_results(moved_ref)
            catch err
                err
            end
            @test e_ref isa ArgumentError
            @test occursin("pinned", e_ref.msg)
            # A result file of a tag outside SUITE_CASES is refused
            stray = write_suite_folder(mkpath(joinpath(tmpdir, "stray")), raw_suite)
            cp(joinpath(stray, "result_s65_25_025.csv"), joinpath(stray, "result_s33.csv"))
            e_stray = try
                SE.read_suite_results(stray)
            catch err
                err
            end
            @test e_stray isa ArgumentError
            @test occursin("s33", e_stray.msg)
            rm(joinpath(folder, "result_s129_25_025_k.csv"))
            @test_throws ArgumentError SE.read_suite_results(folder)
            usage = capture_porous_shell_stdio(() -> SE.run_suite_export_cli(String[]))
            @test usage.value == 2
            @test SE.run_suite_export_cli(["--x", folder]; src_commit=commit) == 2
            @test_throws DomainError SE.write_suite_record(
                folder; record_path=out, src_commit="abc"
            )
        end
    end

    record_path = joinpath(@__DIR__, "..", "benchmarks", "records", "porous_shell.toml")
    raw_record_all = TOML.parsefile(record_path)
    recorded = Dict(
        String(m["id"]) => Float64(m["value"]) for m in raw_record_all["metrics"]
    )

    @testset "Reproducibility check: a metric without a limit or a reference fails" begin
        # A shell depth at Nx = 129, outside NX_LIST: the code has no limit and no
        # reference for it, so its recorded tolerance cannot be checked
        extra_rec = deepcopy(raw_record_all)
        extra = deepcopy(extra_rec["metrics"][1])
        extra["id"], extra["nx"] = "shell_depth_nx129", 129
        push!(extra_rec["metrics"], extra)
        vals = Dict(String(m["id"]) => Float64(m["value"]) for m in extra_rec["metrics"])
        fake_bench = (; runs=Dict{Int,Any}(), metrics=vals)
        mktempdir() do tmpdir
            path = write_porous_shell_toml(joinpath(tmpdir, "extra.toml"), extra_rec)
            ok_x, errs_x = PSV.check_validation_reproducibility(
                fake_bench, path; check_recorded_values=false, d_shell_ref=4590.0
            )
            @test !ok_x
            @test any(
                occursin("'shell_depth_nx129'", e) && occursin("no limit", e) for
                e in errs_x
            )
            @test any(
                occursin("'shell_depth_nx129'", e) &&
                    occursin("no recomputed reference", e) for e in errs_x
            )
            # The 3 metrics of the code raise neither error
            @test !any(
                occursin("nx65", e) && occursin(r"no (limit|recomputed)", e) for e in errs_x
            )
        end
    end

    bench = PSV.evaluate_porous_shell_benchmark()

    @testset "Benchmark run: every metric, physical bounds of the runs" begin
        @test Set(keys(bench.metrics)) == Set(PSV.EXPECTED_METRIC_IDS)
        for n in PSV.NX_LIST
            r = bench.runs[n]
            # A sintered core under an unsintered shell, in both codes
            @test 0.0 < r.d_shell_ref < PSV.R_PLANET
            @test 0.0 < r.d_shell < PSV.R_PLANET
            @test PSV.T_SINTER < r.T_peak < 1500.0
            @test isapprox(r.T_peak, r.T_peak_ref; rtol=0.02)
        end
        # The shell lies within one cell of 750 m or 15 % of the 1-D depth
        r65 = bench.runs[65]
        @test abs(r65.d_shell - r65.d_shell_ref) <= max(750.0, 0.15 * r65.d_shell_ref)
        @test bench.metrics["energy_lse_residual_nx65"] <= 1.0e-9
        # The 1-D shell of the CI case is 4590 m deep (459 cells of 10 m), and the
        # recorded Erebus shell is 5250 m deep (14 bins of 375 m from the surface)
        @test isapprox(r65.d_shell_ref, 4590.0; atol=1e-6)
        @test isapprox(recorded["shell_depth_nx65"], 5250.0; atol=1e-9)
        @test isapprox(r65.d_shell, recorded["shell_depth_nx65"]; atol=375.0)
    end

    @testset "Validation record schema and drift guards" begin
        @test isfile(record_path)
        raw_record = TOML.parsefile(record_path)
        valid_ok, valid_errors = PSV.validate_record_dict(
            raw_record; require_doc_pages=true
        )
        @test valid_ok
        @test isempty(valid_errors)

        rec_ids = [String(m["id"]) for m in raw_record["metrics"]]
        @test rec_ids == collect(PSV.EXPECTED_METRIC_IDS)
        @test issubset(keys(PSV.PHYSICS_REGRESSION_BANDS), PSV.EXPECTED_METRIC_IDS)
        @test all(m["case"] == PSV.CASE_ID for m in raw_record["metrics"])

        missing_key = deepcopy(raw_record)
        delete!(missing_key["config"], "phi_init")
        ok_missing, errs_missing = PSV.validate_record_dict(
            missing_key; require_doc_pages=false, check_doc_sync=false
        )
        @test !ok_missing
        @test any(occursin("phi_init", e) for e in errs_missing)

        # The energy residual has no band
        @test sort!(collect(keys(PSV.PHYSICS_REGRESSION_BANDS))) ==
            ["shell_depth_nx65", "tmax_rel_l2_nx65"]
        for (mid, value) in bench.metrics
            haskey(PSV.PHYSICS_REGRESSION_BANDS, mid) || continue
            band = PSV.PHYSICS_REGRESSION_BANDS[mid]
            @test isapprox(value, recorded[mid]; atol=band.atol, rtol=band.rtol)
        end
        ok_full, errs_full = PSV.check_physics_regression_bands(bench.metrics, record_path)
        @test ok_full
        @test isempty(errs_full)

        target = "tmax_rel_l2_nx65"
        moved = copy(bench.metrics)
        moved[target] += 2 * PSV.PHYSICS_REGRESSION_BANDS[target].atol + 1.0e-6
        ok_moved, errs_moved = PSV.check_physics_regression_bands(moved, record_path)
        @test !ok_moved
        @test length(errs_moved) == 1
        @test occursin(target, only(errs_moved))
        missing_metric = copy(bench.metrics)
        delete!(missing_metric, target)
        ok_gone, errs_gone = PSV.check_physics_regression_bands(missing_metric, record_path)
        @test !ok_gone
        @test any(occursin("Measured metrics missing '$target'", e) for e in errs_gone)
        @test_throws DomainError PSV.check_physics_regression_bands(bench.metrics, "")
        @test !first(PSV.check_physics_regression_bands(bench.metrics, record_path * ".x"))
    end

    @testset "Record matches the configuration, the tolerances and the references" begin
        raw_record = TOML.parsefile(record_path)
        cfg = PSV.make_case_config()
        sha = PSV.canonical_config_sha256(cfg)
        @test sha == raw_record["metadata"]["config_sha256"]
        @test length(sha) == 64
        @test PSV.canonical_config_sha256(PSV.make_case_config(; seed=1)) != sha
        @test PSV.canonical_config_sha256(PSV.make_case_config(; Nx=33)) != sha

        text = PSV.canonical_config_text(cfg)
        for name in PSV.HASHED_CASE_CONSTANTS
            @test occursin("case.$(name)=", text)
        end
        @test occursin("porosity.model=:henke2012", text)
        @test bytes2hex(SHA.sha256(text)) == sha

        summary = PSV.expected_config_summary(cfg)
        rec_cfg = raw_record["config"]
        @test Set(keys(summary)) == Set(keys(rec_cfg))
        @test Set(keys(summary)) == Set(PSV.REQUIRED_CONFIG_KEYS["porous_shell"])
        for (key, expected) in summary
            if expected isa AbstractFloat
                @test isapprox(Float64(rec_cfg[key]), expected; rtol=1.0e-12, atol=0.0)
            else
                @test rec_cfg[key] == expected
            end
        end

        # The record holds the a priori limits and the 1-D reference of the same run
        d_ref = PSV.bench_shell_reference(bench)
        refs = PSV.expected_metric_references(; d_shell_ref=d_ref)
        limits = Dict(
            "shell_depth_nx65" => (max(750.0, 0.15 * d_ref), "isapprox"),
            "tmax_rel_l2_nx65" => (0.10, "le_atol"),
            "energy_lse_residual_nx65" => (1.0e-9, "le_atol"),
        )
        for m in raw_record["metrics"]
            mid = String(m["id"])
            @test isapprox(Float64(m["reference"]), refs[mid]; rtol=1.0e-12, atol=1.0e-12)
            @test isapprox(Float64(m["atol"]), limits[mid][1]; rtol=1.0e-14)
            @test iszero(Float64(m["rtol"]))
            @test m["comparison"] == limits[mid][2]
        end
        @test isapprox(d_ref, PSV.reference_solution().d_shell; rtol=1e-14)
    end

    @testset "Reproducibility check, CLI table and export round trip" begin
        raw_record = TOML.parsefile(record_path)
        ok_repro, errs_repro = PSV.check_validation_reproducibility(bench, record_path)
        @test ok_repro
        @test isempty(errs_repro)
        @test run_porous_shell_validation_check(record_path; bench=bench)
        @test run_validation_cli(["--check", record_path]; porous_shell_bench=bench) == 0
        @test !run_porous_shell_validation_check(record_path * ".missing"; bench=bench)
        @test_throws DomainError run_porous_shell_validation_check(""; bench=bench)
        @test validation_case("porous_shell").bench_keyword == :porous_shell_bench

        mktempdir() do tmpdir
            empty_bench = (; runs=Dict{Int,Any}(), metrics=Dict{String,Float64}())
            no_arg = capture_porous_shell_stdio(
                () -> run_validation_cli(
                    ["--check"];
                    bench=empty_bench,
                    darcy_bench=empty_bench,
                    compaction_bench=empty_bench,
                    diffusion_bench=empty_bench,
                    momentum_bench=empty_bench,
                    porous_shell_bench=bench,
                ),
            )
            @test no_arg.value == 1
            @test occursin("porous_shell.toml", no_arg.out)
            @test occursin("[reproducibility]", no_arg.err)
            @test last(VALIDATION_CASES).id == "porous_shell"

            # A tampered record: a moved value, a changed configuration and its hash
            tampered = deepcopy(raw_record)
            i_tmax = findfirst(m -> m["id"] == "tmax_rel_l2_nx65", tampered["metrics"])
            tampered["metrics"][i_tmax]["value"] += 0.01
            tampered["config"]["phi_init"] = 0.3
            tampered["metadata"]["config_sha256"] = repeat("0", 64)
            t_path = write_porous_shell_toml(joinpath(tmpdir, "tampered.toml"), tampered)
            ok_t, errs_t = PSV.check_validation_reproducibility(bench, t_path)
            @test !ok_t
            @test any(occursin("config.phi_init", e) for e in errs_t)
            @test any(occursin("config_sha256", e) for e in errs_t)
            @test any(occursin("differs from recorded value", e) for e in errs_t)
            @test first(
                PSV.check_validation_reproducibility(
                    bench, record_path; check_recorded_values=false
                ),
            )
            wrong_ref = deepcopy(raw_record)
            wrong_ref["metrics"][1]["reference"] += 100.0
            wrong_ref["metrics"][1]["value"] = wrong_ref["metrics"][1]["reference"]
            w_path = write_porous_shell_toml(joinpath(tmpdir, "wrong_ref.toml"), wrong_ref)
            ok_w, errs_w = PSV.check_validation_reproducibility(bench, w_path)
            @test !ok_w
            @test any(occursin("recomputed reference", e) for e in errs_w)

            # A record whose tolerance is wider than the limit of the code fails
            for (key, widened) in (("atol", 0.5), ("rtol", 0.1), ("comparison", "isapprox"))
                wide = deepcopy(raw_record)
                wide["metrics"][i_tmax][key] = widened
                wide["metrics"][i_tmax]["reference"] = 0.0
                wide_path = write_porous_shell_toml(joinpath(tmpdir, "wide.toml"), wide)
                ok_wide, errs_wide = PSV.check_validation_reproducibility(bench, wide_path)
                @test !ok_wide
                @test any(
                    occursin("differs from the limit of the code", e) for e in errs_wide
                )
            end
            i_d = findfirst(m -> m["id"] == "shell_depth_nx65", raw_record["metrics"])
            wide_d = deepcopy(raw_record)
            wide_d["metrics"][i_d]["atol"] = 2.0 * wide_d["metrics"][i_d]["atol"]
            wd_path = write_porous_shell_toml(joinpath(tmpdir, "wide_d.toml"), wide_d)
            ok_wd, errs_wd = PSV.check_validation_reproducibility(bench, wd_path)
            @test !ok_wd
            @test any(
                occursin("shell_depth_nx65", e) && occursin("limit of the code", e) for
                e in errs_wd
            )

            # Export round trip
            exported = joinpath(tmpdir, "exported_porous_shell.toml")
            commit = String(raw_record["metadata"]["src_commit"])
            @test PorousShellExport.run_export_cli(
                [exported]; src_commit=commit, bench=bench, check_clean_tree=false
            ) == 0
            raw_exported = TOML.parsefile(exported)
            ok_valid, errs_valid = PSV.validate_record_dict(
                raw_exported; require_doc_pages=false, check_doc_sync=false
            )
            @test ok_valid
            @test isempty(errs_valid)
            @test run_porous_shell_validation_check(exported; bench=bench)
            @test [String(m["id"]) for m in raw_exported["metrics"]] == [String(m["id"]) for m in raw_record["metrics"]]

            # A moved value is refused without the flag and written with it
            target = "tmax_rel_l2_nx65"
            band = PSV.PHYSICS_REGRESSION_BANDS[target]
            moved_metrics = copy(bench.metrics)
            moved_metrics[target] += 2 * band.atol + 1.0e-6
            moved_bench = (; runs=bench.runs, metrics=moved_metrics)
            @test PorousShellExport.run_export_cli(
                [exported]; src_commit=commit, bench=moved_bench, check_clean_tree=false
            ) == 1
            @test PorousShellExport.run_export_cli(
                ["--accept-physics-change", exported];
                src_commit=commit,
                bench=moved_bench,
                check_clean_tree=false,
            ) == 0
            @test PorousShellExport.run_export_cli(
                ["--bad-flag"]; src_commit=commit, bench=bench, check_clean_tree=false
            ) == 2
            @test PorousShellExport.run_export_cli(
                [exported];
                src_commit=commit,
                bench=bench,
                check_clean_tree=true,
                dirty_status=" M src/config.jl",
            ) == 1
            @test_throws DomainError PorousShellExport.run_export_cli(
                [exported]; src_commit="abc", bench=bench, check_clean_tree=false
            )
            # A value that fails its own tolerance is never written
            failing = copy(bench.metrics)
            failing["energy_lse_residual_nx65"] = 1.0e-6
            @test PorousShellExport.run_export_cli(
                ["--accept-physics-change", exported];
                src_commit=commit,
                bench=(; runs=bench.runs, metrics=failing),
                check_clean_tree=false,
            ) == 1

            status = PorousShellExport.metric_status
            cases = (
                (0.01, 0.0, "le_atol"), (0.01, 1.0, "le_atol"), (1.0, 1.001, "isapprox")
            )
            @test isequal(
                [status(v, ref, 0.05, 0.0, c) for (v, ref, c) in cases],
                ["PASS", "FAIL", "PASS"],
            )
            @test isequal(status(1.0, 1.0, 0.0, 0.0, "other"), "FAIL")
            @test occursin(r"^[0-9a-f]{40}$", PorousShellExport.last_src_commit())
            @test first(PorousShellExport.check_clean_source_tree()) isa Bool
            @test PorousShellExport.porous_shell_validation_module() === PSV
            @test_throws ErrorException PorousShellExport.porous_shell_validation_module((
                Module(:Empty),
            ))
        end
    end

    @testset "Checks stop on an interrupt and report other errors with a stack trace" begin
        commit = String(raw_record_all["metadata"]["src_commit"])
        @test_throws InterruptException run_porous_shell_validation_check(
            record_path; bench=PorousShellInterruptingBench()
        )
        crash = capture_porous_shell_stdio(
            () -> run_porous_shell_validation_check(
                record_path; bench=PorousShellCrashingBench()
            ),
        )
        @test crash.value === false
        @test occursin("solver failure inside the benchmark", crash.err)
        @test occursin("Stacktrace", crash.err)
        mktempdir() do tmpdir
            target = joinpath(tmpdir, "record.toml")
            @test_throws InterruptException PorousShellExport.run_export_cli(
                [target];
                src_commit=commit,
                bench=PorousShellInterruptingBench(),
                check_clean_tree=false,
            )
            failed = capture_porous_shell_stdio(
                () -> PorousShellExport.run_export_cli(
                    [target];
                    src_commit=commit,
                    bench=PorousShellCrashingBench(),
                    check_clean_tree=false,
                ),
            )
            @test failed.value == 1
            @test occursin("Stacktrace", failed.err)
            @test !isfile(target)
        end
    end

    @testset "Command lines and malformed records are reported, not thrown" begin
        commit = String(raw_record_all["metadata"]["src_commit"])
        mktempdir() do tmpdir
            p1, p2 = joinpath(tmpdir, "first.toml"), joinpath(tmpdir, "second.toml")
            @test PorousShellExport.run_export_cli(
                [p1, p2]; src_commit=commit, bench=bench, check_clean_tree=false
            ) == 2
            @test !isfile(p1) && !isfile(p2)
            empty_run = capture_porous_shell_stdio(
                () -> PorousShellExport.export_main(
                    [""]; src_commit=commit, bench=bench, check_clean_tree=false
                ),
            )
            @test empty_run.value == 1
            @test occursin("Record path must be non-empty", empty_run.err)
            broken = joinpath(tmpdir, "broken.toml")
            write(broken, "[metadata\nid = ")
            for check in (
                () -> PSV.check_validation_reproducibility(bench, broken),
                () -> PSV.check_physics_regression_bands(bench.metrics, broken),
            )
                ok_parse, errs_parse = check()
                @test !ok_parse
                @test occursin("Cannot parse", only(errs_parse))
            end
            bad_value = deepcopy(raw_record_all)
            bad_value["metrics"][1]["value"] = "not a number"
            bad_path = write_porous_shell_toml(
                joinpath(tmpdir, "bad_value.toml"), bad_value
            )
            ok_value, errs_value = PSV.check_validation_reproducibility(bench, bad_path)
            @test !ok_value
            @test any(occursin("numeric", e) for e in errs_value)
            cli = capture_porous_shell_stdio(
                () -> run_validation_cli([bad_path]; porous_shell_bench=bench)
            )
            @test cli.value == 1
            @test occursin("[schema]", cli.err)
            @test_throws DomainError PSV.check_validation_reproducibility(bench, "")
        end
    end
end
