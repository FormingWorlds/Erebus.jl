using Test
using Erebus

isdefined(@__MODULE__, :PorousShellValidation) ||
    include(joinpath(@__DIR__, "..", "benchmarks", "validation_porous_shell.jl"))
const PSV = PorousShellValidation

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
end
