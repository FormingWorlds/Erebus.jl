using Erebus
using Erebus.Config
using Erebus.Physics
using Erebus.Geometry
using Erebus.Particles
using Test

# Standard erf implementation (Abramowitz & Stegun 7.1.26, error < 1.5e-7)
function _erf(x::Real)
    p = 0.3275911
    a1 = 0.254829592
    a2 = -0.284496736
    a3 = 1.421413741
    a4 = -1.453152027
    a5 = 1.061405429
    sign_x = sign(x)
    abs_x = abs(x)
    t = 1.0 / (1.0 + p * abs_x)
    y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-abs_x * abs_x)
    return sign_x * y
end

# Solve transcendental Stefan eigenvalue equation: λ * exp(λ^2) * erf(λ) = Ste / sqrt(pi)
function solve_stefan_eigenvalue(Ste::Real; tol=1e-10, maxiter=100)
    target = Ste / sqrt(pi)
    lo = 0.0001
    hi = 5.0
    for _ in 1:maxiter
        mid = 0.5 * (lo + hi)
        f_mid = mid * exp(mid^2) * _erf(mid) - target
        if abs(f_mid) < tol
            return mid
        elseif f_mid < 0.0
            lo = mid
        else
            hi = mid
        end
    end
    return 0.5 * (lo + hi)
end

@testset "1D Stefan Moving-Boundary Benchmark" begin
    @testset "Transcendental Eigenvalue and Front Position Scaling" begin
        # Test eigenvalue solver for known Stefan numbers
        λ_05 = solve_stefan_eigenvalue(0.5)
        λ_10 = solve_stefan_eigenvalue(1.0)
        @test isapprox(λ_05 * exp(λ_05^2) * _erf(λ_05), 0.5 / sqrt(pi); atol=1e-8)
        @test isapprox(λ_10 * exp(λ_10^2) * _erf(λ_10), 1.0 / sqrt(pi); atol=1e-8)
        @test λ_10 > λ_05

        # Diffusive scaling: front position s(t) = 2*λ*sqrt(κ*t)
        kappa = 1.0e-6 # m^2/s
        s_t1 = 2.0 * λ_05 * sqrt(kappa * 1000.0)
        s_t4 = 2.0 * λ_05 * sqrt(kappa * 4000.0)
        # s(4t) / s(t) must equal 2.0 exactly
        @test isapprox(s_t4 / s_t1, 2.0; rtol=1e-10)

        # Front velocity ds/dt = λ*sqrt(κ/t) scales as t^(-1/2)
        v_t1 = λ_05 * sqrt(kappa / 1000.0)
        v_t4 = λ_05 * sqrt(kappa / 4000.0)
        @test isapprox(v_t1 / v_t4, 2.0; rtol=1e-10)
    end

    @testset "1D Column Equilibrium Moving-Boundary Simulation" begin
        # Set up 1D column on a staggered grid
        Ny_cells = 65
        Nx_cells = 5
        H_domain = 0.1 # 0.1 m column
        W_domain = 0.01 # 0.01 m width
        cfg_grid = GridConfig(; Nx=Nx_cells, Ny=Ny_cells, xsize=W_domain, ysize=H_domain)
        coords = GridCoordinates(cfg_grid)

        # Rock thermophysical properties
        rho = 3000.0
        cp = 1000.0
        k_therm = 3.0
        kappa = k_therm / (rho * cp) # 1.0e-6 m^2/s

        # Dehydration reaction thermochemical parameters
        dH = 40_000.0
        dS = 60.0
        T_eq = dH / dS # 666.67 K at Pf = 0
        T_wall = 800.0 # Top heating boundary

        # Effective latent heat: L_eff = ΔXW * dH / (MD + MH2O)
        MD_mol = 0.120
        MH2O_mol = 0.018
        delta_XW = 0.90
        L_eff = delta_XW * dH / (MD_mol + MH2O_mol)

        # Stefan number and transcendental eigenvalue
        Ste = cp * (T_wall - T_eq) / L_eff
        λ_ana = solve_stefan_eigenvalue(Ste)

        # Set up markers in 1D column
        marknum = 200
        props = setup_marker_properties(marknum, coords)
        (xm, ym, tm, tkm, sxxm, sxym, etavpm, phim, phinewm, pfm0, XWsolidm, XWsolidm0) =
            props

        tm .= 1
        xm .= W_domain / 2
        dy_marker = H_domain / marknum
        for m in 1:marknum
            ym[m] = (m - 0.5) * dy_marker
        end
        phim .= 0.15
        phinewm .= 0.15
        pfm0 .= 1.0e5
        XWsolidm0 .= 0.95
        XWsolidm .= 0.95

        t_elapsed = 200.0 # 200 s
        s_ana = 2.0 * λ_ana * sqrt(kappa * t_elapsed)

        # Impose Stefan similarity temperature profile on grid and markers
        tk2 = fill(T_eq, coords.Ny1, coords.Nx1)
        pf = fill(1.0e5, coords.Ny1, coords.Nx1)
        for j in 1:coords.Nx1, i in 1:coords.Ny1
            y_node = (i - 1) * coords.dy
            if y_node <= s_ana
                tk2[i, j] =
                    T_wall -
                    (T_wall - T_eq) * _erf(y_node / (2.0 * sqrt(kappa * t_elapsed))) /
                    _erf(λ_ana)
            else
                tk2[i, j] = T_eq
            end
        end

        for m in 1:marknum
            if ym[m] <= s_ana
                tkm[m] =
                    T_wall -
                    (T_wall - T_eq) * _erf(ym[m] / (2.0 * sqrt(kappa * t_elapsed))) /
                    _erf(λ_ana)
            else
                tkm[m] = T_eq
            end
        end

        DMP = zeros(Float64, coords.Ny1, coords.Nx1)
        DHP = zeros(Float64, coords.Ny1, coords.Nx1)
        DQPF = zeros(Float64, coords.Ny1, coords.Nx1)
        DMPSUM = zeros(Float64, coords.Ny1, coords.Nx1)
        DHPSUM = zeros(Float64, coords.Ny1, coords.Nx1)
        DQPFSUM = zeros(Float64, coords.Ny1, coords.Nx1)
        WTPSUM = zeros(Float64, coords.Ny1, coords.Nx1)

        cfg_react = ReactionConfig(;
            delta_H=dH,
            delta_S=dS,
            dtreaction_dehydration=1.0, # Fast equilibrium reaction limit
            dtreaction_hydration=1.0,
        )

        perform_thermochemical_reaction!(
            DMP,
            DHP,
            DMPSUM,
            DHPSUM,
            WTPSUM,
            pf,
            tk2,
            tm,
            xm,
            ym,
            XWsolidm0,
            XWsolidm,
            phim,
            phinewm,
            pfm0,
            marknum,
            50.0,
            2,
            1;
            coords=coords,
            DQPF=DQPF,
            DQPFSUM=DQPFSUM,
            cfg=cfg_react,
        )

        # 1. Markers behind front (y < s_ana) dehydrate; markers ahead remain hydrous
        dehydrated_indices = findall(m -> XWsolidm[m] < 0.5, 1:marknum)
        hydrated_indices = findall(m -> XWsolidm[m] >= 0.5, 1:marknum)
        @test !isempty(dehydrated_indices)
        @test !isempty(hydrated_indices)

        # 2. Numerical front position from marker state matches analytical Stefan front within marker spacing
        s_num = maximum(ym[dehydrated_indices])
        @test abs(s_num - s_ana) <= 2.0 * dy_marker

        # 3. Latent heat sink in dehydrated zone: DHP < 0
        @test minimum(DHP) < 0.0
        @test maximum(DHP) <= 1e-12

        # 4. Fluid expulsion in dehydrated zone: DQPF > 0
        @test maximum(DQPF) > 0.0
        @test minimum(DQPF) >= -1e-12
    end
end
