@testset "Physics" begin
    @testset "distance(): metric axioms and invariants" begin
        # 1. Identity of indiscernibles
        @test Erebus.distance(0.0, 0.0, 0.0, 0.0) ≈ 0.0 atol=1e-12
        @test Erebus.distance(42.5, -17.2, 42.5, -17.2) ≈ 0.0 atol=1e-12

        # 2. Symmetry: d(p1, p2) == d(p2, p1)
        p1 = (12.3, -45.6)
        p2 = (-78.9, 10.1)
        d12 = Erebus.distance(p1..., p2...)
        d21 = Erebus.distance(p2..., p1...)
        @test d12 ≈ d21 rtol=1e-12

        # 3. Positivity: d(p1, p2) > 0 for distinct points
        @test d12 > 0.0

        # 4. Triangle inequality: d(A, C) <= d(A, B) + d(B, C)
        pA = (0.0, 0.0)
        pB = (300.0, 400.0)
        pC = (700.0, 100.0)
        dAC = Erebus.distance(pA..., pC...)
        dAB = Erebus.distance(pA..., pB...)
        dBC = Erebus.distance(pB..., pC...)
        @test dAC <= dAB + dBC + 1e-12
        # Collinear points: equality holds
        pMid = (150.0, 200.0)
        @test Erebus.distance(pA..., pB...) ≈
            (Erebus.distance(pA..., pMid...) + Erebus.distance(pMid..., pB...)) rtol=1e-12

        # 5. Translation invariance: d(A + v, B + v) == d(A, B)
        v = (1000.0, -500.0)
        @test Erebus.distance(pA[1] + v[1], pA[2] + v[2], pB[1] + v[1], pB[2] + v[2]) ≈ dAB rtol=1e-12

        # 6. Pinned analytical Pythagorean triples with 3-class discrimination guards
        d34 = Erebus.distance(0.0, 0.0, 3000.0, 4000.0)
        @test isapprox(d34, 5000.0; rtol=1e-12)
        # Exponent / formula error guard: Manhattan L1 distance is 7000.0
        @test abs(d34 - 7000.0) > 1000.0
        # Sign guard
        @test d34 > 0.0
        # Scale guard: order of magnitude in meters
        @test 1000.0 < d34 < 10000.0

        d512 = Erebus.distance(0.0, 0.0, 5000.0, 12000.0)
        @test isapprox(d512, 13000.0; rtol=1e-12)
    end # testset "distance()"

    @testset "total(): two-phase mixture invariants" begin
        s = 3.0e6  # Solid property
        f = 4.2e6  # Fluid property

        # 1. Pure solid limit: phi = 0.0 -> returns solid property
        @test isapprox(Erebus.total(s, f, 0.0), s; rtol=1e-12)

        # 2. Pure fluid limit: phi = 1.0 -> returns fluid property
        @test isapprox(Erebus.total(s, f, 1.0), f; rtol=1e-12)

        # 3. Linear partition of unity midpoint
        @test isapprox(Erebus.total(s, f, 0.5), 0.5 * (s + f); rtol=1e-12)

        # 4. Strict monotonicity across porosity
        val_low = Erebus.total(s, f, 0.2)
        val_high = Erebus.total(s, f, 0.8)
        @test val_low < val_high

        # 5. Physical bounds: min(s, f) <= total <= max(s, f)
        for phi_test in [0.05, 0.25, 0.5, 0.75, 0.95]
            val = Erebus.total(s, f, phi_test)
            @test s <= val <= f
        end

        # 6. Error contracts: unphysical porosities must throw DomainError
        @test_throws DomainError Erebus.total(s, f, -0.1)
        @test_throws DomainError Erebus.total(s, f, 1.1)
        @test_throws DomainError Erebus.total(s, f, NaN)
    end # testset "total()"

    @testset "dot4(): vector space inner product axioms" begin
        # 1. Positivity: dot4(v, v) > 0 for non-zero vectors, == 0 for zero vector
        z = zeros(4)
        @test isapprox(Erebus.dot4(z, z), 0.0; atol=1e-12)
        v = [1.0, -2.0, 3.0, -4.0]
        @test Erebus.dot4(v, v) > 0.0
        @test isapprox(Erebus.dot4(v, v), 1.0 + 4.0 + 9.0 + 16.0; rtol=1e-12)

        # 2. Symmetry: dot4(u, v) == dot4(v, u)
        u = [-5.0, 6.0, 7.0, -8.0]
        @test isapprox(Erebus.dot4(u, v), Erebus.dot4(v, u); rtol=1e-12)

        # 3. Linearity: dot4(a*u + b*w, v) == a*dot4(u, v) + b*dot4(w, v)
        w = [2.0, 1.0, -3.0, 4.0]
        a, b = 2.5, -1.8
        comb = a .* u .+ b .* w
        @test isapprox(
            Erebus.dot4(comb, v), a * Erebus.dot4(u, v) + b * Erebus.dot4(w, v); rtol=1e-12
        )

        # 4. Cauchy-Schwarz inequality: |u . v|^2 <= (u . u) * (v . v)
        dot_uv = Erebus.dot4(u, v)
        @test (dot_uv^2) <= Erebus.dot4(u, u) * Erebus.dot4(v, v) + 1e-12

        # 5. Standard basis orthogonality: e_i . e_j = delta_ij
        for i in 1:4, j in 1:4
            ei = zeros(4)
            ei[i] = 1.0
            ej = zeros(4)
            ej[j] = 1.0
            expected = (i == j) ? 1.0 : 0.0
            @test isapprox(Erebus.dot4(ei, ej), expected; atol=1e-12)
        end
    end # testset "dot4()"

    @testset "ktotal(): physical limits and non-linear mixture invariants" begin
        ks = 3.0  # Matrix rock conductivity [W/m/K]
        kf = 0.6  # Pore water conductivity [W/m/K]

        # 1. Pure solid matrix limit: phi = 0 -> ktotal = ks
        @test isapprox(Erebus.ktotal(ks, kf, 0.0), ks; rtol=1e-10)

        # 2. Pure pore fluid limit: phi = 1 -> ktotal = kf
        @test isapprox(Erebus.ktotal(ks, kf, 1.0), kf; rtol=1e-10)

        # 3. Degenerate homogeneous phase: ks == kf -> ktotal == ks for all phi
        k_homo = 2.5
        @test isapprox(Erebus.ktotal(k_homo, k_homo, 0.1), k_homo; rtol=1e-10)
        @test isapprox(Erebus.ktotal(k_homo, k_homo, 0.5), k_homo; rtol=1e-10)
        @test isapprox(Erebus.ktotal(k_homo, k_homo, 0.9), k_homo; rtol=1e-10)

        # 4. Strict physical monotonicity: ktotal decreases as fluid porosity increases (for ks > kf)
        k_phi20 = Erebus.ktotal(ks, kf, 0.20)
        k_phi40 = Erebus.ktotal(ks, kf, 0.40)
        k_phi60 = Erebus.ktotal(ks, kf, 0.60)
        @test ks > k_phi20 > k_phi40 > k_phi60 > kf

        # 5. Strict bounding: min(ks, kf) <= ktotal <= max(ks, kf)
        for phi in [0.05, 0.2, 0.5, 0.8, 0.95]
            val = Erebus.ktotal(ks, kf, phi)
            @test kf <= val <= ks
        end

        # 6. Pinned reference benchmark with 3-class discrimination guards
        # At ks = 3.0, kf = 0.6, phi = 0.20:
        # k_val = sqrt(3*0.6/2 + ((3*(3*0.2-2) + 0.6*(1-3*0.2))^2)/16) - 0.25*(3*(3*0.2-2) + 0.6*(1-3*0.2))
        # term = 3*(-1.4) + 0.6*(0.4) = -4.2 + 0.24 = -3.96
        # sqrt(0.9 + 3.96^2 / 16) - 0.25*(-3.96) = sqrt(0.9 + 0.9801) + 0.99 = sqrt(1.8801) + 0.99 = 1.371167 + 0.99 = 2.361167
        k_bench = Erebus.ktotal(ks, kf, 0.20)
        expected_bench = sqrt(0.9 + (3.96^2) / 16.0) + 0.99
        @test isapprox(k_bench, expected_bench; rtol=1e-10)
        # Exponent / formula guard: simple linear average is (1-0.2)*3.0 + 0.2*0.6 = 2.520
        @test abs(k_bench - 2.520) > 0.10
        # Sign guard: thermal conductivity must be positive
        @test k_bench > 0.0
        # Scale guard: order of magnitude in W/(m K)
        @test 1.0 < k_bench < 4.0

        # 7. Error contracts: invalid physical domains must throw DomainError
        @test_throws DomainError Erebus.ktotal(-3.0, kf, 0.2)
        @test_throws DomainError Erebus.ktotal(ks, -0.6, 0.2)
        @test_throws DomainError Erebus.ktotal(ks, kf, -0.05)
        @test_throws DomainError Erebus.ktotal(ks, kf, 1.05)
        @test_throws DomainError Erebus.ktotal(ks, kf, NaN)
    end # testset "ktotal()"

    @testset "kphi(): Kozeny-Carman permeability invariants" begin
        k0 = 1.0e-13  # Reference permeability [m^2]
        phi_ref = phim0  # Global reference porosity (0.2)

        # 1. Impermeable limit: as phi -> 0, kphi -> 0
        @test isapprox(Erebus.kphi(k0, 1.0e-6), 0.0; atol=1.0e-25)

        # 2. Reference consistency: at phi == phi_ref, kphi == k0
        @test isapprox(Erebus.kphi(k0, phi_ref), k0; rtol=1e-10)

        # 3. Strict monotonicity: permeability increases strictly monotonically with porosity
        phi_vals = [0.05, 0.10, 0.20, 0.30, 0.50]
        k_vals = [Erebus.kphi(k0, p) for p in phi_vals]
        for i in 1:(length(k_vals) - 1)
            @test k_vals[i] < k_vals[i + 1]
        end

        # 4. Pinned analytical value with 3-class discrimination guards
        # At phi = 0.10, phi_ref = 0.20:
        # k = k0 * (0.1/0.2)^3 * ((1-0.2)/(1-0.1))^2 = k0 * 0.125 * (0.8/0.9)^2 = k0 * 0.125 * (64/81)
        expected_k = k0 * 0.125 * (64.0 / 81.0)
        k_pin = Erebus.kphi(k0, 0.10)
        @test isapprox(k_pin, expected_k; rtol=1e-10)
        # Exponent guard: linear porosity scaling yields 0.5 * k0
        @test abs(k_pin - 0.5 * k0) > 0.3 * k0
        # Sign guard: permeability must be positive
        @test k_pin > 0.0
        # Scale guard: order of magnitude in m^2
        @test 1.0e-16 < k_pin < 1.0e-13

        # 5. Error contracts: invalid porosity and negative reference permeability must throw DomainError
        @test_throws DomainError Erebus.kphi(-1.0e-13, 0.2)
        @test_throws DomainError Erebus.kphi(k0, -0.05)
        @test_throws DomainError Erebus.kphi(k0, 1.0)
        @test_throws DomainError Erebus.kphi(k0, 1.2)
        @test_throws DomainError Erebus.kphi(k0, NaN)
    end # testset "kphi()"

    @testset "ηᶠcur_inv_kᵠ(): Darcy drag resistance invariants" begin
        k0 = 1.0e-13  # Reference permeability [m^2]
        eta = 1.0e-3  # Water viscosity [Pa s]

        # 1. Physical equivalence: ηᶠcur_inv_kᵠ == eta / kphi
        for phi in [0.05, 0.15, 0.25, 0.5]
            expected_drag = eta / Erebus.kphi(k0, phi)
            actual_drag = Erebus.ηᶠcur_inv_kᵠ(k0, phi, eta)
            @test isapprox(actual_drag, expected_drag; rtol=1e-10)
        end

        # 2. Linearity in fluid viscosity: doubling eta doubles drag resistance
        drag_1 = Erebus.ηᶠcur_inv_kᵠ(k0, 0.2, eta)
        drag_2 = Erebus.ηᶠcur_inv_kᵠ(k0, 0.2, 2.0 * eta)
        @test isapprox(drag_2, 2.0 * drag_1; rtol=1e-12)

        # 3. Monotonic decrease with porosity: wider pores reduce flow resistance
        drag_low_phi = Erebus.ηᶠcur_inv_kᵠ(k0, 0.1, eta)
        drag_high_phi = Erebus.ηᶠcur_inv_kᵠ(k0, 0.4, eta)
        @test drag_low_phi > drag_high_phi > 0.0

        # 4. Error contracts: unphysical parameters throw DomainError
        @test_throws DomainError Erebus.ηᶠcur_inv_kᵠ(-k0, 0.2, eta)
        @test_throws DomainError Erebus.ηᶠcur_inv_kᵠ(k0, -0.1, eta)
        @test_throws DomainError Erebus.ηᶠcur_inv_kᵠ(k0, 1.0, eta)
        @test_throws DomainError Erebus.ηᶠcur_inv_kᵠ(k0, 0.2, -eta)
    end # testset "ηᶠcur_inv_kᵠ()"

    @testset "etatotal_rocks(): phase transition and viscosity bounds" begin
        # Sub-solidus temperatures: rock and ice are solid
        t_cold = 200.0  # Below water melting point (273 K) and rock solidus (1416 K)
        for type in 1:2
            eta_cold = Erebus.etatotal_rocks(t_cold, type)
            @test eta_cold >= etamin
            @test eta_cold >= etasolidm[type]
            @test eta_cold >= etafluidm[type]
        end

        # Hydrothermal liquid regime: rock solid, pore water liquid
        t_hydro = 350.0  # Above water melting, below rock solidus
        for type in 1:2
            eta_hydro = Erebus.etatotal_rocks(t_hydro, type)
            @test eta_hydro >= etamin
            @test eta_hydro >= etasolidm[type]
            @test eta_hydro >= etafluidmm[type]
        end

        # Magma regime: above rock solidus
        t_magma = 1600.0
        for type in 1:2
            eta_magma = Erebus.etatotal_rocks(t_magma, type)
            @test eta_magma >= etamin
            @test eta_magma >= etasolidmm[type]
            @test eta_magma >= etafluidmm[type]
        end

        # Error contract: absolute temperature must be positive
        @test_throws DomainError Erebus.etatotal_rocks(-10.0, 1)
        @test_throws DomainError Erebus.etatotal_rocks(0.0, 1)
    end # testset "etatotal_rocks()"

    @testset "Q_radiogenic(): half-life and conservation closure" begin
        f = f_al
        ratio = ratio_al
        E = E_al
        tau = tau_al
        t_half = t_half_al

        # 1. Initial heating power: Q(0) == f * ratio * E / tau
        Q0 = Erebus.Q_radiogenic(f, ratio, E, tau, 0.0)
        expected_Q0 = f * ratio * E * inv(tau)
        @test isapprox(Q0, expected_Q0; rtol=1e-12)

        # 2. Half-life decay precision: at t = t_half, Q(t) == 0.5 * Q0
        Q_half = Erebus.Q_radiogenic(f, ratio, E, tau, t_half)
        @test isapprox(Q_half, 0.5 * Q0; rtol=1e-10)

        # 3. Two half-lives: at t = 2 * t_half, Q(2*t_half) == 0.25 * Q0
        Q_quarter = Erebus.Q_radiogenic(f, ratio, E, tau, 2.0 * t_half)
        @test isapprox(Q_quarter, 0.25 * Q0; rtol=1e-10)

        # 4. Strict monotonic decay: dQ/dt < 0
        @test Q0 > Q_half > Q_quarter > 0.0

        # 5. Integrated total energy closure: \int_0^\infty Q(t) dt = Q0 * tau
        total_energy_per_kg = Q0 * tau
        expected_energy = f * ratio * E
        @test isapprox(total_energy_per_kg, expected_energy; rtol=1e-12)

        # 6. Discrimination guards
        # Prefactor guard: linear decay (Q0 * (1 - t/tau)) differs substantially at 1 half-life
        wrong_linear = Q0 * (1.0 - t_half * inv(tau))
        @test abs(Q_half - wrong_linear) > 0.1 * Q0
        # Sign guard
        @test Q0 > 0.0
        # Scale guard: ~1e-7 W/kg for 26Al at CAI formation
        @test 1.0e-9 < Q0 < 1.0e-5

        # 7. Error contract: negative decay time, non-positive lifetime, negative abundance/ratio throw DomainError
        @test_throws DomainError Erebus.Q_radiogenic(f, ratio, E, tau, -1.0)
        @test_throws DomainError Erebus.Q_radiogenic(f, ratio, E, -tau, 1000.0)
        @test_throws DomainError Erebus.Q_radiogenic(f, ratio, E, 0.0, 1000.0)
        @test_throws DomainError Erebus.Q_radiogenic(-f, ratio, E, tau, 1000.0)
        @test_throws DomainError Erebus.Q_radiogenic(f, -ratio, E, tau, 1000.0)
    end # testset "Q_radiogenic()"

    @testset "calculate_radioactive_heating(): isotope activity and density scaling" begin
        # 1. Inactive isotopes: returns zero heating vectors
        v_zero = @SVector [0.0, 0.0, 0.0]
        heat_solid_off, heat_fluid_off = Erebus.calculate_radioactive_heating(
            false, false, 1000.0
        )
        @test heat_solid_off == v_zero
        @test heat_fluid_off == v_zero

        # 2. Aluminum-26 active: solid phase carries radiogenic power
        heat_s_al, heat_f_al = Erebus.calculate_radioactive_heating(true, false, 0.0)
        @test heat_s_al[1] > 0.0
        @test heat_s_al[2] > 0.0
        @test heat_s_al[3] ≈ 0.0 atol=1e-12  # Sticky air has zero heating
        @test heat_f_al == v_zero

        # Density scaling invariant: Q_vol = Q_mass * rho
        Q_al_mass = Erebus.Q_radiogenic(f_al, ratio_al, E_al, tau_al, 0.0)
        @test isapprox(heat_s_al[1], Q_al_mass * rhosolidm[1]; rtol=1e-10)
        @test isapprox(heat_s_al[2], Q_al_mass * rhosolidm[2]; rtol=1e-10)

        # 3. Iron-60 active: fluid phase carries heating
        heat_s_fe, heat_f_fe = Erebus.calculate_radioactive_heating(false, true, 0.0)
        @test heat_s_fe == v_zero
        @test heat_f_fe[1] > 0.0
        @test heat_f_fe[2] ≈ 0.0 atol=1e-12
        @test heat_f_fe[3] ≈ 0.0 atol=1e-12

        # 4. Temporal decay: heating decreases over time
        heat_s_late, _ = Erebus.calculate_radioactive_heating(true, false, t_half_al)
        @test heat_s_late[1] < heat_s_al[1]
        @test isapprox(heat_s_late[1], 0.5 * heat_s_al[1]; rtol=1e-8)
    end # testset "calculate_radioactive_heating()"

    @testset "compute_gibbs_free_energy(): chemical equilibrium and relaxation" begin
        # 1. Equilibrium limit: when dt >= dtreaction, dGWD == 0.0 identically
        t_hot = 600.0
        pf_val = 5.0e7
        xd = 0.5
        xw = 0.5
        dtr = 1.0e4
        @test iszero(Erebus.compute_gibbs_free_energy(t_hot, pf_val, xd, xw, dtr, dtr))
        @test iszero(
            Erebus.compute_gibbs_free_energy(t_hot, pf_val, xd, xw, 2.0 * dtr, dtr)
        )

        # 2. Incomplete reaction: linear relaxation with timestep ratio (1 - dt/dtr)
        dt_quarter = 0.25 * dtr
        dt_half = 0.50 * dtr
        dt_threequarter = 0.75 * dtr
        dG_quarter = Erebus.compute_gibbs_free_energy(
            t_hot, pf_val, xd, xw, dt_quarter, dtr
        )
        dG_half = Erebus.compute_gibbs_free_energy(t_hot, pf_val, xd, xw, dt_half, dtr)
        dG_threequarter = Erebus.compute_gibbs_free_energy(
            t_hot, pf_val, xd, xw, dt_threequarter, dtr
        )

        @test isapprox(dG_half, (1.0 - 0.50) / (1.0 - 0.25) * dG_quarter; rtol=1e-12)
        @test isapprox(
            dG_threequarter, (1.0 - 0.75) / (1.0 - 0.25) * dG_quarter; rtol=1e-12
        )

        # 3. Chemical equilibrium composition: dG == 0 even for dt = 0 when composition matches equilibrium
        # dG0 = dHWD - T*dSWD + pf*dVWD + RG*T*log(xd/xw) == 0
        ratio_eq = exp(-(ΔHWD - t_hot * ΔSWD + pf_val * ΔVWD) / (RG * t_hot))
        xw_eq = 1.0 / (1.0 + ratio_eq)
        xd_eq = 1.0 - xw_eq
        dG_eq = Erebus.compute_gibbs_free_energy(t_hot, pf_val, xd_eq, xw_eq, 0.0, dtr)
        @test isapprox(dG_eq, 0.0; atol=1e-10)

        # 4. Le Chatelier monotonicity: increasing dry product increases dG (opposes further dehydration)
        dG_low_xd = Erebus.compute_gibbs_free_energy(t_hot, pf_val, 0.2, 0.8, 0.0, dtr)
        dG_high_xd = Erebus.compute_gibbs_free_energy(t_hot, pf_val, 0.8, 0.2, 0.0, dtr)
        @test dG_low_xd < dG_high_xd

        # 5. Discrimination guards on pinned point
        # At T = 600 K, pf = 5e7 Pa, xd = xw = 0.5, log(xd/xw) = 0:
        # dG0 = ΔHWD - 600 * ΔSWD + 5e7 * ΔVWD
        # At dt = 0.5*dtr: dG = dG0 * 0.5
        expected_pinned = 0.5 * (ΔHWD - 600.0 * ΔSWD + 5.0e7 * ΔVWD)
        dG_pinned = Erebus.compute_gibbs_free_energy(
            600.0, 5.0e7, 0.5, 0.5, 5000.0, 10000.0
        )
        @test isapprox(dG_pinned, expected_pinned; rtol=1e-12)
        # Exponent / formula error guard: non-relaxed value is 2 * expected_pinned
        @test abs(dG_pinned - 2.0 * expected_pinned) > 500.0
        # Sign guard: positive driving force opposes spontaneous reaction
        @test dG_pinned > 0.0
        # Scale guard: order of magnitude in J/mol
        @test 500.0 < dG_pinned < 50000.0

        # 6. Error contracts
        @test_throws DomainError Erebus.compute_gibbs_free_energy(
            -100.0, pf_val, xd, xw, 1.0, dtr
        )
        @test_throws DomainError Erebus.compute_gibbs_free_energy(
            t_hot, -1.0, xd, xw, 1.0, dtr
        )
        @test_throws DomainError Erebus.compute_gibbs_free_energy(
            t_hot, pf_val, 0.0, xw, 1.0, dtr
        )
        @test_throws DomainError Erebus.compute_gibbs_free_energy(
            t_hot, pf_val, xd, 0.0, 1.0, dtr
        )
        @test_throws DomainError Erebus.compute_gibbs_free_energy(
            t_hot, pf_val, xd, xw, -1.0, dtr
        )
        @test_throws DomainError Erebus.compute_gibbs_free_energy(
            t_hot, pf_val, xd, xw, 1.0, 0.0
        )
    end # testset "compute_gibbs_free_energy()"

    @testset "compute_relative_enthalpy(): thermodynamic endothermicity and bounds" begin
        # 1. Zero hydrated fraction limit: no hydrated mineral present -> zero reaction enthalpy
        @test isapprox(Erebus.compute_relative_enthalpy(0.8, 0.0), 0.0; atol=1e-12)
        # 2. Zero solid matrix limit: pure pore fluid -> zero reaction enthalpy
        @test isapprox(Erebus.compute_relative_enthalpy(0.0, 0.5), 0.0; atol=1e-12)

        # 3. Endothermic convention: relative enthalpy is non-positive (release of heat during formation)
        for xs in [0.2, 0.5, 0.8], xw in [0.1, 0.5, 0.9]
            @test Erebus.compute_relative_enthalpy(xs, xw) <= 0.0
        end

        # 4. Bilinearity: doubling xs doubles enthalpy magnitude; doubling xw doubles enthalpy magnitude
        h_base = Erebus.compute_relative_enthalpy(0.4, 0.4)
        h_double_xs = Erebus.compute_relative_enthalpy(0.8, 0.4)
        h_double_xw = Erebus.compute_relative_enthalpy(0.4, 0.8)
        @test isapprox(h_double_xs, 2.0 * h_base; rtol=1e-12)
        @test isapprox(h_double_xw, 2.0 * h_base; rtol=1e-12)

        # 5. Pinned benchmark with 3-class discrimination guards
        # At xs = 1.0, xw = 1.0: H = -ΔHWD / (MD + MH₂O)
        # ΔHWD = 74400 J/mol, MD = 0.120 kg/mol, MH₂O = 0.018 kg/mol -> MD + MH₂O = 0.138 kg/mol
        # expected = -74400 / 0.138 = -539130.4348 J/kg
        h_max = Erebus.compute_relative_enthalpy(1.0, 1.0)
        expected_h_max = -ΔHWD / (MD + MH₂O)
        @test isapprox(h_max, expected_h_max; rtol=1e-10)
        # Exponent guard: missing denominator would give -74400
        @test abs(h_max - (-ΔHWD)) > 1.0e5
        # Sign guard: negative enthalpy
        @test h_max < 0.0
        # Scale guard: order of magnitude in J/kg
        @test -1.0e6 < h_max < -1.0e5

        # 6. Numerical tolerance against machine epsilon roundoff
        @test isapprox(Erebus.compute_relative_enthalpy(-1.0e-14, 0.5), 0.0; atol=1e-12)
        @test isapprox(
            Erebus.compute_relative_enthalpy(0.5, 1.0 + 1.0e-14),
            Erebus.compute_relative_enthalpy(0.5, 1.0);
            rtol=1e-12,
        )

        # 7. Error contracts
        @test_throws DomainError Erebus.compute_relative_enthalpy(-0.1, 0.5)
        @test_throws DomainError Erebus.compute_relative_enthalpy(1.1, 0.5)
        @test_throws DomainError Erebus.compute_relative_enthalpy(0.5, -0.1)
        @test_throws DomainError Erebus.compute_relative_enthalpy(0.5, 1.1)
    end # testset "compute_relative_enthalpy()"

    @testset "compute_reaction_constant(): Van 't Hoff equilibrium and monotonicity" begin
        # 1. Strict positivity: K > 0 for all physical inputs
        @test Erebus.compute_reaction_constant(500.0, 1.0e7, 0.0) > 0.0
        @test Erebus.compute_reaction_constant(800.0, 5.0e7, 1.0e4) > 0.0

        # 2. Van 't Hoff equilibrium identity: when driving force equals dG, K == 1.0 identically
        t_val = 650.0
        pf_val = 3.0e7
        dG_drive = ΔHWD - t_val * ΔSWD + ΔVWD * pf_val
        @test isapprox(
            Erebus.compute_reaction_constant(t_val, pf_val, dG_drive), 1.0; rtol=1e-12
        )

        # 3. Monotonicity in chemical affinity: increasing dG promotes forward reaction (increases K)
        k_low = Erebus.compute_reaction_constant(t_val, pf_val, 0.0)
        k_high = Erebus.compute_reaction_constant(t_val, pf_val, 1.0e4)
        @test k_low < k_high

        # 4. Van 't Hoff temperature dependence for endothermic reaction (ΔHWD > 0):
        # Equilibrium constant increases with temperature when dG = 0
        k_cold = Erebus.compute_reaction_constant(500.0, 1.0e7, 0.0)
        k_hot = Erebus.compute_reaction_constant(700.0, 1.0e7, 0.0)
        @test k_hot > k_cold

        # 5. Pinned analytical calculation with 3-class discrimination guards
        # At T = 500 K, pf = 0, dG = ΔHWD - 500*ΔSWD: exponent is 0 -> K = 1.0
        @test isapprox(
            Erebus.compute_reaction_constant(500.0, 0.0, ΔHWD - 500.0 * ΔSWD),
            1.0;
            rtol=1e-12,
        )
        # Sign guard: positive constant
        @test k_hot > 0.0
        # Scale guard: physical equilibrium constant
        @test 1.0e-30 < k_cold < 1.0e10

        # 6. Error contracts
        @test_throws DomainError Erebus.compute_reaction_constant(-300.0, 1.0e7, 0.0)
        @test_throws DomainError Erebus.compute_reaction_constant(500.0, -1.0e6, 0.0)
    end # testset "compute_reaction_constant()"

    @testset "compute_rhocpfluidm(): phase transitions and latent heat spike" begin
        # 1. Positivity across all physical phases
        for t in [150.0, 250.0, 273.0, 300.0, 400.0, 500.0]
            @test Erebus.compute_rhocpfluidm(t, 1) > 0.0
            @test Erebus.compute_rhocpfluidm(t, 9) > 0.0
        end

        # 2. Apparent heat capacity method: latent heat spike in [tmfluidphase - 5, tmfluidphase + 5]
        # At melting point, apparent volumetric heat capacity must exceed both sub-freezing and liquid plateau values
        rhocp_ice = Erebus.compute_rhocpfluidm(tmfluidphase - 10.0, 1)
        rhocp_melt_solid = Erebus.compute_rhocpfluidm(tmfluidphase - 2.0, 1)
        rhocp_melt_liquid = Erebus.compute_rhocpfluidm(tmfluidphase + 2.0, 1)
        rhocp_liquid = Erebus.compute_rhocpfluidm(tmfluidphase + 10.0, 1)

        @test rhocp_melt_solid > rhocp_ice
        @test rhocp_melt_liquid > rhocp_liquid
        # Latent heat contribution: 0.1 * Lf * rho
        @test rhocp_melt_liquid - rhocp_liquid ≈ ρH₂Oᶠ * 0.1 * Lᶠ rtol=1e-10

        # 3. Liquid water plateau: heat capacity is constant 4200 J/(kg K) * rho between (tmfluidphase + 5) and 410 K
        expected_liquid_plateau = ρH₂Oᶠ * 4200.0
        @test isapprox(
            Erebus.compute_rhocpfluidm(290.0, 1), expected_liquid_plateau; rtol=1e-12
        )
        @test isapprox(
            Erebus.compute_rhocpfluidm(350.0, 1), expected_liquid_plateau; rtol=1e-12
        )
        @test isapprox(
            Erebus.compute_rhocpfluidm(400.0, 1), expected_liquid_plateau; rtol=1e-12
        )

        # 4. Mode 9: constant property independent of temperature
        const_val = rhocpfluidm[1]
        for t in [100.0, 273.0, 500.0]
            @test Erebus.compute_rhocpfluidm(t, 9) ≈ const_val
        end

        # 5. Discrimination guards
        # Sign guard
        @test rhocp_liquid > 0.0
        # Scale guard: ~1e6 to 1e8 J/(m³ K)
        @test 1.0e6 < rhocp_liquid < 1.0e8

        # 6. Error contracts
        @test_throws DomainError Erebus.compute_rhocpfluidm(-50.0, 1)
        @test_throws DomainError Erebus.compute_rhocpfluidm(0.0, 1)
        @test_throws ArgumentError Erebus.compute_rhocpfluidm(300.0, 2)
        @test_throws ArgumentError Erebus.compute_rhocpfluidm(300.0, -1)
    end # testset "compute_rhocpfluidm()"

    @testset "compute_ksolidm(): phonon scattering limits and asymptotics" begin
        # 1. Strict monotonic decrease with temperature: Umklapp phonon scattering reduces lattice conductivity
        ts = [150.0, 250.0, 350.0, 500.0, 800.0, 1200.0]
        ks_vals = [Erebus.compute_ksolidm(t, 1) for t in ts]
        for i in 1:(length(ks_vals) - 1)
            @test ks_vals[i] > ks_vals[i + 1]
        end

        # 2. High-temperature asymptotic floor: ks -> 0.73 W/(m K) as T -> infty
        k_very_hot = Erebus.compute_ksolidm(1.0e8, 1)
        @test isapprox(k_very_hot, 0.73; atol=1e-4)
        @test k_very_hot > 0.73

        # 3. Room temperature silicate reference: T = 300 K -> ks ~ 4.16 W/(m K) (typical olivine/dunite)
        k_300 = Erebus.compute_ksolidm(300.0, 1)
        expected_300 = 0.73 + 1293.0 / (300.0 + 77.0)
        @test isapprox(k_300, expected_300; rtol=1e-12)
        @test 3.5 < k_300 < 5.0

        # 4. Mode 9: constant property independent of temperature
        for t in [200.0, 500.0, 1000.0]
            @test Erebus.compute_ksolidm(t, 9) ≈ ksolidm[1]
        end

        # 5. Discrimination guards
        # Sign guard
        @test k_300 > 0.0
        # Exponent guard: linear temperature dependence differs strongly
        @test abs(k_300 - (0.73 + 1293.0 / 300.0)) > 0.5

        # 6. Error contracts
        @test_throws DomainError Erebus.compute_ksolidm(-10.0, 1)
        @test_throws DomainError Erebus.compute_ksolidm(0.0, 1)
        @test_throws ArgumentError Erebus.compute_ksolidm(300.0, 2)
    end # testset "compute_ksolidm()"

    @testset "compute_kfluidm(): phase regime bounds and water conductivity" begin
        # 1. Positivity across all physical phases
        for t in [150.0, 250.0, 273.0, 320.0, 400.0, 600.0]
            @test Erebus.compute_kfluidm(t, 1) > 0.0
            @test Erebus.compute_kfluidm(t, 9) > 0.0
        end

        # 2. Ice phase monotonicity: thermal conductivity of crystalline ice decreases with temperature
        k_ice_cold = Erebus.compute_kfluidm(150.0, 1)
        k_ice_warm = Erebus.compute_kfluidm(250.0, 1)
        @test k_ice_cold > k_ice_warm

        # 3. Liquid water thermal conductivity magnitude: between 0.5 and 0.8 W/(m K)
        k_water_room = Erebus.compute_kfluidm(293.15, 1)
        k_water_hot = Erebus.compute_kfluidm(373.15, 1)
        @test 0.5 < k_water_room < 0.8
        @test 0.5 < k_water_hot < 0.8

        # 4. Mode 9: constant property independent of temperature
        for t in [100.0, 300.0, 700.0]
            @test Erebus.compute_kfluidm(t, 9) ≈ kfluidm[1]
        end

        # 5. Discrimination guards
        # Sign guard
        @test k_water_room > 0.0
        # Scale guard: water thermal conductivity is well within [0.1, 10.0] W/(m K)
        @test 0.1 < k_water_room < 10.0

        # 6. Error contracts
        @test_throws DomainError Erebus.compute_kfluidm(-50.0, 1)
        @test_throws DomainError Erebus.compute_kfluidm(0.0, 1)
        @test_throws ArgumentError Erebus.compute_kfluidm(300.0, 5)
    end # testset "compute_kfluidm()"

    @testset "compute_Δtreaction(): kinetics and activation timescale scaling" begin
        # 1. Porosity scaling invariant: reaction timescale is inversely proportional to porosity
        # (higher porosity -> more fluid-mineral interface -> faster completion)
        t_mid = 500.0
        phi_a = 0.1
        phi_b = 0.2
        for m in 1:3
            dtr_a = Erebus.compute_Δtreaction(t_mid, phi_a, m)
            dtr_b = Erebus.compute_Δtreaction(t_mid, phi_b, m)
            @test isapprox(dtr_a, 2.0 * dtr_b; rtol=1e-12)
            @test dtr_a > 0.0
            @test dtr_b > 0.0
        end

        # 2. Mode 1 (Gaussian rate): optimal kinetics (minimum reaction time) at T = c_I
        t_opt = c_I
        dtr_opt = Erebus.compute_Δtreaction(t_opt, 0.1, 1)
        dtr_below = Erebus.compute_Δtreaction(t_opt - 50.0, 0.1, 1)
        dtr_above = Erebus.compute_Δtreaction(t_opt + 50.0, 0.1, 1)
        @test dtr_opt < dtr_below
        @test dtr_opt < dtr_above
        # Symmetric Gaussian response
        @test isapprox(dtr_below, dtr_above; rtol=1e-12)

        # 3. Mode 2 & 3 (Arrhenius rates): monotonic decrease of reaction time with temperature
        for m in 2:3
            dtr_cold = Erebus.compute_Δtreaction(350.0, 0.1, m)
            dtr_warm = Erebus.compute_Δtreaction(450.0, 0.1, m)
            dtr_hot = Erebus.compute_Δtreaction(550.0, 0.1, m)
            @test dtr_cold > dtr_warm > dtr_hot > 0.0
        end

        # 4. Mode 9: constant property independent of temperature and porosity
        @test Erebus.compute_Δtreaction(200.0, 0.05, 9) ≈ Δtreaction
        @test Erebus.compute_Δtreaction(600.0, 0.50, 9) ≈ Δtreaction

        # 5. Discrimination guards
        # Sign guard
        @test dtr_opt > 0.0
        # Scale guard: reaction time is positive
        @test dtr_opt > 1.0

        # 6. Error contracts
        @test_throws DomainError Erebus.compute_Δtreaction(-100.0, 0.1, 1)
        @test_throws DomainError Erebus.compute_Δtreaction(500.0, -0.1, 1)
        @test_throws DomainError Erebus.compute_Δtreaction(500.0, 0.0, 1)
        @test_throws DomainError Erebus.compute_Δtreaction(500.0, 1.5, 1)
        @test_throws ArgumentError Erebus.compute_Δtreaction(500.0, 0.1, 4)
    end # testset "compute_Δtreaction()"

    @testset "perform_thermochemical_reaction!(): phase conservation and bounds" begin
        marknum = start_marknum
        DMP = zeros(Ny1, Nx1)
        DHP = zeros(Ny1, Nx1)
        DMPSUM = zeros(Ny1, Nx1)
        DHPSUM = zeros(Ny1, Nx1)
        WTPSUM = zeros(Ny1, Nx1)
        pf = fill(3.0e7, Ny1, Nx1)
        tk2 = fill(650.0, Ny1, Nx1)

        # Marker setup: mantle rock (tm=1), crust rock (tm=2), and sticky air (tm=3)
        tm = rand(rgen, 1:3, marknum)
        xm = rand(rgen, 0.0:0.1:xsize, marknum)
        ym = rand(rgen, 0.0:0.1:ysize, marknum)
        XWˢm₀ = fill(0.5, marknum)
        XWˢm = copy(XWˢm₀)
        phim = fill(0.1, marknum)
        phinewm = copy(phim)
        pfm₀ = fill(3.0e7, marknum)
        dt = 1.0e4

        Erebus.perform_thermochemical_reaction!(
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
            XWˢm₀,
            XWˢm,
            phim,
            phinewm,
            pfm₀,
            marknum,
            dt,
            1,
            3,
        )

        # 1. Porosity and composition bounds: physical intervals [0, 1]
        for m in 1:marknum
            if tm[m] < 3
                @test 0.0 <= phinewm[m] <= 1.0
                @test 0.0 <= XWˢm[m] <= 1.0
                @test pfm₀[m] >= 0.0
            end
        end

        # 2. Partition of unity: at every node where markers contributed, WTPSUM > 0
        for j in 1:Nx1, i in 1:Ny1
            if WTPSUM[i, j] > 0.0
                @test isfinite(DMP[i, j])
                @test isfinite(DHP[i, j])
            else
                @test isapprox(DMP[i, j], 0.0; atol=1e-12)
                @test isapprox(DHP[i, j], 0.0; atol=1e-12)
            end
        end

        # 3. Latent heat term sign: DHP > 0 strictly pins enthalpy transfer sign
        @test any(WTPSUM .> 0.0)
        @test any(DHP .> 0.0)
        @test !any(DHP .< 0.0)

        # 4. Timestep evolution: for timestep > 1, XWˢm₀ and phim remain at pre-reaction values
        XW_prev = copy(XWˢm₀)
        phi_prev = copy(phim)
        Erebus.perform_thermochemical_reaction!(
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
            XWˢm₀,
            XWˢm,
            phim,
            phinewm,
            pfm₀,
            marknum,
            dt,
            2,
            3,
        )
        reacted = [m for m in 1:marknum if tm[m] < 3 && phinewm[m] != phi_prev[m]]
        @test !isempty(reacted)
        for m in reacted
            # Pre-reaction arrays preserved on timestep > 1
            @test XWˢm₀[m] ≈ XW_prev[m]
            @test phim[m] ≈ phi_prev[m]
            # Reacted arrays updated
            @test phinewm[m] != phim[m]
        end

        # 4. Sticky air markers (tm = 3) do not react
        # Create a pure air marker run
        marknum_air = 10
        tm_air = fill(3, marknum_air)
        xm_air = fill(xsize / 2, marknum_air)
        ym_air = fill(ysize / 2, marknum_air)
        XW_air0 = fill(0.0, marknum_air)
        XW_air = copy(XW_air0)
        phi_air = fill(phimin, marknum_air)
        phinew_air = copy(phi_air)
        pf_air = fill(0.0, marknum_air)
        DMP_air = zeros(Ny1, Nx1)
        DHP_air = zeros(Ny1, Nx1)
        DMPSUM_air = zeros(Ny1, Nx1)
        DHPSUM_air = zeros(Ny1, Nx1)
        WTPSUM_air = zeros(Ny1, Nx1)

        Erebus.perform_thermochemical_reaction!(
            DMP_air,
            DHP_air,
            DMPSUM_air,
            DHPSUM_air,
            WTPSUM_air,
            pf,
            tk2,
            tm_air,
            xm_air,
            ym_air,
            XW_air0,
            XW_air,
            phi_air,
            phinew_air,
            pf_air,
            marknum_air,
            dt,
            2,
            1,
        )
        @test all(iszero, DMP_air)
        @test all(iszero, DHP_air)
        @test all(iszero, WTPSUM_air)
    end # testset "perform_thermochemical_reaction!()"

    @testset "compute_shear_heating!(): Second Law non-negativity and symmetry" begin
        HS = zeros(Ny1, Nx1)
        ETA = rand(rgen, Ny, Nx) .+ 1.0e20
        SXY = rand(rgen, Ny, Nx) .* 1.0e6
        ETAP = rand(rgen, Ny1, Nx1) .+ 1.0e20
        SXX = rand(rgen, Ny1, Nx1) .* 1.0e6
        RX = rand(rgen, Ny1, Nx1) .+ 1.0e10
        RY = rand(rgen, Ny1, Nx1) .+ 1.0e10
        qxD = rand(rgen, Ny1, Nx1) .* 1.0e-7
        qyD = rand(rgen, Ny1, Nx1) .* 1.0e-7
        PHI = fill(0.1, Ny1, Nx1)
        ETAPHI = rand(rgen, Ny1, Nx1) .+ 1.0e20
        pr = rand(rgen, Ny1, Nx1) .* 1.0e7
        pf = rand(rgen, Ny1, Nx1) .* 1.0e7

        Erebus.compute_shear_heating!(
            HS, ETA, SXY, ETAP, SXX, RX, RY, qxD, qyD, PHI, ETAPHI, pr, pf
        )

        # 1. Second Law of Thermodynamics: viscous dissipation is strictly non-negative everywhere
        for j in 2:Nx, i in 2:Ny
            @test HS[i, j] >= 0.0
        end

        # 2. Quiescent limit: zero stress and zero Darcy flux -> identically zero dissipation
        HS_zero = zeros(Ny1, Nx1)
        Erebus.compute_shear_heating!(
            HS_zero,
            ETA,
            zeros(Ny, Nx),
            ETAP,
            zeros(Ny1, Nx1),
            RX,
            RY,
            zeros(Ny1, Nx1),
            zeros(Ny1, Nx1),
            PHI,
            ETAPHI,
            pr,
            pr,  # pr == pf
        )
        for j in 2:Nx, i in 2:Ny
            @test isapprox(HS_zero[i, j], 0.0; atol=1e-12)
        end

        # 3. Quadratic symmetry: dissipation is invariant under sign inversion of stresses and fluxes
        HS_neg = zeros(Ny1, Nx1)
        Erebus.compute_shear_heating!(
            HS_neg, ETA, -SXY, ETAP, -SXX, RX, RY, -qxD, -qyD, PHI, ETAPHI, pf, pr
        )
        for j in 2:Nx, i in 2:Ny
            @test isapprox(HS[i, j], HS_neg[i, j]; rtol=1e-12)
        end

        # 4. Strict positivity with non-zero shear stress
        @test any(HS[2:Ny, 2:Nx] .> 0.0)
    end # testset "compute_shear_heating!()"

    @testset "compute_adiabatic_heating!(): thermodynamic work and linearity" begin
        HA = zeros(Ny1, Nx1)
        tk1 = fill(300.0, Ny1, Nx1)
        ALPHA = fill(3.0e-5, Ny1, Nx1)
        ALPHAF = fill(2.0e-4, Ny1, Nx1)
        PHI = fill(0.1, Ny1, Nx1)
        vx = fill(1.0e-8, Ny1, Nx1)
        vy = fill(1.0e-8, Ny1, Nx1)
        vxf = fill(1.0e-8, Ny1, Nx1)
        vyf = fill(1.0e-8, Ny1, Nx1)

        # 1. Isobaric limit: spatially uniform pressures -> identically zero adiabatic heating
        ps_const = fill(5.0e7, Ny1, Nx1)
        pf_const = fill(3.0e7, Ny1, Nx1)
        Erebus.compute_adiabatic_heating!(
            HA, tk1, ALPHA, ALPHAF, PHI, vx, vy, vxf, vyf, ps_const, pf_const
        )
        for j in 2:Nx, i in 2:Ny
            @test isapprox(HA[i, j], 0.0; atol=1e-12)
        end

        # 2. Incompressible limit: zero thermal expansion -> identically zero adiabatic heating
        ps_grad = [
            5.0e7 + 1.0e4 * j + 5.0e2 * j^2 + 3.0e3 * i + 2.0e2 * i^2 for
            i in 1:Ny1, j in 1:Nx1
        ]
        pf_grad = [
            3.0e7 + 1.0e4 * j + 5.0e2 * j^2 + 3.0e3 * i + 2.0e2 * i^2 for
            i in 1:Ny1, j in 1:Nx1
        ]
        Erebus.compute_adiabatic_heating!(
            HA,
            tk1,
            zeros(Ny1, Nx1),
            zeros(Ny1, Nx1),
            PHI,
            vx,
            vy,
            vxf,
            vyf,
            ps_grad,
            pf_grad,
        )
        for j in 2:Nx, i in 2:Ny
            @test isapprox(HA[i, j], 0.0; atol=1e-12)
        end

        # 3. Upwind stencil discrimination on non-linear pressure field
        # For quadratic field, forward difference strictly differs from backward difference:
        #   (ps[i, j+1] - ps[i, j])/dx != (ps[i, j] - ps[i, j-1])/dx
        # Positive velocity chooses forward stencil for solid, backward stencil for fluid
        HA_pos = zeros(Ny1, Nx1)
        Erebus.compute_adiabatic_heating!(
            HA_pos, tk1, ALPHA, ALPHAF, PHI, vx, vy, vxf, vyf, ps_grad, pf_grad
        )
        # Verify exact forward/backward difference evaluation at interior point (Ny÷2, Nx÷2)
        i_mid, j_mid = Ny ÷ 2, Nx ÷ 2
        dpsdx_fwd = (ps_grad[i_mid, j_mid + 1] - ps_grad[i_mid, j_mid]) / dx
        dpsdx_bwd = (ps_grad[i_mid, j_mid] - ps_grad[i_mid, j_mid - 1]) / dx
        @test abs(dpsdx_fwd - dpsdx_bwd) > 1.0e-3  # Second derivative non-vanishing
        @test HA_pos[i_mid, j_mid] > 0.0

        # Negative velocity chooses backward stencil for solid, forward stencil for fluid
        HA_neg = zeros(Ny1, Nx1)
        Erebus.compute_adiabatic_heating!(
            HA_neg, tk1, ALPHA, ALPHAF, PHI, -vx, -vy, -vxf, -vyf, ps_grad, pf_grad
        )
        for j in 2:Nx, i in 2:Ny
            @test HA_neg[i, j] < 0.0
            # Under non-linear pressure gradient, magnitude of forward and backward steps differs
            @test !isapprox(abs(HA_pos[i, j]), abs(HA_neg[i, j]); rtol=1e-4)
        end

        # 4. Linearity in thermal expansion: doubling expansivity doubles heating
        HA_double = zeros(Ny1, Nx1)
        Erebus.compute_adiabatic_heating!(
            HA_double,
            tk1,
            2.0 .* ALPHA,
            2.0 .* ALPHAF,
            PHI,
            vx,
            vy,
            vxf,
            vyf,
            ps_grad,
            pf_grad,
        )
        for j in 2:Nx, i in 2:Ny
            @test isapprox(HA_double[i, j], 2.0 * HA_pos[i, j]; rtol=1e-12)
        end
    end # testset "compute_adiabatic_heating!()"

    @testset "poroelastic constitutive relations" begin
        # Incompressible baseline limits
        betaphi = 1.0e-11
        phi = 0.1
        bd0 = Erebus.compute_drained_compressibility(betaphi, phi, 0.0)
        @test bd0 ≈ betaphi / (1.0 - phi) rtol=1e-12
        @test Erebus.compute_biot_willis_coefficient(bd0, 0.0) ≈ 1.0 rtol=1e-12
        @test Erebus.compute_skempton_coefficient(bd0, phi, 0.0, 0.0) ≈ 1.0 rtol=1e-12

        # Physical rock and fluid values
        betasolid = 2.5e-11 # silicate matrix [1/Pa]
        betafluid = 4.0e-10 # liquid water [1/Pa]
        bd = Erebus.compute_drained_compressibility(betaphi, phi, betasolid)
        @test bd ≈ (betaphi + betasolid) / (1.0 - phi) rtol=1e-12
        @test bd > betasolid

        kbw = Erebus.compute_biot_willis_coefficient(bd, betasolid)
        @test 0.0 < kbw < 1.0
        @test kbw ≈ 1.0 - betasolid / bd rtol=1e-12

        ksk = Erebus.compute_skempton_coefficient(bd, phi, betasolid, betafluid)
        @test 0.0 < ksk < 1.0
        num = bd - betasolid
        denom = num + phi * (betafluid - betasolid)
        @test ksk ≈ num / denom rtol=1e-12

        # Parameter sweep across full production porosity range and first-step betaphi = 0.0
        for p in [phimin, 0.01, 0.1, 0.5, phimax], bp in [0.0, 1.0e-11, 1.0e-10]
            bd_sweep = Erebus.compute_drained_compressibility(bp, p, betasolid)
            @test bd_sweep >= betasolid
            kbw_sweep = Erebus.compute_biot_willis_coefficient(bd_sweep, betasolid)
            @test 0.0 <= kbw_sweep <= 1.0
            ksk_sweep = Erebus.compute_skempton_coefficient(
                bd_sweep, p, betasolid, betafluid
            )
            @test 0.0 <= ksk_sweep <= 1.0
        end
    end # testset "poroelastic constitutive relations"

    @testset "thermal buoyancy and compute_rhofluid()" begin
        rho0 = 1000.0
        alpha = 2.0e-4
        T0 = 273.15

        # Sub-freezing / baseline: no thermal expansion
        @test Erebus.compute_rhofluid(250.0, rho0, alpha, T0) ≈ rho0
        @test Erebus.compute_rhofluid(T0, rho0, alpha, T0) ≈ rho0

        # When thermal_buoyancy = false, always returns rho0
        @test Erebus.compute_rhofluid(500.0, rho0, alpha, T0; thermal_buoyancy=false) ≈ rho0

        # Positive thermal expansion for T > T0
        T1 = 373.15 # ΔT = 100 K
        expected_rho1 = rho0 * (1.0 - alpha * 100.0) # 1000 * (1 - 0.02) = 980.0
        @test Erebus.compute_rhofluid(T1, rho0, alpha, T0) ≈ expected_rho1 rtol=1e-12
        @test Erebus.compute_rhofluid(T1, rho0, alpha, T0) < rho0

        # Extreme temperature clamping (clamped to 0.1 * rho0)
        T_extreme = 10000.0
        @test Erebus.compute_rhofluid(T_extreme, rho0, alpha, T0) ≈ 0.1 * rho0

        # Zero or negative alpha
        @test Erebus.compute_rhofluid(400.0, rho0, 0.0, T0) ≈ rho0
        @test Erebus.compute_rhofluid(400.0, rho0, -1e-4, T0) ≈ rho0

        # Verification of thermal buoyancy driving force: Δρ = -ρ0 * α * ΔT
        g = 9.81
        delta_T = 50.0 # K
        rho_hot = Erebus.compute_rhofluid(T0 + delta_T, rho0, alpha, T0)
        buoyancy_force = (rho0 - rho_hot) * g
        expected_buoyancy = rho0 * alpha * delta_T * g
        @test buoyancy_force ≈ expected_buoyancy rtol=1e-12
    end

    @testset "temperature-dependent fluid viscosity and compute_fluid_viscosity()" begin
        # Sub-freezing ice viscosity
        @test Erebus.compute_fluid_viscosity(200.0, 1) ≈ 1.0e12 rtol=1e-12
        @test Erebus.compute_fluid_viscosity(273.0, 2) ≈ 1.0e12 rtol=1e-12

        # Sticky air
        @test Erebus.compute_fluid_viscosity(300.0, 3) ≈ 1.0e-3 rtol=1e-12
        @test Erebus.compute_fluid_viscosity(100.0, 3) ≈ 1.0e-3 rtol=1e-12

        # Reference temperature T0 = 293.15 K -> eta0 = 1.0e-3 Pa s
        @test Erebus.compute_fluid_viscosity(293.15, 2) ≈ 1.0e-3 rtol=1e-12

        # Boiling water T = 373.15 K -> drops by factor ~3.7
        eta_373 = Erebus.compute_fluid_viscosity(
            373.15, 2; Ea=15.0e3, T0=293.15, eta0=1.0e-3
        )
        @test 2.5e-4 < eta_373 < 3.0e-4

        # Supercritical hydrothermal fluid T = 600 K -> drops below 1e-4 Pa s
        eta_600 = Erebus.compute_fluid_viscosity(
            600.0, 2; Ea=15.0e3, T0=293.15, eta0=1.0e-3
        )
        @test 3.0e-5 < eta_600 < 1.0e-4
        @test eta_600 < eta_373 < 1.0e-3

        # Mode :constant
        @test Erebus.compute_fluid_viscosity(500.0, 2; mode=:constant) ≈ 1.0e-3 rtol=1e-12

        # Clamp against extreme temperatures
        @test Erebus.compute_fluid_viscosity(10000.0, 2; etamin=1.0e-5) ≈ 1.0e-5 rtol=1e-12

        # NaN / non-finite handling (safe fallback to ice viscosity)
        @test Erebus.compute_fluid_viscosity(NaN, 2) ≈ 1.0e12 rtol=1e-12
        @test Erebus.compute_fluid_viscosity(Inf, 2) ≈ 1.0e12 rtol=1e-12

        # Invalid mode throws ArgumentError
        @test_throws ArgumentError Erebus.compute_fluid_viscosity(300.0, 2; mode=:bogus)
    end

    @testset "dynamic hydrofracturing and compute_hydrofracture_permeability()" begin
        k0 = 1.0e-15 # m²
        sigma_t = 1.0e7 # 10 MPa
        kappa_frac = 1.0e3

        # Case 1: Compressive / sub-tensile regime (Peff = 5 MPa > -10 MPa) -> no hydrofracture
        Peff_comp = 5.0e6
        @test Erebus.compute_hydrofracture_factor(
            Peff_comp, sigma_t; kappa_frac=kappa_frac
        ) ≈ 1.0 rtol=1e-12
        @test Erebus.compute_hydrofracture_permeability(
            k0, Peff_comp, sigma_t; kappa_frac=kappa_frac
        ) ≈ k0 rtol=1e-12

        # Case 2: Exact tensile limit (Peff = -10 MPa) -> overpressure = 0 -> factor = 1.0
        Peff_limit = -1.0e7
        @test Erebus.compute_hydrofracture_factor(
            Peff_limit, sigma_t; kappa_frac=kappa_frac
        ) ≈ 1.0 rtol=1e-12
        @test Erebus.compute_hydrofracture_permeability(
            k0, Peff_limit, sigma_t; kappa_frac=kappa_frac
        ) ≈ k0 rtol=1e-12

        # Case 3: Overpressured regime (Peff = -20 MPa, overpressure = 10 MPa = 1.0 * sigma_t)
        # factor = 1.0 + 1000.0 * (1.0)^1.0 = 1001.0
        Peff_over = -2.0e7
        expected_factor = 1.0 + kappa_frac * 1.0
        @test Erebus.compute_hydrofracture_factor(
            Peff_over, sigma_t; kappa_frac=kappa_frac
        ) ≈ expected_factor rtol=1e-12
        @test Erebus.compute_hydrofracture_permeability(
            k0, Peff_over, sigma_t; kappa_frac=kappa_frac
        ) ≈ k0 * expected_factor rtol=1e-12

        # Case 4: Power-law scaling with gamma = 2.0
        # normalized overpressure = 2.0 (Peff = -30 MPa)
        # factor = 1.0 + 1000.0 * (2.0)^2.0 = 4001.0
        Peff_over2 = -3.0e7
        @test Erebus.compute_hydrofracture_factor(
            Peff_over2, sigma_t; kappa_frac=kappa_frac, gamma=2.0
        ) ≈ 4001.0 rtol=1e-12

        # Case 5: Inactive (active=false)
        @test Erebus.compute_hydrofracture_factor(Peff_over, sigma_t; active=false) ≈ 1.0 rtol=1e-12
        @test Erebus.compute_hydrofracture_permeability(
            k0, Peff_over, sigma_t; active=false
        ) ≈ k0 rtol=1e-12

        # Case 6: Permeability ceiling (kmax)
        kmax = 1.0e-13
        @test Erebus.compute_hydrofracture_permeability(
            k0, -1.0e8, sigma_t; kappa_frac=kappa_frac, kmax=kmax
        ) ≈ kmax rtol=1e-12

        # Case 7: Non-finite / corrupt inputs
        @test Erebus.compute_hydrofracture_factor(NaN, sigma_t) ≈ 1.0 rtol=1e-12
        @test Erebus.compute_hydrofracture_factor(Peff_over, NaN) ≈ 1.0 rtol=1e-12
        @test Erebus.compute_hydrofracture_factor(Peff_over, -1.0e7) ≈ 1.0 rtol=1e-12
        @test Erebus.compute_hydrofracture_permeability(k0, NaN, sigma_t) ≈ k0 rtol=1e-12
    end
end
