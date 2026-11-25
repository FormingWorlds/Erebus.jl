module MeltingShellValidation

# Melting shell case: the porous shell of a planetesimal (dry porosity laws of Henke et
# al. 2012, 26Al heating) with silicate melting and latent heat, conduction only. The
# melt has the density of the solid (no melt buoyancy) and soft turbulence is off.

using Erebus
using StaticArrays

isdefined(@__MODULE__, :PorousShellValidation) ||
    include(joinpath(@__DIR__, "validation_porous_shell.jl"))

const PSV = PorousShellValidation

const L_MELT = 4.0e5
const T_SOL = 1400.0
const T_LIQ = 1800.0

"""
    melt_fraction_ref(T; T_sol=T_SOL, T_liq=T_LIQ) -> Float64

Melt fraction of rock at temperature `T` [K]: linear between `T_sol` and `T_liq`, clamped
to [0, 1] (no pressure shift).
"""
function melt_fraction_ref(T; T_sol=T_SOL, T_liq=T_LIQ)
    return clamp((T - T_sol) / (T_liq - T_sol), 0.0, 1.0)
end

"""
    rhocp_app_ref(T, phi; L_melt=L_MELT, T_sol=T_SOL, T_liq=T_LIQ) -> Float64

Apparent volumetric heat capacity [J/(m^3 K)] of rock with dry porosity `phi` at `T` [K]:
`(1 - phi) (RHOCP_SOL + RHO_SOL L_melt / (T_liq - T_sol))` for `T_sol < T < T_liq`, the
sensible term alone outside.
"""
function rhocp_app_ref(T, phi; L_melt=L_MELT, T_sol=T_SOL, T_liq=T_LIQ)
    c = PSV.RHOCP_SOL + (T_sol < T < T_liq ? PSV.RHO_SOL * L_melt / (T_liq - T_sol) : 0.0)
    return (1.0 - phi) * c
end

"""
    enthalpy_s_ref(T; L_melt=L_MELT, T_sol=T_SOL, T_liq=T_LIQ, T0=PSV.T_SPACE) -> Float64

Enthalpy per unit volume of dense rock [J/m^3] from `T0` to `T`: sensible heat plus
`RHO_SOL L_melt F(T)`.
"""
function enthalpy_s_ref(T; L_melt=L_MELT, T_sol=T_SOL, T_liq=T_LIQ, T0=PSV.T_SPACE)
    F = melt_fraction_ref(T; T_sol=T_sol, T_liq=T_liq)
    return PSV.RHOCP_SOL * (T - T0) + PSV.RHO_SOL * L_melt * F
end

"""
    make_melting_shell_config(; Nx, phi_init, t_form_yr, n_steps, dt_yr,
                              L_melt=L_MELT, T_sol=T_SOL, T_liq=T_LIQ, kwargs...)

Configuration of the porous shell case (`make_porous_shell_config`, which takes the other
keywords) with `melting.active`, latent heat `L_melt` [J/kg] over the rock melting
interval `T_sol` to `T_liq` [K], the melt density equal to the rock density and soft
turbulence off.

# Throws

    - `ArgumentError`: if `T_sol >= T_liq`, or the errors of `make_porous_shell_config`
"""
function make_melting_shell_config(;
    L_melt::Float64=L_MELT, T_sol::Float64=T_SOL, T_liq::Float64=T_LIQ, kwargs...
)
    T_sol < T_liq || throw(ArgumentError("T_sol must be below T_liq, got $T_sol, $T_liq"))
    melting = MeltingConfig(;
        active=true,
        soft_turbulence=false,
        L_melt=L_melt,
        T_solidus=SVector{3,Float64}([T_sol, T_sol, NaN]),
        T_liquidus=SVector{3,Float64}([T_liq, T_liq, NaN]),
        rho_melt=PSV.RHO_SOL,
    )
    return PSV.make_porous_shell_config(; melting=melting, kwargs...)
end

"""
    melt_radius(xm, ym, T_max; xc, yc, R_p, dr, T_sol=T_SOL, nsector=16) -> Float64

Radius [m] of the molten centre from the peak temperatures `T_max` [K] of rock markers at
`(xm, ym)`. The markers with `r < R_p` about `(xc, yc)` fall into `nsector` equal angular
sectors and centre-anchored radial bins `[(b - 1) dr, min(b dr, R_p))`. In each sector the
scan runs outward over the non-empty bins while the bin median of `T_max` is above
`T_sol`; the sector radius is the outer edge of the last such bin (0 if the first non-empty
bin is not molten). The result is the median of the sector radii.

# Throws

    - `ArgumentError`: if `xm`, `ym` and `T_max` differ in length, or no marker is inside `R_p`
    - `DomainError`: if `dr` or `R_p` is not positive and finite, or `nsector < 1`
"""
function melt_radius(xm, ym, T_max; xc, yc, R_p, dr, T_sol=T_SOL, nsector::Int=16)
    length(xm) == length(ym) == length(T_max) ||
        throw(ArgumentError("xm, ym and T_max must have the same length"))
    (isfinite(dr) && dr > 0.0) || throw(DomainError(dr, "dr must be positive and finite"))
    (isfinite(R_p) && R_p > 0.0) ||
        throw(DomainError(R_p, "R_p must be positive and finite"))
    nsector >= 1 || throw(DomainError(nsector, "nsector must be at least 1"))
    nb = ceil(Int, R_p / dr)
    bins = [[Float64[] for _ in 1:nb] for _ in 1:nsector]
    for m in eachindex(xm)
        dx, dy = xm[m] - xc, ym[m] - yc
        r = hypot(dx, dy)
        r < R_p || continue
        s = min(floor(Int, mod(atan(dy, dx), 2pi) / (2pi / nsector)) + 1, nsector)
        push!(bins[s][min(floor(Int, r / dr) + 1, nb)], T_max[m])
    end
    any(sec -> any(!isempty, sec), bins) || throw(ArgumentError("no marker is inside R_p"))
    radii = Float64[]
    for sec in bins
        any(!isempty, sec) || continue
        rs = 0.0
        for b in 1:nb
            isempty(sec[b]) && continue
            PSV._median(sec[b]) > T_sol || break
            rs = min(b * dr, R_p)
        end
        push!(radii, rs)
    end
    return PSV._median(radii)
end

end # module MeltingShellValidation
