
"""
Set up geodesic and physical properties of the set of markers.

$(SIGNATURES)

# Details

    - randomized: fill in random values for grid properties instead of zeros

# Returns

    - xm : horizontal marker coordinate [m]
    - ym : vertical marker coordinate [m]
    - tm : marker material type
    - tkm : marker temperature [K]
    - sxxm : marker σ′xx [Pa]
    - sxym : marker σxy [Pa]
    - etavpm : marker viscoplastic viscosity [Pa]
    - phim : marker porosity
"""
function setup_marker_properties(
    marknum; randomized=false, coords::Union{Nothing,GridCoordinates}=nothing
)
    dx_val = coords === nothing ? dx : coords.dx
    dy_val = coords === nothing ? dy : coords.dy
    xsize_val = coords === nothing ? xsize : coords.xsize
    ysize_val = coords === nothing ? ysize : coords.ysize
    # horizontal marker coordinate [m]
    xm = if randomized
        rand(rgen, (-dx_val):0.1:(xsize_val + dx_val), marknum)
    else
        zeros(marknum)
    end
    # vertical marker coordinate [m]
    ym = if randomized
        rand(rgen, (-dy_val):0.1:(ysize_val + dy_val), marknum)
    else
        zeros(marknum)
    end
    # marker material type
    tm = randomized ? rand(rgen, 1:3, marknum) : zeros(Int, marknum)
    # marker temperature [K]
    tkm = randomized ? rand(rgen, 273:300, marknum) : zeros(marknum)
    # marker σ′xx [Pa]
    sxxm = randomized ? rand(rgen, marknum) : zeros(marknum)
    # marker σxy [Pa]
    sxym = randomized ? rand(rgen, marknum) : zeros(marknum)
    # marker viscoplastic viscosity [Pa]
    etavpm = randomized ? rand(rgen, marknum) : zeros(marknum)
    # marker porosity
    phim = randomized ? rand(rgen, marknum) : zeros(marknum)
    # reacted marker porosity
    phinewm = randomized ? rand(rgen, marknum) : zeros(marknum)
    # previous marker fluid pressure
    pfm0 = randomized ? rand(rgen, marknum) : zeros(marknum)
    # marker melt molar fraction
    XWsolidm = randomized ? rand(rgen, marknum) : zeros(marknum)
    # previous marker melt molar fraction
    XWsolidm0 = randomized ? rand(rgen, marknum) : zeros(marknum)
    return (xm, ym, tm, tkm, sxxm, sxym, etavpm, phim, phinewm, pfm0, XWsolidm, XWsolidm0)
end # function setup_marker_properties()

function setup_marker_properties(marknum, coords::GridCoordinates; randomized=false)
    return setup_marker_properties(marknum; randomized=randomized, coords=coords)
end

"""
Set up additional helper marker properties to facility comptuations.

$(SIGNATURES)

# Details

    - randomized: fill in random values for grid properties instead of zeros

# Returns

    - rhototalm: total density of markers
    - rhocptotalm : total volumetric heat capacity of markers
    - etatotalm: total viscosity of markers
    - hrtotalm: total radiogenic heat production of markers
    - ktotalm: total thermal conductivity of markers
    - tkm_rhocptotalm: total thermal energy of markers
    - etafluidcur_inv_kphim: fluid viscosity over permeability of markers
    - inv_gggtotalm: inverse of total shear modulus of markers
    - fricttotalm: total friction coefficient of markers
    - cohestotalm: total compressive strength of markers
    - tenstotalm: total tensile strength of markers
    - rhofluidcur: fluid density of markers
    - alphasolidcur: solid thermal expansion coefficient of markers
    - alphafluidcur: fluid thermal expansion coefficient of markers
"""
function setup_marker_properties_helpers(marknum; randomized=false)
    # marker total density
    rhototalm = randomized ? rand(rgen, marknum) : zeros(marknum)
    # marker total volumetric heat capacity
    rhocptotalm = randomized ? rand(rgen, marknum) : zeros(marknum)
    # marker total viscosity
    etatotalm = randomized ? rand(rgen, marknum) : zeros(marknum)
    # marker total radiogenic heat production
    hrtotalm = randomized ? rand(rgen, marknum) : zeros(marknum)
    # marker total thermal conductivity
    ktotalm = randomized ? rand(rgen, marknum) : zeros(marknum)
    # marker total thermal energy
    tkm_rhocptotalm = randomized ? rand(rgen, marknum) : zeros(marknum)
    # marker fluid viscosity over permeability
    etafluidcur_inv_kphim = randomized ? rand(rgen, marknum) : zeros(marknum)
    # marker inverse of total shear modulus
    inv_gggtotalm = randomized ? rand(rgen, marknum) : zeros(marknum)
    # marker total friction coefficient
    fricttotalm = randomized ? rand(rgen, marknum) : zeros(marknum)
    # marker total compressive strength
    cohestotalm = randomized ? rand(rgen, marknum) : zeros(marknum)
    # marker total tensile strength
    tenstotalm = randomized ? rand(rgen, marknum) : zeros(marknum)
    # marker fluid density
    rhofluidcur = randomized ? rand(rgen, marknum) : zeros(marknum)
    # marker solid thermal expansion coefficient
    alphasolidcur = randomized ? rand(rgen, marknum) : zeros(marknum)
    # marker fluid thermal expansion coefficient
    alphafluidcur = randomized ? rand(rgen, marknum) : zeros(marknum)
    return (
        rhototalm,
        rhocptotalm,
        etatotalm,
        hrtotalm,
        ktotalm,
        tkm_rhocptotalm,
        etafluidcur_inv_kphim,
        inv_gggtotalm,
        fricttotalm,
        cohestotalm,
        tenstotalm,
        rhofluidcur,
        alphasolidcur,
        alphafluidcur,
    )
end # function setup_marker_properties_helpers()

"""
Set up additional marker geometry helpers to facilitate marker handling.

$(SIGNATURES)

# Details

    - nothing

# Returns

    - mdis: minimum distance of marker launch anchor points to nearest marker
    - mnum: number of marker nearest to marker launch anchor positions
"""
function setup_marker_geometry_helpers(Nxm::Int=Nxm, Nym::Int=Nym)
    mdis = fill(mdis_init, Nym, Nxm)
    mnum = zeros(Int, Nym, Nxm)
    return mdis, mnum
end
function setup_marker_geometry_helpers(coords::GridCoordinates)
    return setup_marker_geometry_helpers(coords.Nxm, coords.Nym)
end

"""
Define initial set of markers according to model parameters

$(SIGNATURES)

# Details

    - xm: x coordinates of markers
    - ym: y coordinates of markers
    - tm: material type of markers
    - phim: porosity of markers
    - etavpm: matrix viscosity of markers
    - rhototalm: total density of markers
    - rhocptotalm: total volumetric heat capacity of markers
    - etatotalm: total viscosity of markers
    - hrtotalm: total radiogenic heat production of markers
    - ktotalm: total thermal conductivity of markers
    - tkm: temperature of markers 
    - inv_gggtotalm: inverse of total shear modulus of markers
    - fricttotalm: total friction coefficient of markers
    - cohestotalm: total compressive strength of markers
    - tenstotalm: total tensile strength of markers
    - rhofluidcur: fluid density of markers
    - alphasolidcur: solid thermal expansion coefficient of markers
    - alphafluidcur: fluid thermal expansion coefficient of markers
    - XWsolidm0: previous wet solid molar fraction of markers
    - randomized: uniformly random-distribute marker x/y positions within cells
                  and randomly set initial marker porosity 

# Returns

    - nothing
"""
function define_markers!(
    xm,
    ym,
    tm,
    phim,
    etavpm,
    rhototalm,
    rhocptotalm,
    etatotalm,
    hrtotalm,
    ktotalm,
    tkm,
    inv_gggtotalm,
    fricttotalm,
    cohestotalm,
    tenstotalm,
    rhofluidcur,
    alphasolidcur,
    alphafluidcur,
    XWsolidm0;
    randomized=random_markers,
    coords=nothing,
    xcenter_val=xcenter,
    ycenter_val=ycenter,
    rplanet_val=rplanet,
    rcrust_val=rcrust,
)
    Nxm_val = coords === nothing ? Nxm : coords.Nxm
    Nym_val = coords === nothing ? Nym : coords.Nym
    dxm_val = coords === nothing ? dxm : coords.dxm
    dym_val = coords === nothing ? dym : coords.dym
    for jm in 1:1:Nxm_val, im in 1:1:Nym_val
        # calculate marker counter
        m = (jm-1) * Nym_val + im
        # define marker coordinates
        xm[m] = dxm_val/2 + (jm-1) * dxm_val
        ym[m] = dym_val/2 + (im-1) * dym_val
        # random marker position within cell
        if randomized
            xm[m] += (rand(rgen)-0.5) * dxm_val
            ym[m] += (rand(rgen)-0.5) * dym_val
        end
        # primary marker properties 
        rmark = distance(xm[m], ym[m], xcenter_val, ycenter_val)
        if rmark < rplanet_val
            # planet
            tm[m] = ifelse(rmark>rcrust_val, 2, 1)
            # porosity
            phim[m] = phim0
            if randomized
                phim[m] += phim0 * (rand(rgen)-0.5)
            end
            # matrix viscosity
            etavpm[m] = etasolidm[tm[m]] # * exp(-αη*phim[m]) # ∇! CHANGE!!!
            # wet solid molar fraction
            XWsolidm0[m] = XWsolidm_init[tm[m]]
            if randomized
                XWsolidm0[m] += XWsolidm_init[tm[m]] * (rand(rgen)-0.5)
            end
        else
            # sticky space ("air") [to have internal free surface]
            tm[m] = 3
            # porosity
            phim[m] = phimin
            # matrix viscosity
            etavpm[m] = etasolidm[tm[m]]
            # static properties for air markers
            rhototalm[m] = rhosolidm[tm[m]]
            rhocptotalm[m] = rhocpsolidm[tm[m]]
            etatotalm[m] = etasolidm[tm[m]]
            hrtotalm[m] = start_hrsolidm[tm[m]]
            ktotalm[m] = ksolidm[tm[m]]
        end
        # common initialisations for all marker types
        tkm[m] = tkm0[tm[m]]
        inv_gggtotalm[m] = inv(gggsolidm[tm[m]])
        fricttotalm[m] = frictsolidm[tm[m]]
        cohestotalm[m] = cohessolidm[tm[m]]
        tenstotalm[m] = tenssolidm[tm[m]]
        rhofluidcur[m] = rhofluidm[tm[m]]
        alphasolidcur[m] = alphasolidm[tm[m]]
        alphafluidcur[m] = alphafluidm[tm[m]]
    end
    return nothing
end

"""
Compute properties of given marker and save them to corresponding arrays.

$(SIGNATURES)

# Details

    - m: marker number
    - tm: type of markers
    - tkm: temperature of markers
    - rhototalm: total density of markers
    - rhocptotalm: total volumetric heat capacity of markers
    - etatotalm: total viscosity of markers
    - hrtotalm: total radiogenic heat production of markers
    - ktotalm: total thermal conductivity of markers
    - tkm_rhocptotalm: total thermal energy of markers
    - etafluidcur_inv_kphim: (fluid viscosity)/permeability of markers
    - hrsolidm: vector of radiogenic heat production of solid materials
    - hrfluidm: vector of radiogenic heat production of fluid materials
    - phim: porosity of markers
    - XWˢm₀: previous wet solid molar fraction of markers
    - mode: marker property computation mode
        - 1: dynamic, based on (Touloukian, 1970; Hobbs, 1974;
             Travis and Schubert, 2005)
        - 2: constant parameter rhocpfluidm

# Returns

    - nothing
"""
function compute_marker_properties!(
    m,
    tm,
    tkm,
    rhototalm,
    rhocptotalm,
    etatotalm,
    hrtotalm,
    ktotalm,
    tkm_rhocptotalm,
    etafluidcur_inv_kphim,
    hrsolidm,
    hrfluidm,
    phim,
    XWˢm₀,
    mode,
    rhofluidcur=nothing;
    thermal_buoyancy::Bool=true,
    alphafluid=nothing,
    tmfluidphase_val::Real=tmfluidphase,
    fluid_viscosity_mode::Symbol=:arrhenius,
    fluid_viscosity_Ea::Real=15.0e3,
    fluid_viscosity_T0::Real=293.15,
    fluid_viscosity_eta0::Real=1.0e-3,
)
    # @timeit to "compute_marker_properties!" begin
    if tm[m] < 3
        # rocks
        XDˢm₀ = 1.0 - XWˢm₀[m]
        rhosolidm0 = (MD + MH₂O*XWˢm₀[m]) / (VDˢ*XDˢm₀ + VWˢ*XWˢm₀[m]) # (16.161)
        alpha_val = if alphafluid === nothing
            alphafluidm[tm[m]]
        else
            (alphafluid isa Real ? alphafluid : alphafluid[tm[m]])
        end
        rhofluidm0 = ifelse(
            tkm[m] > tmfluidphase_val,
            compute_rhofluid(
                tkm[m],
                ρH₂Oᶠ,
                alpha_val,
                tmfluidphase_val;
                thermal_buoyancy=thermal_buoyancy,
            ),
            ρH₂Oᶠⁱ,
        ) # (16.162)
        if rhofluidcur !== nothing
            rhofluidcur[m] = rhofluidm0
        end
        rhototalm[m] = total(rhosolidm0, rhofluidm0, phim[m])
        rhocptotalm[m] = total(
            rhocpsolidm[tm[m]], compute_rhocpfluidm(tkm[m], mode), phim[m]
        )
        etasolidcur = ifelse(tkm[m]>tmsolidphase, etasolidmm[tm[m]], etasolidm[tm[m]])
        etafluidcur = compute_fluid_viscosity(
            tkm[m],
            tm[m];
            mode=fluid_viscosity_mode,
            eta0=fluid_viscosity_eta0,
            eta_ice=etafluidm[tm[m]],
            eta_air=etafluidm[tm[m]],
            Ea=fluid_viscosity_Ea,
            T0=fluid_viscosity_T0,
            tmfluidphase=tmfluidphase_val,
        )
        etatotalm[m] = max(etamin, etasolidcur, etafluidcur)
        hrtotalm[m] = total(hrsolidm[tm[m]], hrfluidm[tm[m]], phim[m])
        ktotalm[m] = ktotal(
            compute_ksolidm(tkm[m], mode), compute_kfluidm(tkm[m], mode), phim[m]
        )
    else
        # sticky air
        etafluidcur = etafluidm[tm[m]]
        if rhofluidcur !== nothing
            rhofluidcur[m] = rhofluidm[tm[m]]
        end
    end
    # common for rocks and air
    tkm_rhocptotalm[m] = tkm[m] * rhocptotalm[m]
    # kphim[m] = kphi(kphim0[tm[m]], phim[m])
    # etafluidcur_inv_kphim[m] = etafluidcur[m] * inv(kphim[m])
    etafluidcur_inv_kphim[m] = ηᶠcur_inv_kᵠ(kphim0[tm[m]], phim[m], etafluidcur)
    # end # @timeit to "compute_marker_properties!"
    return nothing
end # function compute_marker_properties!

"""
Update marker viscoplastic viscosity based on basic node yielding status and
marker temperature- and material-based viscosity.

$(SIGNATURES)

# Details

    - m: marker number
    - xm: x-coordinates of markers
    - ym: y-coordinates of markers
    - tm: type of markers
    - tkm: temperature of markers
    - etatotalm: total viscosity of markers
    - etavpm: matrix viscosity of markers
    - YNY: plastic yielding node property at basic nodes
    - YNY_inv_ETA: inverse viscoplastic viscosity at yielding basic nodes

# Returns

    -nothing
"""
function update_marker_viscosity!(
    m,
    xm,
    ym,
    tm,
    tkm,
    etatotalm,
    etavpm,
    YNY,
    YNY_inv_ETA;
    coords::Union{Nothing,GridCoordinates}=nothing,
)
    x_val = coords === nothing ? x : coords.x
    y_val = coords === nothing ? y : coords.y
    dx_val = coords === nothing ? dx : coords.dx
    dy_val = coords === nothing ? dy : coords.dy
    jmin_basic_val = coords === nothing ? jmin_basic : coords.jmin_basic
    jmax_basic_val = coords === nothing ? jmax_basic : coords.jmax_basic
    imin_basic_val = coords === nothing ? imin_basic : coords.imin_basic
    imax_basic_val = coords === nothing ? imax_basic : coords.imax_basic
    @inbounds i, j, weights = fix_weights(
        xm[m],
        ym[m],
        x_val,
        y_val,
        dx_val,
        dy_val,
        jmin_basic_val,
        jmax_basic_val,
        imin_basic_val,
        imax_basic_val,
    )
    @inbounds if tm[m] < 3
        # rocks: update etatotalm[m] based on current marker temperature
        @inbounds etatotalm[m] = etatotal_rocks(tkm[m], tm[m]) # * exp(-αη*phim[m]) # ∇! CHANGE!!!
        # else
        # air: constant etatotalm[m]=etasolidm[tm[m]] as initialized
        # pass
    end
    if any(grid_vector(i, j, YNY))
        interpolate_to_marker!(m, i, j, weights, etavpm, YNY_inv_ETA)
        @inbounds etavpm[m] = inv(etavpm[m])
        @inbounds etavpm[m] = ifelse(etavpm[m]>etatotalm[m], etatotalm[m], etavpm[m])
    else
        @inbounds etavpm[m] = etatotalm[m]
    end
    return nothing
end

"""
Set up properties to be interpolated from markers to staggered grid.

$(SIGNATURES)

# Details

    - nothing

# Returns

    - ETA0SUM: interpolation of ETA0 at basic nodes
    - ETASUM: interpolation of ETA at basic nodes
    - GGGSUM: interpolation of GGG at basic nodes
    - SXYSUM: interpolation of SXY at basic nodes
    - COHSUM: interpolation of COH at basic nodes
    - TENSUM: interpolation of TEN at basic nodes
    - FRISUM: interpolation of FRI at basic nodes
    - WTSUM: interpolation weights at basic nodes
    - RHOXSUM: interpolation of RHOX at Vx nodes
    - RHOFXSUM: interpolation of RHOFX at Vx nodes
    - KXSUM: interpolation of KX at Vx nodes
    - PHIXSUM: interpolation of PHIX at Vx nodes
    - RXSUM: interpolation of RX at Vx nodes
    - WTXSUM: interpolation weights at Vx nodes
    - RHOYSUM: interpolation of RHOY at Vy nodes
    - RHOFYSUM: interpolation of RHOFY at Vy nodes
    - KYSUM: interpolation of KY at Vy nodes
    - PHIYSUM: interpolation of PHIX at Vy nodes
    - RYSUM: interpolation of RY at Vy nodes
    - WTYSUM: interpolation weights at Vy nodes
    - RHOSUM: interpolation of RHO at P nodes
    - RHOCPSUM: interpolation of RHOCP at P nodes
    - ALPHASUM: interpolation of ALPHA at P nodes
    - ALPHAFSUM: interpolation of ALPHAF at P nodes
    - HRSUM: interpolation of HR at P nodes
    - GGGPSUM: interpolation of GGGP at P nodes
    - SXXSUM: interpolation of SXX at P nodes
    - TKSUM: interpolation of TK at P nodes
    - PHISUM: interpolation of PHI at P nodes
    - DMPSUM: interpolation of DMP at P nodes
    - DHPSUM: interpolation of DHP at P nodes
    - XWSSUM: interpolation of XWS at P nodes
    - WTPSUM: interpolation weights at P nodes
"""
function setup_interpolated_properties(Nx::Int=Nx, Ny::Int=Ny)
    Nx1 = Nx + 1
    Ny1 = Ny + 1
    # basic nodes
    ETA0SUM = zeros(Ny, Nx)
    ETASUM = zeros(Ny, Nx)
    GGGSUM = zeros(Ny, Nx)
    SXYSUM = zeros(Ny, Nx)
    COHSUM = zeros(Ny, Nx)
    TENSUM = zeros(Ny, Nx)
    FRISUM = zeros(Ny, Nx)
    WTSUM = zeros(Ny, Nx)
    # Vx nodes
    RHOXSUM = zeros(Ny1, Nx1)
    RHOFXSUM = zeros(Ny1, Nx1)
    KXSUM = zeros(Ny1, Nx1)
    PHIXSUM = zeros(Ny1, Nx1)
    RXSUM = zeros(Ny1, Nx1)
    WTXSUM = zeros(Ny1, Nx1)
    # Vy nodes
    RHOYSUM = zeros(Ny1, Nx1)
    RHOFYSUM = zeros(Ny1, Nx1)
    KYSUM = zeros(Ny1, Nx1)
    PHIYSUM = zeros(Ny1, Nx1)
    RYSUM = zeros(Ny1, Nx1)
    WTYSUM = zeros(Ny1, Nx1)
    # P Nodes
    RHOSUM = zeros(Ny1, Nx1)
    RHOCPSUM = zeros(Ny1, Nx1)
    ALPHASUM = zeros(Ny1, Nx1)
    ALPHAFSUM = zeros(Ny1, Nx1)
    HRSUM = zeros(Ny1, Nx1)
    GGGPSUM = zeros(Ny1, Nx1)
    SXXSUM = zeros(Ny1, Nx1)
    TKSUM = zeros(Ny1, Nx1)
    PHISUM = zeros(Ny1, Nx1)
    DMPSUM = zeros(Ny1, Nx1)
    DHPSUM = zeros(Ny1, Nx1)
    XWSSUM = zeros(Ny1, Nx1)
    WTPSUM = zeros(Ny1, Nx1)
    return (
        ETA0SUM,
        ETASUM,
        GGGSUM,
        SXYSUM,
        COHSUM,
        TENSUM,
        FRISUM,
        WTSUM,
        RHOXSUM,
        RHOFXSUM,
        KXSUM,
        PHIXSUM,
        RXSUM,
        WTXSUM,
        RHOYSUM,
        RHOFYSUM,
        KYSUM,
        PHIYSUM,
        RYSUM,
        WTYSUM,
        RHOSUM,
        RHOCPSUM,
        ALPHASUM,
        ALPHAFSUM,
        HRSUM,
        GGGPSUM,
        SXXSUM,
        TKSUM,
        PHISUM,
        DMPSUM,
        DHPSUM,
        XWSSUM,
        WTPSUM,
    )
end
function setup_interpolated_properties(coords::GridCoordinates)
    return setup_interpolated_properties(coords.Nx, coords.Ny)
end

"""
Per-thread buffers for thread-safe parallel marker interpolation.
"""
struct ThreadInterpolationBuffers
    ETA0SUM::Matrix{Float64}
    ETASUM::Matrix{Float64}
    GGGSUM::Matrix{Float64}
    SXYSUM::Matrix{Float64}
    COHSUM::Matrix{Float64}
    TENSUM::Matrix{Float64}
    FRISUM::Matrix{Float64}
    WTSUM::Matrix{Float64}
    RHOXSUM::Matrix{Float64}
    RHOFXSUM::Matrix{Float64}
    KXSUM::Matrix{Float64}
    PHIXSUM::Matrix{Float64}
    RXSUM::Matrix{Float64}
    WTXSUM::Matrix{Float64}
    RHOYSUM::Matrix{Float64}
    RHOFYSUM::Matrix{Float64}
    KYSUM::Matrix{Float64}
    PHIYSUM::Matrix{Float64}
    RYSUM::Matrix{Float64}
    WTYSUM::Matrix{Float64}
    RHOSUM::Matrix{Float64}
    RHOCPSUM::Matrix{Float64}
    ALPHASUM::Matrix{Float64}
    ALPHAFSUM::Matrix{Float64}
    HRSUM::Matrix{Float64}
    GGGPSUM::Matrix{Float64}
    SXXSUM::Matrix{Float64}
    TKSUM::Matrix{Float64}
    PHISUM::Matrix{Float64}
    WTPSUM::Matrix{Float64}
end

"""
Allocate a vector of per-thread interpolation buffers.
"""
function allocate_thread_interpolation_buffers(
    nthreads::Int=max(Threads.nthreads(), Threads.maxthreadid()), Nx::Int=Nx, Ny::Int=Ny
)
    Nx1 = Nx + 1
    Ny1 = Ny + 1
    return [
        ThreadInterpolationBuffers(
            zeros(Ny, Nx),
            zeros(Ny, Nx),
            zeros(Ny, Nx),
            zeros(Ny, Nx),
            zeros(Ny, Nx),
            zeros(Ny, Nx),
            zeros(Ny, Nx),
            zeros(Ny, Nx),
            zeros(Ny1, Nx1),
            zeros(Ny1, Nx1),
            zeros(Ny1, Nx1),
            zeros(Ny1, Nx1),
            zeros(Ny1, Nx1),
            zeros(Ny1, Nx1),
            zeros(Ny1, Nx1),
            zeros(Ny1, Nx1),
            zeros(Ny1, Nx1),
            zeros(Ny1, Nx1),
            zeros(Ny1, Nx1),
            zeros(Ny1, Nx1),
            zeros(Ny1, Nx1),
            zeros(Ny1, Nx1),
            zeros(Ny1, Nx1),
            zeros(Ny1, Nx1),
            zeros(Ny1, Nx1),
            zeros(Ny1, Nx1),
            zeros(Ny1, Nx1),
            zeros(Ny1, Nx1),
            zeros(Ny1, Nx1),
            zeros(Ny1, Nx1),
        ) for _ in 1:nthreads
    ]
end
function allocate_thread_interpolation_buffers(nthreads::Int, coords::GridCoordinates)
    return allocate_thread_interpolation_buffers(nthreads, coords.Nx, coords.Ny)
end
function allocate_thread_interpolation_buffers(coords::GridCoordinates)
    return allocate_thread_interpolation_buffers(
        max(Threads.nthreads(), Threads.maxthreadid()), coords.Nx, coords.Ny
    )
end

"""
Reset all per-thread interpolation buffers to zero.
"""
function reset_thread_buffers!(buffers::Vector{ThreadInterpolationBuffers})
    Threads.@threads for b in buffers
        b.ETA0SUM .= 0.0
        b.ETASUM .= 0.0
        b.GGGSUM .= 0.0
        b.SXYSUM .= 0.0
        b.COHSUM .= 0.0
        b.TENSUM .= 0.0
        b.FRISUM .= 0.0
        b.WTSUM .= 0.0
        b.RHOXSUM .= 0.0
        b.RHOFXSUM .= 0.0
        b.KXSUM .= 0.0
        b.PHIXSUM .= 0.0
        b.RXSUM .= 0.0
        b.WTXSUM .= 0.0
        b.RHOYSUM .= 0.0
        b.RHOFYSUM .= 0.0
        b.KYSUM .= 0.0
        b.PHIYSUM .= 0.0
        b.RYSUM .= 0.0
        b.WTYSUM .= 0.0
        b.RHOSUM .= 0.0
        b.RHOCPSUM .= 0.0
        b.ALPHASUM .= 0.0
        b.ALPHAFSUM .= 0.0
        b.HRSUM .= 0.0
        b.GGGPSUM .= 0.0
        b.SXXSUM .= 0.0
        b.TKSUM .= 0.0
        b.PHISUM .= 0.0
        b.WTPSUM .= 0.0
    end
    return nothing
end

"""
Reduce (accumulate) per-thread interpolation buffers into the master grid arrays.
"""
function reduce_thread_buffers!(
    ETA0SUM,
    ETASUM,
    GGGSUM,
    SXYSUM,
    COHSUM,
    TENSUM,
    FRISUM,
    WTSUM,
    RHOXSUM,
    RHOFXSUM,
    KXSUM,
    PHIXSUM,
    RXSUM,
    WTXSUM,
    RHOYSUM,
    RHOFYSUM,
    KYSUM,
    PHIYSUM,
    RYSUM,
    WTYSUM,
    RHOSUM,
    RHOCPSUM,
    ALPHASUM,
    ALPHAFSUM,
    HRSUM,
    GGGPSUM,
    SXXSUM,
    TKSUM,
    PHISUM,
    WTPSUM,
    buffers::Vector{ThreadInterpolationBuffers},
)
    for b in buffers
        ETA0SUM .+= b.ETA0SUM
        ETASUM .+= b.ETASUM
        GGGSUM .+= b.GGGSUM
        SXYSUM .+= b.SXYSUM
        COHSUM .+= b.COHSUM
        TENSUM .+= b.TENSUM
        FRISUM .+= b.FRISUM
        WTSUM .+= b.WTSUM
        RHOXSUM .+= b.RHOXSUM
        RHOFXSUM .+= b.RHOFXSUM
        KXSUM .+= b.KXSUM
        PHIXSUM .+= b.PHIXSUM
        RXSUM .+= b.RXSUM
        WTXSUM .+= b.WTXSUM
        RHOYSUM .+= b.RHOYSUM
        RHOFYSUM .+= b.RHOFYSUM
        KYSUM .+= b.KYSUM
        PHIYSUM .+= b.PHIYSUM
        RYSUM .+= b.RYSUM
        WTYSUM .+= b.WTYSUM
        RHOSUM .+= b.RHOSUM
        RHOCPSUM .+= b.RHOCPSUM
        ALPHASUM .+= b.ALPHASUM
        ALPHAFSUM .+= b.ALPHAFSUM
        HRSUM .+= b.HRSUM
        GGGPSUM .+= b.GGGPSUM
        SXXSUM .+= b.SXXSUM
        TKSUM .+= b.TKSUM
        PHISUM .+= b.PHISUM
        WTPSUM .+= b.WTPSUM
    end
    return nothing
end

"""
Reset properties to be interpolated from markers to staggered grid.

$(SIGNATURES)

# Details   

    - ETA0SUM: interpolation of ETA0 at basic nodes
    - ETASUM: interpolation of ETA at basic nodes
    - GGGSUM: interpolation of GGG at basic nodes
    - SXYSUM: interpolation of SXY at basic nodes
    - COHSUM: interpolation of COH at basic nodes
    - TENSUM: interpolation of TEN at basic nodes
    - FRISUM: interpolation of FRI at basic nodes
    - WTSUM: interpolation weights at basic nodes
    - RHOXSUM: interpolation of RHOX at Vx nodes
    - RHOFXSUM: interpolation of RHOFX at Vx nodes
    - KXSUM: interpolation of KX at Vx nodes
    - PHIXSUM: interpolation of PHIX at Vx nodes
    - RXSUM: interpolation of RX at Vx nodes
    - WTXSUM: interpolation weights at Vx nodes
    - RHOYSUM: interpolation of RHOY at Vy nodes
    - RHOFYSUM: interpolation of RHOFY at Vy nodes
    - KYSUM: interpolation of KY at Vy nodes
    - PHIYSUM: interpolation of PHIX at Vy nodes
    - RYSUM: interpolation of RY at Vy nodes
    - WTYSUM: interpolation weights at Vy nodes
    - RHOSUM: interpolation of RHO at P nodes
    - RHOCPSUM: interpolation of RHOCP at P nodes
    - ALPHASUM: interpolation of ALPHA at P nodes
    - ALPHAFSUM: interpolation of ALPHAF at P nodes
    - HRSUM: interpolation of HR at P nodes
    - GGGPSUM: interpolation of GGGP at P nodes
    - SXXSUM: interpolation of SXX at P nodes
    - TKSUM: interpolation of TK at P nodes
    - PHISUM: interpolation of PHI at P nodes
    - WTPSUM: interpolation weights at P nodes

# Returns

    - nothing
"""
function reset_interpolated_properties!(
    ETA0SUM,
    ETASUM,
    GGGSUM,
    SXYSUM,
    COHSUM,
    TENSUM,
    FRISUM,
    WTSUM,
    RHOXSUM,
    RHOFXSUM,
    KXSUM,
    PHIXSUM,
    RXSUM,
    WTXSUM,
    RHOYSUM,
    RHOFYSUM,
    KYSUM,
    PHIYSUM,
    RYSUM,
    WTYSUM,
    RHOSUM,
    RHOCPSUM,
    ALPHASUM,
    ALPHAFSUM,
    HRSUM,
    GGGPSUM,
    SXXSUM,
    TKSUM,
    PHISUM,
    WTPSUM,
)
    # basic nodes
    ETA0SUM .= zero(0.0)
    ETASUM .= zero(0.0)
    GGGSUM .= zero(0.0)
    SXYSUM .= zero(0.0)
    COHSUM .= zero(0.0)
    TENSUM .= zero(0.0)
    FRISUM .= zero(0.0)
    WTSUM .= zero(0.0)
    # Vx nodes
    RHOXSUM .= zero(0.0)
    RHOFXSUM .= zero(0.0)
    KXSUM .= zero(0.0)
    PHIXSUM .= zero(0.0)
    RXSUM .= zero(0.0)
    WTXSUM .= zero(0.0)
    # Vy nodes
    RHOYSUM .= zero(0.0)
    RHOFYSUM .= zero(0.0)
    KYSUM .= zero(0.0)
    PHIYSUM .= zero(0.0)
    RYSUM .= zero(0.0)
    WTYSUM .= zero(0.0)
    # P Nodes
    RHOSUM .= zero(0.0)
    RHOCPSUM .= zero(0.0)
    ALPHASUM .= zero(0.0)
    ALPHAFSUM .= zero(0.0)
    HRSUM .= zero(0.0)
    GGGPSUM .= zero(0.0)
    SXXSUM .= zero(0.0)
    TKSUM .= zero(0.0)
    PHISUM .= zero(0.0)
    WTPSUM .= zero(0.0)
    return nothing
end

"""
Reset thermochemical properties to be interpolated to staggered grid.

$(SIGNATURES)

# Details   

    - DMPSUM: interpolation of DMP at P nodes
    - DHPSUM: interpolation of DHP at P nodes
    - WTPSUM: interpolation weights at P nodes
    
# Returns

    - nothing
"""
function reset_thermochemical_properties!(DMPSUM, DHPSUM, WTPSUM, DQPFSUM=nothing)
    # P Nodes
    DMPSUM .= zero(0.0)
    DHPSUM .= zero(0.0)
    WTPSUM .= zero(0.0)
    if DQPFSUM !== nothing
        DQPFSUM .= zero(0.0)
    end
    return nothing
end

"""
Compute top and left grid nodes indices and bilinear interpolation weigths to
nearest four grid nodes for given (x, y) position and grid axes.

$(SIGNATURES)

# Details

    - x: x-position [m]
    - y: y-position [m]
    - x_axis: x-grid reference axis array [m]
    - y_axis: y-grid reference axis array [m]
    - dx: x-grid axis mesh width [m]
    - dy: y-grid axis mesh width [m]
    - jmin: minimum assignable index on x-grid axis (basic/Vx/Vy/P)
    - jmax: maximum assignable index on x-grid axis (basic/Vx/Vy/P)
    - imin: minimum assignable index on y-grid axis (basic/Vx/Vy/P)
    - imax: maximum assignable index on y-grid axis (basic/Vx/Vy/P)

# Returns
    - i: top (with reference to y) node index on y-grid axis
    - j: left (with reference to x) node index on x-grid axis
    - bilinear_weights: vector of 4 bilinear interpolation weights to
      nearest four grid nodes:
        [wtmij  : i  , j   node,
        wtmi1j : i+1, j   node,
        wtmij1 : i  , j+1 node,
        wtmi1j1: i+1, j+1 node]
"""
function fix_weights(x, y, x_axis, y_axis, dx, dy, jmin, jmax, imin, imax)
    i, j, dxmj, dymi = fix_distances(x, y, x_axis, y_axis, dx, dy, jmin, jmax, imin, imax)
    return i,
    j,
    SVector(
        (1.0-dymi/dy) * (1.0-dxmj/dx),
        (dymi/dy) * (1.0-dxmj/dx),
        (1.0-dymi/dy) * (dxmj/dx),
        (dymi/dy) * (dxmj/dx),
    )
end # function fix_weights

"""
Compute top and left grid nodes indices (i, j) and x- and y-distances to that 
grid node (i, j) for given (x, y) position and grid axes.

$(SIGNATURES)

# Details

    - x: x-position [m]
    - y: y-position [m]
    - x_axis: x-grid reference axis array [m]
    - y_axis: y-grid reference axis array [m]
    - dx: x-grid axis mesh width [m]
    - dy: y-grid axis mesh width [m]
    - jmin: minimum assignable index on x-grid axis (basic/Vx/Vy/P)
    - jmax: maximum assignable index on x-grid axis (basic/Vx/Vy/P)
    - imin: minimum assignable index on y-grid axis (basic/Vx/Vy/P)
    - imax: maximum assignable index on y-grid axis (basic/Vx/Vy/P)

# Returns
    - i: top (with reference to y) node index on y-grid axis
    - j: left (with reference to x) node index on x-grid axis
    - dxmj: x-distance from (x, y) point to (i, j) node
    - dymi: y-distance from (x, y) point to (i, j) node
"""
function fix_distances(x, y, x_axis, y_axis, dx, dy, jmin, jmax, imin, imax)
    @inbounds begin
        i, j = fix(x, y, x_axis, y_axis, dx, dy, jmin, jmax, imin, imax)
        dxmj = x - x_axis[j]
        dymi = y - y_axis[i]
    end # @inbounds
    return i, j, dxmj, dymi
end # function fix_distances

"""
Compute top and left grid nodes indices (i, j) and x- and y-distances to that 
grid node (i, j) for given (x, y) position and grid axes.

$(SIGNATURES)

# Details

    - x: x-position [m]
    - y: y-position [m]
    - x_axis: x-grid reference axis array [m]
    - y_axis: y-grid reference axis array [m]
    - dx: x-grid axis mesh width [m]
    - dy: y-grid axis mesh width [m]
    - jmin: minimum assignable index on x-grid axis (basic/Vx/Vy/P)
    - jmax: maximum assignable index on x-grid axis (basic/Vx/Vy/P)
    - imin: minimum assignable index on y-grid axis (basic/Vx/Vy/P)
    - imax: maximum assignable index on y-grid axis (basic/Vx/Vy/P)

# Returns
    - i: top (with reference to y) node index on y-grid axis
    - j: left (with reference to x) node index on x-grid axis
"""
function fix(x, y, x_axis, y_axis, dx, dy, jmin, jmax, imin, imax)
    @inbounds j = unsafe_trunc(Int, (x-x_axis[1])*inv(dx)) + 1
    @inbounds i = unsafe_trunc(Int, (y-y_axis[1])*inv(dy)) + 1
    if j < jmin
        j = jmin
    elseif j > jmax
        j = jmax
    end
    if i < imin
        i = imin
    elseif i > imax
        i = imax
    end
    return i, j
end # function fix

"""
Reduce a 3D (i, j, k) array along its third (k) axis by addition and write the result
into (i, j, 1) without reallocating the array's memory.

$(SIGNATURES)

# Details

    - A: 3D array [i, j, k]

# Returns

    - nothing
"""
function reduce_add_3darray!(A)
    ii = axes(A, 1)
    ij = axes(A, 2)
    ik = axes(A, 3)
    for k in ik[(begin + 1):end], j in ij, i in ii
        @inbounds A[i, j, 1] += A[i, j, k]
    end
    return nothing
end

"""
Interpolate a property to neareast four nodes on a given grid location
using given bilinear interpolation weights.

# Details

    - i: top (with reference to y) node index on vertical y-grid axis
    - j: left (with reference to x) node index on horizontal x-grid axis
    - weights: vector of 4 bilinear interpolation weights to
      nearest four grid nodes:
        [wtmij  : i  , j   node,
        wtmi1j : i+1, j   node,
        wtmij1 : i  , j+1 node,
        wtmi1j1: i+1, j+1 node]
    - property: property to be interpolated to grid using weights
    - grid: threaded grid array on which to interpolate property

# Returns

    - nothing
"""
function interpolate_add_to_grid!(i, j, weights, property, grid)
    # @timeit to "interpolate_add_to_grid!" begin
    @inbounds grid[i, j] += property * weights[1]
    @inbounds grid[i + 1, j] += property * weights[2]
    @inbounds grid[i, j + 1] += property * weights[3]
    @inbounds grid[i + 1, j + 1] += property * weights[4]
    # end # @timeit to "interpolate_add_to_grid!"
    return nothing
end # function interpolate_add_to_grid!

"""
Interpolate a property from nearest the four nodes on a given grid to a marker.

# Details

    - m: number of marker to interpolate to
    - i: top (with reference to y) node index on vertical y-grid axis
    - j: left (with reference to x) node index on horizontal x-grid axis
    - weights: vector of 4 bilinear interpolation weights to
    nearest four grid nodes:
        [wtmij  : i  , j   node,
        wtmi1j : i+1, j   node,
        wtmij1 : i  , j+1 node,
        wtmi1j1: i+1, j+1 node]
    - marker_property: marker property array into which to interpolate
    - property_grid: grid whose property is to be interpolated to marker
    - m: marker 
"""
function interpolate_to_marker!(m, i, j, weights, marker_property, grid)
    @inbounds marker_property[m] = dot4(grid_vector(i, j, grid), weights)
    return nothing
end # function interpolate_to_marker

"""
Interpolate a property from nearest the four nodes on a given grid to a marker
and add it to the markers property.

# Details

    - m: number of marker to interpolate to
    - i: top (with reference to y) node index on vertical y-grid axis
    - j: left (with reference to x) node index on horizontal x-grid axis
    - weights: vector of 4 bilinear interpolation weights to
    nearest four grid nodes:
        [wtmij  : i  , j   node,
        wtmi1j : i+1, j   node,
        wtmij1 : i  , j+1 node,
        wtmi1j1: i+1, j+1 node]
    - marker_property: marker property array into which to interpolate and add
    - property_grid: grid whose property is to be interpolated to marker
    - m: marker 
"""
function interpolate_add_to_marker!(m, i, j, weights, marker_property, grid)
    @inbounds marker_property[m] += dot4(grid_vector(i, j, grid), weights)
    return nothing
end # function interpolate_add_to_marker

"""
Interpolate selected marker properties to basic nodes.

$(SIGNATURES)

# Details

    - m: marker number
    - xmm: marker's x-position [m]
    - ymm: marker's y-position [m]
    - etatotalm: total viscosity of markers
    - etavpm: matrix viscosity of markers
    - inv_gggtotalm: inverse of total shear modulus of markers
    - sxym: marker σxy [Pa]
    - cohestotalm: total compressive strength of markers
    - tenstotalm: total tensile strength of markers 
    - fricttotalm: total friction coefficient of markers
    - ETA0SUM: viscous viscosity interpolated to basic nodes
    - ETASUM: viscoplastic viscosity interpolated to basic nodes
    - GGGSUM: shear modulus interpolated to basic nodes
    - SXYSUM: σxy shear stress interpolated to basic nodes 
    - COHSUM: compressive strength interpolated to basic nodes
    - TENSUM: tensile strength interpolated to basic nodes
    - FRISUM: friction interpolated to basic nodes 
    - WTSUM: weight array for bilinear interpolation to basic nodes

# Returns

    - nothing
"""
function marker_to_basic_nodes!(
    m,
    xmm,
    ymm,
    etatotalm,
    etavpm,
    inv_gggtotalm,
    sxym,
    cohestotalm,
    tenstotalm,
    fricttotalm,
    ETA0SUM,
    ETASUM,
    GGGSUM,
    SXYSUM,
    COHSUM,
    TENSUM,
    FRISUM,
    WTSUM;
    coords=nothing,
)
    x_val = coords === nothing ? x : coords.x
    y_val = coords === nothing ? y : coords.y
    dx_val = coords === nothing ? dx : coords.dx
    dy_val = coords === nothing ? dy : coords.dy
    jmin_basic_val = coords === nothing ? jmin_basic : coords.jmin_basic
    jmax_basic_val = coords === nothing ? jmax_basic : coords.jmax_basic
    imin_basic_val = coords === nothing ? imin_basic : coords.imin_basic
    imax_basic_val = coords === nothing ? imax_basic : coords.imax_basic
    i, j, weights = fix_weights(
        xmm,
        ymm,
        x_val,
        y_val,
        dx_val,
        dy_val,
        jmin_basic_val,
        jmax_basic_val,
        imin_basic_val,
        imax_basic_val,
    )
    @inbounds begin
        interpolate_add_to_grid!(i, j, weights, etatotalm[m], ETA0SUM)
        interpolate_add_to_grid!(i, j, weights, etavpm[m], ETASUM)
        interpolate_add_to_grid!(i, j, weights, inv_gggtotalm[m], GGGSUM)
        interpolate_add_to_grid!(i, j, weights, sxym[m], SXYSUM)
        interpolate_add_to_grid!(i, j, weights, cohestotalm[m], COHSUM)
        interpolate_add_to_grid!(i, j, weights, tenstotalm[m], TENSUM)
        interpolate_add_to_grid!(i, j, weights, fricttotalm[m], FRISUM)
        interpolate_add_to_grid!(i, j, weights, one(1.0), WTSUM)
    end # @inbounds
    return nothing
end

"""
Interpolate selected marker properties to Vx nodes.

$(SIGNATURES)

# Details

    - m: marker number
    - xmm: marker's x-position [m]
    - ymm: marker's y-position [m]
    - rhototalm: total density of markers
    - rhofluidcur: fluid density of markers
    - ktotalm: total thermal conductivity of markers
    - phim: porosity of markers
    - etafluidcur_inv_kphim: fluid viscosity over permeability of markers
    - RHOXSUM: density interpolated to Vx nodes
    - RHOFXSUM: fluid density interpolated to Vx nodes
    - KXSUM: thermal conductivity interpolated to Vx nodes
    - PHIXSUM: porosity interpolated to Vx nodes
    - RXSUM: ηfluid/kϕ interpolated to Vx nodes
    - WTXSUM: weight for bilinear interpolation to Vx nodes

# Returns

    - nothing
"""
function marker_to_vx_nodes!(
    m,
    xmm,
    ymm,
    rhototalm,
    rhofluidcur,
    ktotalm,
    phim,
    etafluidcur_inv_kphim,
    RHOXSUM,
    RHOFXSUM,
    KXSUM,
    PHIXSUM,
    RXSUM,
    WTXSUM;
    coords=nothing,
)
    xvx_val = coords === nothing ? xvx : coords.xvx
    yvx_val = coords === nothing ? yvx : coords.yvx
    dx_val = coords === nothing ? dx : coords.dx
    dy_val = coords === nothing ? dy : coords.dy
    jmin_vx_val = coords === nothing ? jmin_vx : coords.jmin_vx
    jmax_vx_val = coords === nothing ? jmax_vx : coords.jmax_vx
    imin_vx_val = coords === nothing ? imin_vx : coords.imin_vx
    imax_vx_val = coords === nothing ? imax_vx : coords.imax_vx
    i, j, weights = fix_weights(
        xmm,
        ymm,
        xvx_val,
        yvx_val,
        dx_val,
        dy_val,
        jmin_vx_val,
        jmax_vx_val,
        imin_vx_val,
        imax_vx_val,
    )
    @inbounds begin
        interpolate_add_to_grid!(i, j, weights, rhototalm[m], RHOXSUM)
        interpolate_add_to_grid!(i, j, weights, rhofluidcur[m], RHOFXSUM)
        interpolate_add_to_grid!(i, j, weights, ktotalm[m], KXSUM)
        interpolate_add_to_grid!(i, j, weights, phim[m], PHIXSUM)
        interpolate_add_to_grid!(i, j, weights, etafluidcur_inv_kphim[m], RXSUM)
        interpolate_add_to_grid!(i, j, weights, one(1.0), WTXSUM)
    end # @inbounds
    return nothing
end

"""
Interpolate selected marker properties to Vy nodes.

$(SIGNATURES)

# Details

    - m: marker number
    - xmm: marker's x-position [m]
    - ymm: marker's y-position [m]
    - rhototalm: total density of markers
    - rhofluidcur: fluid density of markers
    - ktotalm: total thermal conductivity of markers
    - phim: porosity of markers
    - etafluidcur_inv_kphim: fluid viscosity over permeability of markers
    - RHOYSUM: density interpolated to Vy nodes
    - RHOFYSUM: fluid density interpolated to Vy nodes
    - KYSUM: thermal conductivity interpolated to Vy nodes
    - PHIYSUM: porosity interpolated to Vy nodes
    - RYSUM: ηfluid/kϕ interpolated to Vy nodes
    - WTYSUM: weight for bilinear interpolation to Vy nodes

# Returns

    - nothing
"""
function marker_to_vy_nodes!(
    m,
    xmm,
    ymm,
    rhototalm,
    rhofluidcur,
    ktotalm,
    phim,
    etafluidcur_inv_kphim,
    RHOYSUM,
    RHOFYSUM,
    KYSUM,
    PHIYSUM,
    RYSUM,
    WTYSUM;
    coords=nothing,
)
    xvy_val = coords === nothing ? xvy : coords.xvy
    yvy_val = coords === nothing ? yvy : coords.yvy
    dx_val = coords === nothing ? dx : coords.dx
    dy_val = coords === nothing ? dy : coords.dy
    jmin_vy_val = coords === nothing ? jmin_vy : coords.jmin_vy
    jmax_vy_val = coords === nothing ? jmax_vy : coords.jmax_vy
    imin_vy_val = coords === nothing ? imin_vy : coords.imin_vy
    imax_vy_val = coords === nothing ? imax_vy : coords.imax_vy
    i, j, weights = fix_weights(
        xmm,
        ymm,
        xvy_val,
        yvy_val,
        dx_val,
        dy_val,
        jmin_vy_val,
        jmax_vy_val,
        imin_vy_val,
        imax_vy_val,
    )
    @inbounds begin
        interpolate_add_to_grid!(i, j, weights, rhototalm[m], RHOYSUM)
        interpolate_add_to_grid!(i, j, weights, rhofluidcur[m], RHOFYSUM)
        interpolate_add_to_grid!(i, j, weights, ktotalm[m], KYSUM)
        interpolate_add_to_grid!(i, j, weights, phim[m], PHIYSUM)
        interpolate_add_to_grid!(i, j, weights, etafluidcur_inv_kphim[m], RYSUM)
        interpolate_add_to_grid!(i, j, weights, one(1.0), WTYSUM)
    end # @inbounds
    return nothing
end

"""
Interpolate selected marker properties to P nodes.

$(SIGNATURES)

# Details

    - m: marker number
    - xmm: marker x-position [m]
    - ymm: marker y-position [m]
	- inv_gggtotalm: inverse of total shear modulus of markers
	- sxxm: marker σ′xx [Pa]
	- rhototalm: total density of markers
	- rhocptotalm: total volumetric heat capacity of markers 
	- alphasolidcur: solid thermal expansion coefficient of markers 
	- alphafluidcur: fluid thermal expansion coefficient of markers
	- hrtotalm: total radiogenic heat production of markers
	- phim:  marker porosity
	- tkm_rhocptotalm: total thermal energy of markers
    - GGGPSUM: shear modulus interpolated to P nodes
    - SXXSUM: σ'xx interpolated to P nodes
    - RHOSUM: density interpolated to P nodes
    - RHOCPSUM: volumetric heat capacity interpolated to P nodes
    - ALPHASUM: thermal expansion interpolated to P nodes
    - ALPHAFSUM: fluid thermal expansion interpolated to P nodes
    - HRSUM: radioactive heating interpolated to P nodes
    - PHISUM: porosity interpolated to P nodes
    - TKSUM: heat capacity interpolated to P nodes
    - WTPSUM: weight for bilinear interpolation to P nodes

# Returns

    - nothing
"""
function marker_to_p_nodes!(
    m,
    xmm,
    ymm,
    inv_gggtotalm,
    sxxm,
    rhototalm,
    rhocptotalm,
    alphasolidcur,
    alphafluidcur,
    hrtotalm,
    phim,
    tkm_rhocptotalm,
    GGGPSUM,
    SXXSUM,
    RHOSUM,
    RHOCPSUM,
    ALPHASUM,
    ALPHAFSUM,
    HRSUM,
    PHISUM,
    TKSUM,
    WTPSUM;
    coords=nothing,
)
    xp_val = coords === nothing ? xp : coords.xp
    yp_val = coords === nothing ? yp : coords.yp
    dx_val = coords === nothing ? dx : coords.dx
    dy_val = coords === nothing ? dy : coords.dy
    jmin_p_val = coords === nothing ? jmin_p : coords.jmin_p
    jmax_p_val = coords === nothing ? jmax_p : coords.jmax_p
    imin_p_val = coords === nothing ? imin_p : coords.imin_p
    imax_p_val = coords === nothing ? imax_p : coords.imax_p
    i, j, weights = fix_weights(
        xmm,
        ymm,
        xp_val,
        yp_val,
        dx_val,
        dy_val,
        jmin_p_val,
        jmax_p_val,
        imin_p_val,
        imax_p_val,
    )
    @inbounds begin
        interpolate_add_to_grid!(i, j, weights, inv_gggtotalm[m], GGGPSUM)
        interpolate_add_to_grid!(i, j, weights, sxxm[m], SXXSUM)
        interpolate_add_to_grid!(i, j, weights, rhototalm[m], RHOSUM)
        interpolate_add_to_grid!(i, j, weights, rhocptotalm[m], RHOCPSUM)
        interpolate_add_to_grid!(i, j, weights, alphasolidcur[m], ALPHASUM)
        interpolate_add_to_grid!(i, j, weights, alphafluidcur[m], ALPHAFSUM)
        interpolate_add_to_grid!(i, j, weights, hrtotalm[m], HRSUM)
        interpolate_add_to_grid!(i, j, weights, phim[m], PHISUM)
        interpolate_add_to_grid!(i, j, weights, tkm_rhocptotalm[m], TKSUM)
        interpolate_add_to_grid!(i, j, weights, one(1.0), WTPSUM)
    end # @inbounds
    return nothing
end

"""
Interpolate marker wet silicate (solid) fraction to P nodes.

$(SIGNATURES)

# Details

    - m: marker number
    - xmm: marker x-position [m]
    - ymm: marker y-position [m]
    - XWsolidm0: previous marker wet silicate (solid) molar fraction
    - XWSSUM: wet silicate (solid) molar fraction interpolated to P nodes
    - WTPSUM: weight for bilinear interpolation to P nodes

# Returns

    - nothing
"""
function molarfraction_marker_to_p_nodes!(
    m, xmm, ymm, XWsolidm0, XWSSUM, WTPSUM; coords=nothing
)
    xp_val = coords === nothing ? xp : coords.xp
    yp_val = coords === nothing ? yp : coords.yp
    dx_val = coords === nothing ? dx : coords.dx
    dy_val = coords === nothing ? dy : coords.dy
    jmin_p_val = coords === nothing ? jmin_p : coords.jmin_p
    jmax_p_val = coords === nothing ? jmax_p : coords.jmax_p
    imin_p_val = coords === nothing ? imin_p : coords.imin_p
    imax_p_val = coords === nothing ? imax_p : coords.imax_p
    i, j, weights = fix_weights(
        xmm,
        ymm,
        xp_val,
        yp_val,
        dx_val,
        dy_val,
        jmin_p_val,
        jmax_p_val,
        imin_p_val,
        imax_p_val,
    )
    @inbounds begin
        interpolate_add_to_grid!(i, j, weights, XWsolidm0[m], XWSSUM)
        interpolate_add_to_grid!(i, j, weights, one(1.0), WTPSUM)
    end # @inbounds
    return nothing
end

"""
Compute P nodes wet silicate (solid) fraction from markers.

$(SIGNATURES)

# Details

    - xm: x-coordinate of markers
    - ym: y-coordinate of markers
    - XWsolidm0: previous marker wet silicate (solid) molar fraction
    - XWS: wet silicate (solid) molar fraction at P nodes
    - XWSSUM: wet silicate (solid) molar fraction interpolated to P nodes
    - WTPSUM: weight for bilinear interpolation to P nodes
    - marknum: current total number of markers in use

# Returns

    - nothing
"""
function update_p_nodes_melt_composition!(
    xm, ym, XWsolidm0, XWS, XWSSUM, WTPSUM, marknum; coords=nothing
)
    XWSSUM .= zero(0.0)
    WTPSUM .= zero(0.0)
    for m in 1:1:marknum
        molarfraction_marker_to_p_nodes!(
            m, xm[m], ym[m], XWsolidm0, XWSSUM, WTPSUM; coords=coords
        )
    end
    compute_molarfraction!(XWSSUM, WTPSUM, XWS)
    return nothing
end

"""
Compute properties of basic nodes based on interpolation arrays.

$(SIGNATURES)

# Details

    - ETA0SUM: ETA0 interpolation array
    - ETASUM: ETA interpolation array
    - GGGSUM: GGG interpolation array
    - SXYSUM: SXY interpolation array
    - COHSUM: COH interpolation array
    - TENSUM: TEN interpolation array
    - FRISUM: FRI interpolation array
    - WTSUM: WT interpolation array
    - ETA0: ETA0 basic node array
    - ETA: ETA basic node array
    - GGG: GGG basic node array
    - SXY0: SXY basic node array
    - COH: COH basic node array
    - TEN: TEN basic node array
    - FRI: FRI basic node array
    - YNY: YNY basic node array

# Returns

    - nothing

"""
function compute_basic_node_properties!(
    ETA0SUM,
    ETASUM,
    GGGSUM,
    SXYSUM,
    COHSUM,
    TENSUM,
    FRISUM,
    WTSUM,
    ETA0,
    ETA,
    GGG,
    SXY0,
    COH,
    TEN,
    FRI,
    YNY,
)
    # @timeit to "compute_basic_node_properties!" begin
    Ny, Nx = size(WTSUM)
    @inbounds begin
        for j in 1:1:Nx, i in 1:1:Ny
            if WTSUM[i, j] > 0.0
                ETA0[i, j] = ETA0SUM[i, j] * inv(WTSUM[i, j])
                ETA[i, j] = ETASUM[i, j] * inv(WTSUM[i, j])
                if ETA[i, j] < ETA0[i, j]
                    YNY[i, j] = true
                end
                GGG[i, j] = inv(GGGSUM[i, j]) * WTSUM[i, j]
                SXY0[i, j] = SXYSUM[i, j] * inv(WTSUM[i, j])
                COH[i, j] = COHSUM[i, j] * inv(WTSUM[i, j])
                TEN[i, j] = TENSUM[i, j] * inv(WTSUM[i, j])
                FRI[i, j] = FRISUM[i, j] * inv(WTSUM[i, j])
            end
        end
    end # @inbounds
    # end # @timeit to "compute_basic_node_properties!"
    return nothing
end # function compute_basic_node_properties!

"""
Compute properties of Vx nodes based on interpolation arrays.

$(SIGNATURES)

# Details

    - RHOXSUM: RHOX interpolation array
    - RHOFXSUM: RHOFX interpolation array
    - KXSUM: KX interpolation array
    - PHIXSUM: PHIX interpolation array
    - RXSUM: RX interpolation array
    - WTXSUM: WTX interpolation array
    - RHOX: RHOX Vx node array
    - RHOFX: RHOFX Vx node array
    - KX: KX Vx node array
    - PHIX: PHIX Vx node array
    - RX: RX Vx node array

# Returns

    - nothing

"""
function compute_vx_node_properties!(
    RHOXSUM, RHOFXSUM, KXSUM, PHIXSUM, RXSUM, WTXSUM, RHOX, RHOFX, KX, PHIX, RX
)
    # @timeit to "compute_vx_node_properties!" begin
    Ny1, Nx1 = size(WTXSUM)
    @inbounds begin
        for j in 1:1:Nx1, i in 1:1:Ny1
            if WTXSUM[i, j] > 0.0
                RHOX[i, j] = RHOXSUM[i, j] * inv(WTXSUM[i, j])
                RHOFX[i, j] = RHOFXSUM[i, j] * inv(WTXSUM[i, j])
                KX[i, j] = KXSUM[i, j] * inv(WTXSUM[i, j])
                PHIX[i, j] = PHIXSUM[i, j] * inv(WTXSUM[i, j])
                RX[i, j] = RXSUM[i, j] * inv(WTXSUM[i, j])
            end
        end
    end # @inbounds
    # end # @timeit to "compute_vx_node_properties!"
    return nothing
end # function compute_vx_node_properties!

"""
Compute properties of Vy nodes based on interpolation arrays.

$(SIGNATURES)

# Details

    - RHOYSUM: RHOY interpolation array
    - RHOFYSUM: RHOFY interpolation array
    - KYSUM: KY interpolation array
    - PHIYSUM: PHIY interpolation array
    - RYSUM: RY interpolation array
    - WTYSUM: WTY interpolation array
    - RHOY: RHOY Vy node array
    - RHOFY: RHOFY Vy node array
    - KY: KY Vy node array
    - PHIY: PHIY Vy node array
    - RY: RY Vy node array

# Returns

    - nothing

"""
function compute_vy_node_properties!(
    RHOYSUM, RHOFYSUM, KYSUM, PHIYSUM, RYSUM, WTYSUM, RHOY, RHOFY, KY, PHIY, RY
)
    # @timeit to "compute_vy_node_properties!" begin
    Ny1, Nx1 = size(WTYSUM)
    @inbounds begin
        for j in 1:1:Nx1, i in 1:1:Ny1
            if WTYSUM[i, j] > 0.0
                RHOY[i, j] = RHOYSUM[i, j] * inv(WTYSUM[i, j])
                RHOFY[i, j] = RHOFYSUM[i, j] * inv(WTYSUM[i, j])
                KY[i, j] = KYSUM[i, j] * inv(WTYSUM[i, j])
                PHIY[i, j] = PHIYSUM[i, j] * inv(WTYSUM[i, j])
                RY[i, j] = RYSUM[i, j] * inv(WTYSUM[i, j])
            end
        end
    end # @inbounds
    # end # @timeit to "compute_vy_node_properties!"
    return nothing
end # function compute_vy_node_properties!

"""
Compute properties of P nodes based on interpolation arrays.

$(SIGNATURES)

# Details

    - GGGPSUM: GGGP interpolation array
    - SXX0SUM: SXX0 interpolation array
    - RHOSUM: RHO interpolation array
    - RHOCPSUM: RHOCP interpolation array
    - ALPHASUM: ALPHA interpolation array
    - ALPHAFSUM: ALPHAF interpolation array
    - HRSUM: HR interpolation array
    - PHISUM: PHI interpolation array
    - TKSUM: TK interpolation array
    - WTPSUM: WTP interpolation array
    - GGGP: GGGP P node array
    - SXX0: SXX0 P node array
    - RHO: RHO P node array
    - RHOCP: RHOCP P node array
    - ALPHA: ALPHA P node array
    - ALPHAF: ALPHAF P node array
    - HR: HR P node array
    - PHI: PHI P node array
    - BETAPHI: BETAPHI P node array
    - tk1: tk1 P node array

# Returns

    - nothing

"""
function compute_p_node_properties!(
    RHOSUM,
    RHOCPSUM,
    ALPHASUM,
    ALPHAFSUM,
    HRSUM,
    GGGPSUM,
    SXXSUM,
    TKSUM,
    PHISUM,
    WTPSUM,
    RHO,
    RHOCP,
    ALPHA,
    ALPHAF,
    HR,
    GGGP,
    SXX0,
    tk1,
    PHI,
    BETAPHI,
)
    # @timeit to "compute_p_node_properties!" begin
    Ny1, Nx1 = size(WTPSUM)
    @inbounds begin
        for j in 1:1:Nx1, i in 1:1:Ny1
            if WTPSUM[i, j] > 0.0
                RHO[i, j] = RHOSUM[i, j] * inv(WTPSUM[i, j])
                RHOCP[i, j] = RHOCPSUM[i, j] * inv(WTPSUM[i, j])
                ALPHA[i, j] = ALPHASUM[i, j] * inv(WTPSUM[i, j])
                ALPHAF[i, j] = ALPHAFSUM[i, j] * inv(WTPSUM[i, j])
                HR[i, j] = HRSUM[i, j] * inv(WTPSUM[i, j])
                GGGP[i, j] = inv(GGGPSUM[i, j]) * WTPSUM[i, j]
                SXX0[i, j] = SXXSUM[i, j] * inv(WTPSUM[i, j])
                tk1[i, j] = TKSUM[i, j] * inv(RHOCPSUM[i, j])
                PHI[i, j] = PHISUM[i, j] * inv(WTPSUM[i, j])
                BETAPHI[i, j] = inv(GGGP[i, j]) * PHI[i, j]
            end
        end
    end # @inbounds
    # end # @timeit to "compute_p_node_properties!"
    return nothing
end # function compute_p_node_properties!

"""
Compute wet solid molar fraction at P nodes based on interpolation arrays.

$(SIGNATURES)

# Details

    - XWSSUM: XWX interpolation array
    - WTPSUM: WTP interpolation array
    - XWS: wet solid molar fraction at P nodes

# Returns

    - nothing
"""
function compute_molarfraction!(XWSSUM, WTPSUM, XWS)
    # @timeit to "compute_molarfraction!" begin
    Ny1, Nx1 = size(WTPSUM)
    @inbounds begin
        for j in 1:1:Nx1, i in 1:1:Ny1
            if WTPSUM[i, j] > 0.0
                XWS[i, j] = XWSSUM[i, j] * inv(WTPSUM[i, j])
            else
                XWS[i, j] = zero(0.0)
            end
        end
    end # @inbounds
    # end # @timeit to "compute_molarfraction!"
    return nothing
end # function compute_molarfraction!

"""
Add a RK4 stage velocity to the RK4 velocity vector.

$(SIGNATURES)

# Details:

    - vrk4: current RK4 velocity vector
    - v: RK4 velocity of stage `rk` to be added to velocity vector
    - rk: RK4 stage number

# Returns

    - vrk4: updated RK4 velocity vector
"""
function add_vrk4(vrk4, v, rk)
    if rk == 1
        return vrk4 + @SVector [v, 0.0, 0.0, 0.0]
    elseif rk == 2
        return vrk4 + @SVector [0.0, v, 0.0, 0.0]
    elseif rk == 3
        return vrk4 + @SVector [0.0, 0.0, v, 0.0]
    elseif rk == 4
        return vrk4 + @SVector [0.0, 0.0, 0.0, v]
    else
        return vrk4
    end
end

""" 
Compute solid velocities, fluid velocities at P nodes.

$(SIGNATURES)

# Details

    - vx: solid vx-velocity at Vx nodes
    - vy: solid vy-velocity at Vy nodes
    - vxf: fluid vx-velocity at Vx nodes
    - vyf: fluid vy-velocity at Vy nodes
    - vxp: solid vx-velocity at P nodes
    - vyp: solid vy-velocity at P nodes
    - vxpf: fluid vx-velocity at P nodes
    - vypf: fluid vy-velocity at P nodes

# Returns

    - nothing
"""
function compute_velocities!(vx, vy, vxf, vyf, vxp, vyp, vxpf, vypf; coords=nothing)
    # @timeit to "compute_velocities!" begin
    Ny1, Nx1 = size(vxp)
    Nx = Nx1 - 1
    Ny = Ny1 - 1
    @inbounds begin
        # compute solid velocities at P nodes
        for j in 2:1:Nx, i in 2:1:Ny
            vxp[i, j] = 0.5 * (vx[i, j] + vx[i, j - 1])
            vyp[i, j] = 0.5 * (vy[i, j] + vy[i - 1, j])
            vxpf[i, j] = 0.5 * (vxf[i, j] + vxf[i, j - 1])
            vypf[i, j] = 0.5 * (vyf[i, j] + vyf[i - 1, j])
        end
        # apply boundary conditions
        # vxp
        # top: free slip
        @views @. vxp[1, 2:(Nx - 1)] = - bctop * vxp[2, 2:(Nx - 1)]
        # bottom: free slip
        @views @. vxp[Ny1, 2:(Nx - 1)] = - bcbottom * vxp[Ny, 2:(Nx - 1)]
        # left
        @views @. vxp[:, 1] = 2.0*vxleft - vxp[:, 2]
        # right
        @views @. vxp[:, Nx1] = 2.0*vxright - vxp[:, Nx]
        # vyp
        # left: free slip
        @views @. vyp[2:(Ny - 1), 1] = - bcleft * vyp[2:(Ny - 1), 2]
        # right: free slip
        @views @. vyp[2:(Ny - 1), Nx1] = - bcright * vyp[2:(Ny - 1), Nx]
        # top
        @views @. vyp[1, :] = 2.0*vytop - vyp[2, :]
        # bottom
        @views @. vyp[Ny1, :] = 2.0*vybottom - vyp[Ny, :]
        # vxpf
        # top: free slip
        @views @. vxpf[1, 2:(Nx - 1)] = - bcftop * vxpf[2, 2:(Nx - 1)]
        # bottom: free slip
        @views @. vxpf[Ny1, 2:(Nx - 1)] = - bcfbottom * vxpf[Ny, 2:(Nx - 1)]
        # left
        @views @. vxpf[:, 1] = 2.0*vxleft - vxpf[:, 2]
        # right
        @views @. vxpf[:, Nx1] = 2.0*vxright - vxpf[:, Nx]
        # vypf
        # left: free slip
        @views @. vypf[2:(Ny - 1), 1] = - bcfleft * vypf[2:(Ny - 1), 2]
        # right: free slip
        @views @. vypf[2:(Ny - 1), Nx1] = - bcfright * vypf[2:(Ny - 1), Nx]
        # top
        @views @. vypf[1, :] = 2.0*vytop - vypf[2, :]
        # bottom
        @views @. vypf[Ny1, :] = 2.0*vybottom - vypf[Ny, :]
    end # @inbounds
    # end # @timeit to "compute_velocities!"
    return nothing
end # function compute_velocities!

"""
Compute rotatation rate in basic nodes based on velocity derivatives
at Vx and Vy nodes.

$(SIGNATURES)

# Details

    - vx: solid vx-velocity at Vx nodes
    - vy: solid vy-velocity at Vy nodes
    - wyx: rotation rate at basic nodes
  
# Returns

    - nothing
"""
function compute_rotation_rate!(vx, vy, wyx; coords=nothing)
    # @timeit to "compute_rotation_rate!" begin
    Ny, Nx = size(wyx)
    dx_val = coords === nothing ? dx : coords.dx
    dy_val = coords === nothing ? dy : coords.dy
    for j in 1:1:Nx, i in 1:1:Ny
        @inbounds wyx[i, j] =
            0.5 *
            ((vy[i, j + 1]-vy[i, j])*inv(dx_val) - (vx[i + 1, j]-vx[i, j])*inv(dy_val))
    end
    # end # @timeit to "compute_rotation_rate!"
    return nothing
end # function compute_rotation_rate!

"""
Move markers using classic Runge-Kutta integration (RK4) taking into account
two-phase flow and solid-fluid temperatures.

$(SIGNATURES)

# Details

    - xm: x-coordinate of markers
    - ym: y-coordinate of markers
    - tm: type of markers 
    - tkm: temperature of markers
    - phim: porosity of markers  
    - sxxm: marker σ′xx
    - sxym: marker σxy
    - vx: solid vx-velocity at Vx nodes
    - vy: solid vy-velocity at Vy nodes
    - vxf: fluid vx-velocity at Vx nodes
    - vyf: fluid vy-velocity at Vy nodes
    - wyx: rotation rate at basic nodes
    - tk2: next temperature at P nodes
    - marknum: number of markers in use
    - dt: computational time step
    - mode: marker property computation mode
        - 1: dynamic, based on (Touloukian, 1970; Hobbs, 1974;
             Travis and Schubert, 2005)
        - 2: constant parameter rhocpfluidm

# Returns

    - nothing
"""
function move_markers_rk4!(
    xm,
    ym,
    tm,
    tkm,
    phim,
    sxym,
    sxxm,
    vx,
    vy,
    vxf,
    vyf,
    wyx,
    tk2,
    marknum,
    dt,
    mode;
    coords=nothing,
)
    xp_val = coords === nothing ? xp : coords.xp
    yp_val = coords === nothing ? yp : coords.yp
    x_val = coords === nothing ? x : coords.x
    y_val = coords === nothing ? y : coords.y
    xvx_val = coords === nothing ? xvx : coords.xvx
    yvx_val = coords === nothing ? yvx : coords.yvx
    xvy_val = coords === nothing ? xvy : coords.xvy
    yvy_val = coords === nothing ? yvy : coords.yvy
    dx_val = coords === nothing ? dx : coords.dx
    dy_val = coords === nothing ? dy : coords.dy
    jmin_p_val = coords === nothing ? jmin_p : coords.jmin_p
    jmax_p_val = coords === nothing ? jmax_p : coords.jmax_p
    imin_p_val = coords === nothing ? imin_p : coords.imin_p
    imax_p_val = coords === nothing ? imax_p : coords.imax_p
    jmin_basic_val = coords === nothing ? jmin_basic : coords.jmin_basic
    jmax_basic_val = coords === nothing ? jmax_basic : coords.jmax_basic
    imin_basic_val = coords === nothing ? imin_basic : coords.imin_basic
    imax_basic_val = coords === nothing ? imax_basic : coords.imax_basic
    jmin_vx_val = coords === nothing ? jmin_vx : coords.jmin_vx
    jmax_vx_val = coords === nothing ? jmax_vx : coords.jmax_vx
    imin_vx_val = coords === nothing ? imin_vx : coords.imin_vx
    imax_vx_val = coords === nothing ? imax_vx : coords.imax_vx
    jmin_vy_val = coords === nothing ? jmin_vy : coords.jmin_vy
    jmax_vy_val = coords === nothing ? jmax_vy : coords.jmax_vy
    imin_vy_val = coords === nothing ? imin_vy : coords.imin_vy
    imax_vy_val = coords === nothing ? imax_vy : coords.imax_vy
    Nx_val = coords === nothing ? Nx : coords.Nx
    Ny_val = coords === nothing ? Ny : coords.Ny

    # @timeit to "move_markers_rk4!" begin
    @inbounds begin
        Threads.@threads :static for m in 1:1:marknum
            xmm = xmrk4 = xm[m]
            ymm = ymrk4 = ym[m]
            i, j, weights = fix_weights(
                xmrk4,
                ymrk4,
                xp_val,
                yp_val,
                dx_val,
                dy_val,
                jmin_p_val,
                jmax_p_val,
                imin_p_val,
                imax_p_val,
            )
            # interpolate local solid temperature from P grid
            tksm₀ = dot4(grid_vector(i, j, tk2), weights)
            i, j, weights = fix_weights(
                xmrk4,
                ymrk4,
                x_val,
                y_val,
                dx_val,
                dy_val,
                jmin_basic_val,
                jmax_basic_val,
                imin_basic_val,
                imax_basic_val,
            )
            # interpolate local rotation rate from basic grid
            ωm = dot4(grid_vector(i, j, wyx), weights)
            # incremental rotation angle
            θ = dt * ωm
            # compute analytic stress rotation using σ′′xx = -σ′′yy
            sxxm₁ = sxxm[m]*cos(θ)^2 - sxxm[m]*sin(θ)^2-sxym[m]*sin(2.0*θ)
            sxym₁ = sxxm[m]*sin(2.0*θ) + sxym[m]*cos(2.0*θ)
            # update stresses
            sxxm[m] = sxxm₁
            sxym[m] = sxym₁
            # setup RK4 scheme 
            # RK4 coordinate positions A, B, C, D
            # RK4 velocities va, vb, vc, vd
            vxrk4 = @SVector zeros(4)
            vyrk4 = @SVector zeros(4)
            # advance marker using RK4 scheme on solid velocity
            for rk in 1:1:4
                # interpolate vx
                i, j, dxmj, dymi = fix_distances(
                    xmrk4,
                    ymrk4,
                    xvx_val,
                    yvx_val,
                    dx_val,
                    dy_val,
                    jmin_vx_val,
                    jmax_vx_val,
                    imin_vx_val,
                    imax_vx_val,
                )
                # compute vx velocity for left and right of current cell
                vxm₁₃ = vx[i, j]*(1.0-dxmj/dx_val) + vx[i, j + 1]*dxmj/dx_val
                vxm₂₄ = vx[i + 1, j]*(1.0-dxmj/dx_val) + vx[i + 1, j + 1]*dxmj/dx_val
                # compute second order vx velocity corrections
                if dxmj/dx_val >= 0.5
                    # in right half of cell but not at right edge of grid
                    if j < Nx_val-1
                        vxm₁₃ +=
                            0.5 *
                            ((dxmj/dx_val-0.5)^2) *
                            (vx[i, j] - 2.0*vx[i, j + 1] + vx[i, j + 2])
                        vxm₂₄ +=
                            0.5 *
                            ((dxmj/dx_val-0.5)^2) *
                            (vx[i + 1, j] - 2.0*vx[i + 1, j + 1] + vx[i + 1, j + 2])
                    end
                else
                    # in left half of cell but not at left edge of grid
                    if j > 1
                        vxm₁₃ +=
                            0.5 *
                            ((dxmj/dx_val-0.5)^2) *
                            (vx[i, j - 1] - 2.0*vx[i, j] + vx[i, j + 1])
                        vxm₂₄ +=
                            0.5 *
                            ((dxmj/dx_val-0.5)^2) *
                            (vx[i + 1, j - 1] - 2.0*vx[i + 1, j] + vx[i + 1, j + 1])
                    end
                end
                # compute current RK step vx
                vxrk4 = add_vrk4(vxrk4, vxm₁₃*(1.0-dymi/dy_val) + vxm₂₄*dymi/dy_val, rk)
                # interpolate vy
                i, j, dxmj, dymi = fix_distances(
                    xmrk4,
                    ymrk4,
                    xvy_val,
                    yvy_val,
                    dx_val,
                    dy_val,
                    jmin_vy_val,
                    jmax_vy_val,
                    imin_vy_val,
                    imax_vy_val,
                )
                # compute vy velocity for top and bottom of current cell
                vym₁₂ = vy[i, j]*(1.0-dymi/dy_val) + vy[i + 1, j]*dymi/dy_val
                vym₃₄ = vy[i, j + 1]*(1.0-dymi/dy_val) + vy[i + 1, j + 1]*dymi/dy_val
                # compute second order vy velocity corrections
                if dymi/dy_val >= 0.5
                    # in bottom half of cell but not at bottom edge of grid
                    if i < Ny_val-1
                        vym₁₂ +=
                            0.5 *
                            ((dymi/dy_val-0.5)^2) *
                            (vy[i, j] - 2.0*vy[i + 1, j] + vy[i + 2, j])
                        vym₃₄ +=
                            0.5 *
                            ((dymi/dy_val-0.5)^2) *
                            (vy[i, j + 1] - 2.0*vy[i + 1, j + 1] + vy[i + 2, j + 1])
                    end
                else
                    # in top half of cell but not at top edge of grid
                    if i > 1
                        vym₁₂ +=
                            0.5 *
                            ((dymi/dy_val-0.5)^2) *
                            (vy[i - 1, j] - 2.0*vy[i, j] + vy[i + 1, j])
                        vym₃₄ +=
                            0.5 *
                            ((dymi/dy_val-0.5)^2) *
                            (vy[i - 1, j + 1] - 2.0*vy[i, j + 1] + vy[i + 1, j + 1])
                    end
                end
                # compute current RK step vy
                vyrk4 = add_vrk4(vyrk4, vym₁₂*(1.0-dxmj/dx_val) + vym₃₄*dxmj/dx_val, rk)
                # calculate next RK step x and y positions if not at final step
                if rk < 4
                    xmrk4 = xmm + dt*crk4[rk]*vxrk4[rk]
                    ymrk4 = ymm + dt*crk4[rk]*vyrk4[rk]
                end
            end # RK4 solid velocity loop
            # advance marker using RK4 solid velocity
            xm[m] += dt * dot4(brk4, vxrk4)
            ym[m] += dt * dot4(brk4, vyrk4)
            # reset RK4 scheme for fluid velocity backtracing
            xmm = xmrk4 = xm[m]
            ymm = ymrk4 = ym[m]
            vxrk4 = @SVector zeros(4)
            vyrk4 = @SVector zeros(4)
            # backtrack marker using RK4 scheme on fluid velocity
            for rk in 1:1:4
                # interpolate vxf
                i, j, dxmj, dymi = fix_distances(
                    xmrk4,
                    ymrk4,
                    xvx_val,
                    yvx_val,
                    dx_val,
                    dy_val,
                    jmin_vx_val,
                    jmax_vx_val,
                    imin_vx_val,
                    imax_vx_val,
                )
                # compute vxf velocity for left and right of current cell
                vxfm₁₃ = vxf[i, j]*(1.0-dxmj/dx_val) + vxf[i, j + 1]*dxmj/dx_val
                vxfm₂₄ = vxf[i + 1, j]*(1.0-dxmj/dx_val) + vxf[i + 1, j + 1]*dxmj/dx_val
                # compute second order vxf velocity corrections
                if dxmj/dx_val >= 0.5
                    # in right half of cell but not at right edge of grid
                    if j < Nx_val-1
                        vxfm₁₃ +=
                            0.5 *
                            ((dxmj/dx_val-0.5)^2) *
                            (vxf[i, j] - 2.0*vxf[i, j + 1] + vxf[i, j + 2])
                        vxfm₂₄ +=
                            0.5 *
                            ((dxmj/dx_val-0.5)^2) *
                            (vxf[i + 1, j] - 2.0*vxf[i + 1, j + 1] + vxf[i + 1, j + 2])
                    end
                else
                    # in left half of cell but not at left edge of grid
                    if j > 1
                        vxfm₁₃ +=
                            0.5 *
                            ((dxmj/dx_val-0.5)^2) *
                            (vxf[i, j - 1] - 2.0*vxf[i, j] + vxf[i, j + 1])
                        vxfm₂₄ +=
                            0.5 *
                            ((dxmj/dx_val-0.5)^2) *
                            (vxf[i + 1, j - 1] - 2.0*vxf[i + 1, j] + vxf[i + 1, j + 1])
                    end
                end
                # compute current RK step vxf
                vxrk4 = add_vrk4(vxrk4, vxfm₁₃*(1.0-dymi/dy_val) + vxfm₂₄*dymi/dy_val, rk)
                # interpolate vyf
                i, j, dxmj, dymi = fix_distances(
                    xmrk4,
                    ymrk4,
                    xvy_val,
                    yvy_val,
                    dx_val,
                    dy_val,
                    jmin_vy_val,
                    jmax_vy_val,
                    imin_vy_val,
                    imax_vy_val,
                )
                # compute vyf velocity for top and bottom of current cell
                vyfm₁₂ = vyf[i, j]*(1.0-dymi/dy_val) + vyf[i + 1, j]*dymi/dy_val
                vyfm₃₄ = vyf[i, j + 1]*(1.0-dymi/dy_val) + vyf[i + 1, j + 1]*dymi/dy_val
                # compute second order vyf velocity corrections
                if dymi/dy_val >= 0.5
                    # in bottom half of cell but not at bottom edge of grid
                    if i < Ny_val-1
                        vyfm₁₂ +=
                            0.5 *
                            ((dymi/dy_val-0.5)^2) *
                            (vyf[i, j] - 2.0*vyf[i + 1, j] + vyf[i + 2, j])
                        vyfm₃₄ +=
                            0.5 *
                            ((dymi/dy_val-0.5)^2) *
                            (vyf[i, j + 1] - 2.0*vyf[i + 1, j + 1] + vyf[i + 2, j + 1])
                    end
                else
                    # in top half of cell but not at top edge of grid
                    if i > 1
                        vyfm₁₂ +=
                            0.5 *
                            ((dymi/dy_val-0.5)^2) *
                            (vyf[i - 1, j] - 2.0*vyf[i, j] + vyf[i + 1, j])
                        vyfm₃₄ +=
                            0.5 *
                            ((dymi/dy_val-0.5)^2) *
                            (vyf[i - 1, j + 1] - 2.0*vyf[i, j + 1] + vyf[i + 1, j + 1])
                    end
                end
                # compute current RK step vyf
                vyrk4 = add_vrk4(vyrk4, vyfm₁₂*(1.0-dxmj/dx_val) + vyfm₃₄*dxmj/dx_val, rk)
                # calculate next RK step x and y positions if not at final step
                if rk < 4
                    xmrk4 = xmm - dt*crk4[rk]*vxrk4[rk]
                    ymrk4 = ymm - dt*crk4[rk]*vyrk4[rk]
                end
            end # RK4 fluid velocity loop
            # backtrace marker using RK4 fluid velocity
            xmrk4 = xmm - dt*dot4(brk4, vxrk4)
            ymrk4 = ymm - dt*dot4(brk4, vyrk4)
            # interpolate fluid temperature at backtraced marker position
            i, j, weights = fix_weights(
                xmrk4,
                ymrk4,
                xp_val,
                yp_val,
                dx_val,
                dy_val,
                jmin_p_val,
                jmax_p_val,
                imin_p_val,
                imax_p_val,
            )
            # interpolate backtraced local fluid temperature from P grid
            tkfm₀ = dot4(grid_vector(i, j, tk2), weights)
            # compute marker fluid-solid temperature difference
            δtkfsm = tkfm₀ - tksm₀
            # correct marker temperature

            tkm[m] =
                (
                    (1.0-phim[m])*tkm[m]*rhocpsolidm[tm[m]] +
                    phim[m]*(tkm[m]+δtkfsm)*compute_rhocpfluidm(tkm[m], mode)
                ) / (
                    (1.0-phim[m])*rhocpsolidm[tm[m]] +
                    phim[m]*compute_rhocpfluidm(tkm[m], mode)
                )
        end # marker loop
    end # @inbounds   
    # end # @timeit to "move_markers_rk4!"
    return nothing
end # function move_markers_rk4!

"""
Backtrack pressure nodes using classic Runge-Kutta integration (RK4) 
to update total, solid, and fluid pressure under consideration of 
two-phase flow velocities.

$(SIGNATURES)

# Details

    - pr: total pressure at P nodes
    - pr0: previous time step total pressure at P nodes
    - ps: solid pressure at P nodes
    - ps0: previous time step solid pressure at P nodes
    - pf: fluid pressure at P nodes
    - pf0: previous time step fluid pressure at P nodes
    - vx: solid vx-velocity at Vx nodes
    - vy: solid vy-velocity at Vy nodes 
    - vxf: fluid vx-velocity at Vx nodes 
    - vyf: fluid vy-velocity at Vy nodes
    - dt: computational time step

# Returns

    - nothing
"""
function backtrace_pressures_rk4!(
    pr, pr0, ps, ps0, pf, pf0, vx, vy, vxf, vyf, dt; coords=nothing
)
    # @timeit to "backtrace_pressures_rk4!" begin
    xp_val = coords === nothing ? xp : coords.xp
    yp_val = coords === nothing ? yp : coords.yp
    xvx_val = coords === nothing ? xvx : coords.xvx
    yvx_val = coords === nothing ? yvx : coords.yvx
    xvy_val = coords === nothing ? xvy : coords.xvy
    yvy_val = coords === nothing ? yvy : coords.yvy
    dx_val = coords === nothing ? dx : coords.dx
    dy_val = coords === nothing ? dy : coords.dy
    jmin_p_val = coords === nothing ? jmin_p : coords.jmin_p
    jmax_p_val = coords === nothing ? jmax_p : coords.jmax_p
    imin_p_val = coords === nothing ? imin_p : coords.imin_p
    imax_p_val = coords === nothing ? imax_p : coords.imax_p
    jmin_vx_val = coords === nothing ? jmin_vx : coords.jmin_vx
    jmax_vx_val = coords === nothing ? jmax_vx : coords.jmax_vx
    imin_vx_val = coords === nothing ? imin_vx : coords.imin_vx
    imax_vx_val = coords === nothing ? imax_vx : coords.imax_vx
    jmin_vy_val = coords === nothing ? jmin_vy : coords.jmin_vy
    jmax_vy_val = coords === nothing ? jmax_vy : coords.jmax_vy
    imin_vy_val = coords === nothing ? imin_vy : coords.imin_vy
    imax_vy_val = coords === nothing ? imax_vy : coords.imax_vy
    Nx_val = coords === nothing ? Nx : coords.Nx
    Ny_val = coords === nothing ? Ny : coords.Ny

    @inbounds begin
        # setup RK4 scheme
        vxm = zeros(4)
        vym = zeros(4)
        # advance pressure generation
        pr0 .= pr
        ps0 .= ps
        # backtrace P nodes: total and solid pressure
        for jj in 2:1:Nx_val, ii in 2:1:Ny_val
            xA = xcur = xp_val[jj]
            yA = ycur = yp_val[ii]
            for rk in 1:1:4
                # interpolate vx
                i, j, dxmj, dymi = fix_distances(
                    xcur,
                    ycur,
                    xvx_val,
                    yvx_val,
                    dx_val,
                    dy_val,
                    jmin_vx_val,
                    jmax_vx_val,
                    imin_vx_val,
                    imax_vx_val,
                )
                # compute vx velocity for left and right of current cell
                vxm₁₃ = vx[i, j]*(1.0-dxmj/dx_val) + vx[i, j + 1]*dxmj/dx_val
                vxm₂₄ = vx[i + 1, j]*(1.0-dxmj/dx_val) + vx[i + 1, j + 1]*dxmj/dx_val
                # compute second order vx velocity corrections
                if dxmj/dx_val >= 0.5
                    if j < Nx_val-1
                        vxm₁₃ +=
                            0.5 *
                            ((dxmj/dx_val-0.5)^2) *
                            (vx[i, j] - 2.0*vx[i, j + 1] + vx[i, j + 2])
                        vxm₂₄ +=
                            0.5 *
                            ((dxmj/dx_val-0.5)^2) *
                            (vx[i + 1, j] - 2.0*vx[i + 1, j + 1] + vx[i + 1, j + 2])
                    end
                else
                    if j > 1
                        vxm₁₃ +=
                            0.5 *
                            ((dxmj/dx_val-0.5)^2) *
                            (vx[i, j - 1] - 2.0*vx[i, j] + vx[i, j + 1])
                        vxm₂₄ +=
                            0.5 *
                            ((dxmj/dx_val-0.5)^2) *
                            (vx[i + 1, j - 1] - 2.0*vx[i + 1, j] + vx[i + 1, j + 1])
                    end
                end
                # compute current RK step vx
                vxm[rk] = (1.0-dymi/dy_val)*vxm₁₃ + (dymi/dy_val)*vxm₂₄
                # interpolate vy
                i, j, dxmj, dymi = fix_distances(
                    xcur,
                    ycur,
                    xvy_val,
                    yvy_val,
                    dx_val,
                    dy_val,
                    jmin_vy_val,
                    jmax_vy_val,
                    imin_vy_val,
                    imax_vy_val,
                )
                # compute vy velocity for left and right of current cell
                vym₁₂ = vy[i, j]*(1.0-dymi/dy_val) + vy[i + 1, j]*dymi/dy_val
                vym₃₄ = vy[i, j + 1]*(1.0-dymi/dy_val) + vy[i + 1, j + 1]*dymi/dy_val
                # compute second order vy velocity corrections
                if dymi/dy_val >= 0.5
                    if i < Ny_val-1
                        vym₁₂ +=
                            0.5 *
                            ((dymi/dy_val-0.5)^2) *
                            (vy[i, j] - 2.0*vy[i + 1, j] + vy[i + 2, j])
                        vym₃₄ +=
                            0.5 *
                            ((dymi/dy_val-0.5)^2) *
                            (vy[i, j + 1] - 2.0*vy[i + 1, j + 1] + vy[i + 2, j + 1])
                    end
                else
                    if i > 1
                        vym₁₂=vym₁₂+1/2*((dymi/dy_val-0.5)^2)*(
                            vy[i - 1, j]-2*vy[i, j]+vy[i + 1, j]
                        )
                        vym₃₄=vym₃₄+1/2*((dymi/dy_val-0.5)^2)*(
                            vy[i - 1, j + 1]-2*vy[i, j + 1]+vy[i + 1, j + 1]
                        )
                    end
                end
                # compute current RK step vy
                vym[rk] = (1.0-dxmj/dx_val)*vym₁₂ + (dxmj/dx_val)*vym₃₄
                # calculate next RK step x and y positions if not at final 
                if rk==1 || rk==2
                    xcur=xA - 0.5*dt*vxm[rk]
                    ycur=yA - 0.5*dt*vym[rk]
                elseif rk==3
                    xcur = xA - dt*vxm[rk]
                    ycur = yA - dt*vym[rk]
                end
            end
            # compute effective velocity using RK4
            xcur = xA - dt*1//6*(vxm[1]+2*vxm[2]+2*vxm[3]+vxm[4])
            ycur = yA - dt*1//6*(vym[1]+2*vym[2]+2*vym[3]+vym[4])
            i, j, weights = fix_weights(
                xcur,
                ycur,
                xp_val,
                yp_val,
                dx_val,
                dy_val,
                jmin_p_val,
                jmax_p_val,
                imin_p_val,
                imax_p_val,
            )
            pr0[ii, jj] = dot4(grid_vector(i, j, pr), weights)
            ps0[ii, jj] = dot4(grid_vector(i, j, ps), weights)
        end # jj, ii total and solid pressure loop
        # backtrace P nodes: fluid pressure
        pf0 .= pf
        for jj in 2:1:Nx_val, ii in 2:1:Ny_val
            # Save initial nodal coordinates
            xA = xcur = xp_val[jj]
            yA = ycur = yp_val[ii]
            for rk in 1:1:4
                # interpolate vx
                i, j, dxmj, dymi = fix_distances(
                    xcur,
                    ycur,
                    xvx_val,
                    yvx_val,
                    dx_val,
                    dy_val,
                    jmin_vx_val,
                    jmax_vx_val,
                    imin_vx_val,
                    imax_vx_val,
                )
                # compute vx velocity for left and right of current cell
                vxm₁₃ = vxf[i, j]*(1.0-dxmj/dx_val) + vxf[i, j + 1]*dxmj/dx_val
                vxm₂₄ = vxf[i + 1, j]*(1.0-dxmj/dx_val) + vxf[i + 1, j + 1]*dxmj/dx_val
                # compute second order vx velocity corrections
                if dxmj/dx_val>=0.5
                    if j < Nx_val-1
                        vxm₁₃ +=
                            0.5 *
                            ((dxmj/dx_val-0.5)^2) *
                            (vxf[i, j] - 2.0*vxf[i, j + 1] + vxf[i, j + 2])
                        vxm₂₄ +=
                            0.5 *
                            ((dxmj/dx_val-0.5)^2) *
                            (vxf[i + 1, j] - 2.0*vxf[i + 1, j + 1] + vxf[i + 1, j + 2])
                    end
                else
                    if j > 1
                        vxm₁₃ +=
                            0.5 *
                            ((dxmj/dx_val-0.5)^2) *
                            (vxf[i, j - 1] - 2.0*vxf[i, j] + vxf[i, j + 1])
                        vxm₂₄ +=
                            0.5 *
                            ((dxmj/dx_val-0.5)^2) *
                            (vxf[i + 1, j - 1] - 2.0*vxf[i + 1, j] + vxf[i + 1, j + 1])
                    end
                end
                # compute current RK step vx
                vxm[rk] = (1.0-dymi/dy_val)*vxm₁₃ + (dymi/dy_val)*vxm₂₄
                # interpolate vy
                i, j, dxmj, dymi = fix_distances(
                    xcur,
                    ycur,
                    xvy_val,
                    yvy_val,
                    dx_val,
                    dy_val,
                    jmin_vy_val,
                    jmax_vy_val,
                    imin_vy_val,
                    imax_vy_val,
                )
                # compute vy velocity for top and bottom of current cell
                vym₁₂ = vyf[i, j]*(1.0-dymi/dy_val) + vyf[i + 1, j]*dymi/dy_val
                vym₃₄ = vyf[i, j + 1]*(1.0-dymi/dy_val) + vyf[i + 1, j + 1]*dymi/dy_val
                # compute second order vy velocity corrections
                if dymi/dy_val >= 0.5
                    if i < Ny_val-1
                        vym₁₂ +=
                            0.5 *
                            ((dymi/dy_val-0.5)^2) *
                            (vyf[i, j] - 2.0*vyf[i + 1, j] + vyf[i + 2, j])
                        vym₃₄ +=
                            0.5 *
                            ((dymi/dy_val-0.5)^2) *
                            (vyf[i, j + 1] - 2.0*vyf[i + 1, j + 1] + vyf[i + 2, j + 1])
                    end
                else
                    if i > 1
                        vym₁₂ +=
                            0.5 *
                            ((dymi/dy_val-0.5)^2) *
                            (vyf[i - 1, j] - 2.0*vyf[i, j] + vyf[i + 1, j])
                        vym₃₄ +=
                            0.5 *
                            ((dymi/dy_val-0.5)^2) *
                            (vyf[i - 1, j + 1] - 2.0*vyf[i, j + 1] + vyf[i + 1, j + 1])
                    end
                end
                # compute current RK step vy
                vym[rk] = (1.0-dxmj/dx_val)*vym₁₂ + (dxmj/dx_val)*vym₃₄
                # calculate next RK step x and y positions if not at final 
                if rk==1 || rk==2
                    xcur = xA - 0.5*dt*vxm[rk]
                    ycur = yA - 0.5*dt*vym[rk]
                elseif rk==3
                    xcur = xA - dt*vxm[rk]
                    ycur = yA - dt*vym[rk]
                end
            end
            # compute effective velocity using RK4
            xcur = xA - dt*1//6*(vxm[1]+2*vxm[2]+2*vxm[3]+vxm[4])
            ycur = yA - dt*1//6*(vym[1]+2*vym[2]+2*vym[3]+vym[4])
            i, j, weights = fix_weights(
                xcur,
                ycur,
                xp_val,
                yp_val,
                dx_val,
                dy_val,
                jmin_p_val,
                jmax_p_val,
                imin_p_val,
                imax_p_val,
            )
            pf0[ii, jj] = dot4(grid_vector(i, j, pf), weights)
        end
    end # inbounds
    # end # @timeit to "backtrace_pressures_rk4!"
    return nothing
end # function backtrace_pressures_rk4!

"""
Update marker population geometry status given a marker number and 
nearest top/left marker grid point.

$(SIGNATURES)

# Details

    - m: marker number
    - i: vertical index of top/left marker grid point
    - j: horizontal index of top/left marker grid point
    - xm: horizontal x-position of markers
	- ym: vertical y-position of markers
    - mdis: minimum distance of marker launch anchor points to nearest marker
    - mnum: number of marker nearest to marker launch anchor positions

# Returns

    - nothing
"""
function update_marker_population_geometry!(
    m, i, j, xm, ym, mdis, mnum, xxm_axis=xxm, yym_axis=yym
)
    # @timeit to "update_marker_population_geometry!" begin
    @inbounds begin
        dismij = distance(xm[m], ym[m], xxm_axis[j], yym_axis[i])
        dismi1j = distance(xm[m], ym[m], xxm_axis[j], yym_axis[i + 1])
        dismij1 = distance(xm[m], ym[m], xxm_axis[j + 1], yym_axis[i])
        dismi1j1 = distance(xm[m], ym[m], xxm_axis[j + 1], yym_axis[i + 1])
        if dismij < mdis[i, j]
            mdis[i, j] = dismij
            mnum[i, j] = m
        end
        if dismi1j < mdis[i + 1, j]
            mdis[i + 1, j] = dismi1j
            mnum[i + 1, j] = m
        end
        if dismij1 < mdis[i, j + 1]
            mdis[i, j + 1] = dismij1
            mnum[i, j + 1] = m
        end
        if dismi1j1 < mdis[i + 1, j + 1]
            mdis[i + 1, j + 1] = dismi1j1
            mnum[i + 1, j + 1] = m
        end
    end # @inbounds
    # end # @timeit to "update_marker_population_geometry!"
    return nothing
end

"""
Add markers to populate currently sparsely filled grid areas.

$(SIGNATURES)

# Details

    - xm: horizontal x-position of markers
	- ym: vertical y-position of markers
	- tm: type of markers
	- tkm: temperature of markers 
	- phim: porosity of markers 
	- sxxm: marker σ′xx of markers
	- sxym: σxy of markers
	- etavpm: viscoplastic viscosity of markers
    - phinewm: reacted porosity of markers
    - pfm0: previous fluid pressure of markers
    - XWsolidm: melt molar fraction of markers
    - XWsolidm0: previous melt molar fraction of markers
    - rhototalm: total density of markers
    - rhocptotalm : total volumetric heat capacity of markers
    - etatotalm: total viscosity of markers
    - hrtotalm: total radiogenic heat production of markers
    - ktotalm: total thermal conductivity of markers
    - inv_gggtotalm: inverse of total shear modulus of markers
    - fricttotalm: total friction coefficient of markers
    - cohestotalm: total compressive strength of markers
    - tenstotalm: total tensile strength of markers
    - rhofluidcur: fluid density of markers
    - alphasolidcur: solid thermal expansion coefficient of markers
    - alphafluidcur: fluid thermal expansion coefficient of markers
    - tkm_rhocptotalm: total thermal energy of markers
    - etafluidcur_inv_kphim: fluid viscosity over permeability of markers
	- mdis: minimum distance of marker launch anchor points to nearest marker 
	- mnum: number of marker nearest to marker launch anchor positions

# Returns

    - marknum: updated number of markers in use
"""
function replenish_markers!(
    xm,
    ym,
    tm,
    tkm,
    phim,
    sxxm,
    sxym,
    etavpm,
    phinewm,
    pfm0,
    XWsolidm,
    XWsolidm0,
    rhototalm,
    rhocptotalm,
    etatotalm,
    hrtotalm,
    ktotalm,
    inv_gggtotalm,
    fricttotalm,
    cohestotalm,
    tenstotalm,
    rhofluidcur,
    alphasolidcur,
    alphafluidcur,
    tkm_rhocptotalm,
    etafluidcur_inv_kphim,
    mdis,
    mnum;
    randomized=random_markers,
    coords::Union{Nothing,GridCoordinates}=nothing,
)
    # @timeit to "replenish_markers!" begin
    Nym_val, Nxm_val = size(mnum)
    xxm_val = coords === nothing ? xxm : coords.xxm
    yym_val = coords === nothing ? yym : coords.yym
    dxm_val = coords === nothing ? dxm : coords.dxm
    dym_val = coords === nothing ? dym : coords.dym
    jmin_m_val = coords === nothing ? jmin_m : coords.jmin_m
    jmax_m_val = coords === nothing ? jmax_m : coords.jmax_m
    imin_m_val = coords === nothing ? imin_m : coords.imin_m
    imax_m_val = coords === nothing ? imax_m : coords.imax_m
    Nxmc_val = coords === nothing ? Nxmc : coords.Nxmc
    Nymc_val = coords === nothing ? Nymc : coords.Nymc
    # reset marker population geometry tracker
    mdis .= mdis_init
    mnum .= 0
    # establish marker distribution
    @inbounds begin
        for m in 1:1:length(xm)
            i, j = fix(
                xm[m],
                ym[m],
                xxm_val,
                yym_val,
                dxm_val,
                dym_val,
                jmin_m_val,
                jmax_m_val,
                imin_m_val,
                imax_m_val,
            )
            update_marker_population_geometry!(
                m, i, j, xm, ym, mdis, mnum, xxm_val, yym_val
            )
        end
        dii = 5 * Nymc_val
        djj = 5 * Nxmc_val
        for j in 1:1:Nxm_val, i in 1:1:Nym_val
            if mnum[i, j] == 0
                for jj in max(j - djj, 1):1:min(j + djj, Nxm_val)
                    for ii in max(i - dii, 1):1:min(i + dii, Nym_val)
                        if mnum[ii, jj] > 0
                            m = mnum[ii, jj]
                            dismij = distance(xm[m], ym[m], xxm_val[j], yym_val[i])
                            if dismij < mdis[i, j]
                                mdis[i, j] = dismij
                                mnum[i, j] = -m
                            end
                        end
                    end
                end
                # add new marker            
                if mnum[i, j] < 0
                    # add marker
                    if randomized
                        # production runs
                        push!(xm, xxm_val[j] + (rand(rgen)-0.5)*dxm_val)
                        push!(ym, yym_val[i] + (rand(rgen)-0.5)*dym_val)
                    else
                        # for testing
                        push!(xm, xxm_val[j])
                        push!(ym, yym_val[i])
                    end
                    # copy marker properties
                    m = -mnum[i, j]
                    push!(tm, tm[m])
                    push!(tkm, tkm[m])
                    push!(phim, phim[m])
                    push!(sxxm, sxxm[m])
                    push!(sxym, sxym[m])
                    push!(etavpm, etavpm[m])
                    push!(phinewm, phinewm[m])
                    push!(pfm0, pfm0[m])
                    push!(XWsolidm, XWsolidm[m])
                    push!(XWsolidm0, XWsolidm0[m])
                    push!(rhototalm, rhototalm[m])
                    push!(rhocptotalm, rhocptotalm[m])
                    push!(etatotalm, etatotalm[m])
                    push!(hrtotalm, hrtotalm[m])
                    push!(ktotalm, ktotalm[m])
                    push!(inv_gggtotalm, inv_gggtotalm[m])
                    push!(fricttotalm, fricttotalm[m])
                    push!(cohestotalm, cohestotalm[m])
                    push!(tenstotalm, tenstotalm[m])
                    push!(rhofluidcur, rhofluidcur[m])
                    push!(alphasolidcur, alphasolidcur[m])
                    push!(alphafluidcur, alphafluidcur[m])
                    push!(tkm_rhocptotalm, tkm_rhocptotalm[m])
                    push!(etafluidcur_inv_kphim, etafluidcur_inv_kphim[m])
                end
            end
        end
    end # @inbounds  
    return length(xm)
    # end # @timeit to "replenish_markers!"
end # function replenish_markers!

"""
Apply xy (basic grid) and xx (P grid) subgrid stress diffusion to markers.

$(SIGNATURES)

# Details

    - xm: marker x-coordinates
    - ym: marker y-coordinates
    - tm: marker type
    - inv_gggtotalm: inverse of total shear modulus of markers
    - sxxm: marker σ′xx [Pa]
    - sxym: marker σxy [Pa]
    - SXX0: σ₀′xx at P nodes [1/s]
    - SXY0: σ₀xy at basic nodes [1/s]
    - DSXX: stress change Δσ′xx at P nodes
    - DSXY: stress change Δσxy at basic nodes
    - SXXSUM: interpolation of SXX at P nodes
    - SXYSUM: interpolation of SXY at basic nodes
    - WTPSUM: interpolation weights at P nodes
    - WTSUM: interpolation weights at basic nodes
    - dt: computational time step
    - marknum: total number of markers in use

# Returns
    
    - nothing
"""
function apply_subgrid_stress_diffusion!(
    xm,
    ym,
    tm,
    inv_gggtotalm,
    sxxm,
    sxym,
    SXX0,
    SXY0,
    DSXX,
    DSXY,
    SXXSUM,
    SXYSUM,
    WTPSUM,
    WTSUM,
    dt,
    marknum;
    coords::Union{Nothing,GridCoordinates}=nothing,
)
    # @timeit to "apply_subgrid_stress_diffusion!" begin
    # only perform subgrid stress diffusion if enabled by dsubgrids > 0
    if dsubgrids == 0.0
        return nothing
    end
    Ny, Nx = size(SXYSUM)
    xp_val = coords === nothing ? xp : coords.xp
    yp_val = coords === nothing ? yp : coords.yp
    x_val = coords === nothing ? x : coords.x
    y_val = coords === nothing ? y : coords.y
    dx_val = coords === nothing ? dx : coords.dx
    dy_val = coords === nothing ? dy : coords.dy
    jmin_p_val = coords === nothing ? jmin_p : coords.jmin_p
    jmax_p_val = coords === nothing ? jmax_p : coords.jmax_p
    imin_p_val = coords === nothing ? imin_p : coords.imin_p
    imax_p_val = coords === nothing ? imax_p : coords.imax_p
    jmin_basic_val = coords === nothing ? jmin_basic : coords.jmin_basic
    jmax_basic_val = coords === nothing ? jmax_basic : coords.jmax_basic
    imin_basic_val = coords === nothing ? imin_basic : coords.imin_basic
    imax_basic_val = coords === nothing ? imax_basic : coords.imax_basic
    # fix etam[tm[m]] RMK: It's a temporary fix, not yet implemented in source
    etam = @SVector ones(3)
    # reset interpolation arrays
    SXXSUM .= 0.0
    WTPSUM .= 0.0
    SXYSUM .= 0.0
    WTSUM .= 0.0
    # iterate over markers
    for m in 1:1:marknum
        @inbounds i_p, j_p, weights_p = fix_weights(
            xm[m],
            ym[m],
            xp_val,
            yp_val,
            dx_val,
            dy_val,
            jmin_p_val,
            jmax_p_val,
            imin_p_val,
            imax_p_val,
        )
        @inbounds i_basic, j_basic, weights_basic = fix_weights(
            xm[m],
            ym[m],
            x_val,
            y_val,
            dx_val,
            dy_val,
            jmin_basic_val,
            jmax_basic_val,
            imin_basic_val,
            imax_basic_val,
        )
        # σ₀′xx at P nodes
        # compute marker-node σ′xx difference
        @inbounds δσxxm₀ = sxxm[m] - dot4(grid_vector(i_p, j_p, SXX0), weights_p)
        # time-relax σ′xx difference
        @inbounds δσxxm₀ *= (exp(-dsubgrids*dt/(etam[tm[m]]*inv_gggtotalm[m])) - 1.0)
        # correct marker stress
        @inbounds sxxm[m] += δσxxm₀
        # update subgrid diffusion on P nodes
        interpolate_add_to_grid!(i_p, j_p, weights_p, δσxxm₀, SXXSUM)
        interpolate_add_to_grid!(i_p, j_p, weights_p, 1.0, WTPSUM)
        # σ₀xy at basic nodes
        # compute marker-node σxy difference
        @inbounds δσxy₀ = sxym[m] - dot4(grid_vector(i_basic, j_basic, SXY0), weights_basic)
        # time-relax σxy difference
        @inbounds δσxy₀ *= (exp(-dsubgrids*dt/(etam[tm[m]]*inv_gggtotalm[m])) - 1.0)
        # correct marker stress
        @inbounds sxym[m] += δσxy₀
        # update subgrid diffusion on basic nodes
        interpolate_add_to_grid!(i_basic, j_basic, weights_basic, δσxy₀, SXYSUM)
        interpolate_add_to_grid!(i_basic, j_basic, weights_basic, 1.0, WTSUM)
    end
    # compute DSXXsubgrid and update DSXX at inner P nodes
    @views @. DSXX[2:Ny, 2:Nx][WTPSUM[2:Ny, 2:Nx] > 0.0] -= (
        SXXSUM[2:Ny, 2:Nx][WTPSUM[2:Ny, 2:Nx] > 0.0] /
        WTPSUM[2:Ny, 2:Nx][WTPSUM[2:Ny, 2:Nx] > 0.0]
    )
    # compute DSXYsubgrid and update DSXY at all basic nodes
    @views @. DSXY[WTSUM[:, :] > 0.0] -=
        SXYSUM[:, :][WTSUM[:, :] > 0.0] / WTSUM[:, :][WTSUM[:, :] > 0.0]
    # end # @timeit to "apply_subgrid_stress_diffusion!"
    return nothing
end # function apply_subgrid_stress_diffusion!

"""
Update marker stress based on xy (basic grid) and xx (P grid) stress changes.

$(SIGNATURES)

# Details

    - xm: x-coordinates of markers
    - ym: y-coordinates of markers
    - sxxm: marker σ′xx [Pa]
    - sxym: marker σxy [Pa]
    - DSXX: stress change Δσ′xx at P nodes
    - DSXY: stress change Δσxy at basic nodes
    - marknum: total number of markers in use

# Returns

    - nothing
"""
function update_marker_stress!(
    xm, ym, sxxm, sxym, DSXX, DSXY, marknum; coords::Union{Nothing,GridCoordinates}=nothing
)
    # @timeit to "update_marker_stress!" begin
    Ny, Nx = size(DSXY)
    xp_val = coords === nothing ? xp : coords.xp
    yp_val = coords === nothing ? yp : coords.yp
    x_val = coords === nothing ? x : coords.x
    y_val = coords === nothing ? y : coords.y
    dx_val = coords === nothing ? dx : coords.dx
    dy_val = coords === nothing ? dy : coords.dy
    @threads :static for m in 1:1:marknum
        @inbounds i_p, j_p, weights_p = fix_weights(
            xm[m], ym[m], xp_val, yp_val, dx_val, dy_val, 2, Nx-1, 2, Ny-1
        )
        @inbounds i_basic, j_basic, weights_basic = fix_weights(
            xm[m], ym[m], x_val, y_val, dx_val, dy_val, 1, Nx-1, 1, Ny-1
        )
        # interpolate updated DSXX, DSXY back to markers
        interpolate_add_to_marker!(m, i_p, j_p, weights_p, sxxm, DSXX)
        interpolate_add_to_marker!(m, i_basic, j_basic, weights_basic, sxym, DSXY)
    end
    # end # @timeit to "update_marker_stress!"
    return nothing
end # function update_marker_stress!

"""
Apply subgrid temperature diffusion to markers.

$(SIGNATURES)

# Details

    - xm: marker x-coordinates 
	- ym: marker y-coordinates
	- tm: marker type
	- tkm: marker temperature
    - phim: marker porosity
    - tk1: current temperature at P nodes
	- DT: temperature change at P nodes
	- TKSUM: interpolation of TK at P nodes 
	- RHOCPSUM: interpolation of RHOCP at P nodes
    - dt: computational time step
	- marknum: total number of markers in use
    - mode: marker property computation mode
        - 1: dynamic, based on (Touloukian, 1970; Hobbs, 1974;
             Travis and Schubert, 2005)
        - 2: constant parameter rhocpfluidm

# Returns

    - nothing
"""
function apply_subgrid_temperature_diffusion!(
    xm,
    ym,
    tm,
    tkm,
    phim,
    tk1,
    DT,
    TKSUM,
    RHOCPSUM,
    dt,
    marknum,
    mode;
    coords::Union{Nothing,GridCoordinates}=nothing,
)
    # @timeit to "apply_subgrid_temperature_diffusion!" begin
    # only perform subgrid temperature diffusion if enabled by dsubgridt > 0
    if dsubgridt == 0.0
        return nothing
    end
    Ny1, Nx1 = size(RHOCPSUM)
    xp_val = coords === nothing ? xp : coords.xp
    yp_val = coords === nothing ? yp : coords.yp
    dx_val = coords === nothing ? dx : coords.dx
    dy_val = coords === nothing ? dy : coords.dy
    jmin_p_val = coords === nothing ? jmin_p : coords.jmin_p
    jmax_p_val = coords === nothing ? jmax_p : coords.jmax_p
    imin_p_val = coords === nothing ? imin_p : coords.imin_p
    imax_p_val = coords === nothing ? imax_p : coords.imax_p
    # reset interpolation arrays
    TKSUM .= 0.0
    RHOCPSUM .= 0.0
    # iterate over markers
    for m in 1:1:marknum
        @inbounds i, j, weights = fix_weights(
            xm[m],
            ym[m],
            xp_val,
            yp_val,
            dx_val,
            dy_val,
            jmin_p_val,
            jmax_p_val,
            imin_p_val,
            imax_p_val,
        )
        # compute marker-node temperature difference
        @inbounds δtkm = tkm[m] - dot4(grid_vector(i, j, tk1), weights)
        # compute marker properties
        @inbounds if tm[m] < 3
            # rocks
            @inbounds rhocptotalm = total(
                rhocpsolidm[tm[m]], compute_rhocpfluidm(tkm[m], mode), phim[m]
            )
            @inbounds ktotalm[m] = ktotal(
                compute_ksolidm(tkm[m], mode), compute_kfluidm(tkm[m], mode), phim[m]
            )
        else
            # sticky air
            @inbounds rhocptotalm = rhocpsolidm[tm[m]]
            @inbounds ktotalm = ksolidm[tm[m]]
        end
        # time-relax δtkm difference
        δtkm *= (
            exp(-dsubgridt*ktotalm*dt/rhocptotalm*(2.0*(inv(dx_val^2)+inv(dy_val^2)))) - 1.0
        )
        # correct marker temperature
        @inbounds tkm[m] += δtkm
        # update subgrid diffusion on P nodes
        interpolate_add_to_grid!(i, j, weights, δtkm*rhocptotalm, TKSUM)
        interpolate_add_to_grid!(i, j, weights, rhocptotalm, RHOCPSUM)
    end
    # compute DTsubgrid=TKSUM/RHOCPSUM and update temperature field at P nodes
    for j in 1:1:Nx1, i in 1:1:Ny1
        if RHOCPSUM[i, j] > 0.0
            DT[i, j] -= TKSUM[i, j] / RHOCPSUM[i, j]
        end
    end
    # end # @timeit to "apply_subgrid_temperature_diffusion!"
    return nothing
end # function apply_subgrid_temperature_diffusion! 

"""
Update marker temperature based on P grid temperature changes.

$(SIGNATURES)

# Details

    - xm: x-coordinates of markers
    - ym: y-coordinates of markers
    - tkm: marker temperature
    - DT: temperature change at P nodes
    - tk2: next temperature at P nodes
    - timestep: current time step
    - marknum: total number of markers in use

# Returns

    - nothing
# """
function update_marker_temperature!(
    xm, ym, tkm, DT, tk2, timestep, marknum; coords::Union{Nothing,GridCoordinates}=nothing
)
    # @timeit to "update_marker_temperature!" begin
    xp_val = coords === nothing ? xp : coords.xp
    yp_val = coords === nothing ? yp : coords.yp
    dx_val = coords === nothing ? dx : coords.dx
    dy_val = coords === nothing ? dy : coords.dy
    jmin_p_val = coords === nothing ? jmin_p : coords.jmin_p
    jmax_p_val = coords === nothing ? jmax_p : coords.jmax_p
    imin_p_val = coords === nothing ? imin_p : coords.imin_p
    imax_p_val = coords === nothing ? imax_p : coords.imax_p
    if timestep == 1
        # interpolate tk2 to markers instead of DT for first time step        
        @threads :static for m in 1:1:marknum
            @inbounds i, j, weights = fix_weights(
                xm[m],
                ym[m],
                xp_val,
                yp_val,
                dx_val,
                dy_val,
                jmin_p_val,
                jmax_p_val,
                imin_p_val,
                imax_p_val,
            )
            interpolate_to_marker!(m, i, j, weights, tkm, tk2)
        end
    else
        # interpolate and apply DT to markers for subsequent time steps
        @threads :static for m in 1:1:marknum
            @inbounds i, j, weights = fix_weights(
                xm[m],
                ym[m],
                xp_val,
                yp_val,
                dx_val,
                dy_val,
                jmin_p_val,
                jmax_p_val,
                imin_p_val,
                imax_p_val,
            )
            interpolate_add_to_marker!(m, i, j, weights, tkm, DT)
        end
    end
    # end # @timeit to "update_marker_temperature!"
    return nothing
end # function update_marker_temperature!

"""
Update marker porosity for compaction based on Dln[(1-ϕ)/ϕ]/Dt at P grid.

$(SIGNATURES)

# Details

    - xm: x-coordinates of markers
    - ym: y-coordinates of markers
    - tm: marker type
    - phim: marker porosity
    - APHI: Dln[(1-ϕ)/ϕ]/Dt at P nodes
    - dt: computational time step
    - marknum: total number of markers in use

# Returns

    - nothing
# """
function update_marker_porosity!(
    xm,
    ym,
    tm,
    phim,
    APHI,
    dt,
    marknum;
    phimin=phimin,
    phimax=phimax,
    coords::Union{Nothing,GridCoordinates}=nothing,
)
    # @timeit to "update_marker_porosity!" begin
    # update porosity for compaction
    xp_val = coords === nothing ? xp : coords.xp
    yp_val = coords === nothing ? yp : coords.yp
    dx_val = coords === nothing ? dx : coords.dx
    dy_val = coords === nothing ? dy : coords.dy
    jmin_p_val = coords === nothing ? jmin_p : coords.jmin_p
    jmax_p_val = coords === nothing ? jmax_p : coords.jmax_p
    imin_p_val = coords === nothing ? imin_p : coords.imin_p
    imax_p_val = coords === nothing ? imax_p : coords.imax_p
    @inbounds begin
        @threads :static for m in 1:1:marknum
            if tm[m] < 3
                # rocks
                i, j, weights = fix_weights(
                    xm[m],
                    ym[m],
                    xp_val,
                    yp_val,
                    dx_val,
                    dy_val,
                    jmin_p_val,
                    jmax_p_val,
                    imin_p_val,
                    imax_p_val,
                )
                # compute Dln[(1-ϕ)/ϕ]/Dt at marker
                aphim = dot4(grid_vector(i, j, APHI), weights)
                # update marker porosity
                phim[m] = max(
                    phimin, min(phimax, phim[m] / ((1.0-phim[m])*exp(aphim*dt) + phim[m]))
                )
            end
        end
    end # @inbounds
    # end # @timeit to "update_marker_porosity!"
    return nothing
end # function update_marker_porosity!
