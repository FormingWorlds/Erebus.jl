#!/usr/bin/env julia
# tools/replay_i2elvis_porosity.jl
# Replays the cold pressing and hot sintering block of the i2elvis code (markmart.c)
# at a set of marker states and prints the table that test/test_porosity_henke.jl
# compares with the Erebus functions.
#
# Usage: julia tools/replay_i2elvis_porosity.jl <path to i2elvis clone> [C compiler]
#
# The i2elvis source is GPL-3, so it is not copied into this repository. The script
# reads the lines from a local clone, wraps them in a C function and compiles it.

"""
    REPLAY_STATES

Marker states `(phi, T [K], P [Pa], phi_init, dt [s])` over every branch of the block:
cold pressing (T < 700 K), hot sintering with phi > 0.1 and with phi <= 0.1, the 0.85 cap
of the starting porosity, and the 700 K threshold. The time steps give a porosity change
of about 1e-3 of the porosity so that 1 explicit step is accurate.
"""
const REPLAY_STATES = (
    (0.80, 500.0, 5.0e4, 0.80, 1.0e10),
    (0.60, 650.0, 1.0e6, 0.75, 1.0e10),
    (0.84, 699.5, 1.0e3, 0.95, 1.0e10),
    (0.38, 900.0, 1.0e5, 0.40, 1.9e7),
    (0.08, 1000.0, 5.0e5, 0.40, 2.0e5),
    (0.30, 1100.0, 2.0e6, 0.50, 1.0e4),
    (0.84, 800.0, 1.0e3, 0.95, 1.8e14),
    (0.40, 700.0, 2.0e5, 0.42, 6.0e12),
)

"""
    extract_lines(lines, first, last, anchor_first, anchor_last) -> String

Return lines `first:last` of `lines` joined with newlines, after checking that the first
and last lines contain the expected anchor text.

# Throws
- `ErrorException` if an anchor is missing (the clone differs from the replayed version).
"""
function extract_lines(lines, first, last, anchor_first, anchor_last)
    occursin(anchor_first, lines[first]) ||
        error("markmart.c line $first does not contain '$anchor_first'")
    occursin(anchor_last, lines[last]) ||
        error("markmart.c line $last does not contain '$anchor_last'")
    return join(lines[first:last], "\n")
end

"""
    build_wrapper(markmart_path) -> String

Build the C source of the replay: the local declarations of `ronurecalc()`, the 2 lines
that set the marker pressure and temperature, and the porosity block, verbatim. The
globals of the block (`markpor`, `markk`, `eps`, `timestep`, `por_init`) are set from
the arguments of the wrapper.
"""
function build_wrapper(markmart_path::AbstractString)
    lines = readlines(markmart_path)
    decls = extract_lines(lines, 195, 222, "long int m1", "double mpb=1.0")
    pt = extract_lines(lines, 282, 283, "mpb=eps[10]", "mtk=(double)(markk[mm1])")
    block = extract_lines(lines, 362, 472, "Effect of cold pressing", "markpor[mm1]>1.000")
    return """
    #include <math.h>
    #include <stdio.h>
    #define MINV(a,b) ((a)<=(b)? (a):(b))
    #define ABSV(a)   ((a)>=0? (a):-(a))
    static float markpor[1], markk[1];
    static double eps[500], timestep, por_init;
    static void replay(double phi_in, double t_in, double p_in, double phi_init_in,
                       double dt_in, double *phi_out, double *change_out)
    {
    $decls
    markpor[0] = (float)phi_in; markk[0] = (float)t_in; eps[10] = p_in;
    por_init = phi_init_in; timestep = dt_in; mm1 = 0; mm2 = 5; por_change = 0.0;
    $pt
    $block
    *phi_out = (double)markpor[0]; *change_out = por_change;
    }
    int main(void)
    {
    double a, b, c, d, e, phi, change;
    while (scanf("%lf %lf %lf %lf %lf", &a, &b, &c, &d, &e) == 5) {
        replay(a, b, c, d, e, &phi, &change);
        printf("%.17e %.17e\\n", phi, change);
    }
    return 0;
    }
    """
end

"""
    run_replay(i2elvis_dir; cc="cc") -> Vector{NTuple{7,Float64}}

Compile the wrapper and run it at `REPLAY_STATES`. Returns rows
`(phi, T, P, phi_init, dt, phi_new, por_change)`.
"""
function run_replay(i2elvis_dir::AbstractString; cc::AbstractString="cc")
    src = build_wrapper(joinpath(i2elvis_dir, "markmart.c"))
    dir = mktempdir()
    cpath = joinpath(dir, "replay.c")
    exe = joinpath(dir, "replay")
    write(cpath, src)
    run(`$cc -O0 -o $exe $cpath -lm`)
    input = join([join(string.(s), " ") for s in REPLAY_STATES], "\n") * "\n"
    out = read(pipeline(`$exe`; stdin=IOBuffer(input)), String)
    rows = NTuple{7,Float64}[]
    for (s, line) in zip(REPLAY_STATES, split(strip(out), "\n"))
        vals = parse.(Float64, split(line))
        push!(rows, (s..., vals[1], vals[2]))
    end
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) &&
        error("usage: julia tools/replay_i2elvis_porosity.jl <i2elvis dir> [cc]")
    cc = length(ARGS) >= 2 ? ARGS[2] : get(ENV, "CC", "cc")
    rows = run_replay(ARGS[1]; cc=cc)
    println("# (phi, T [K], P [Pa], phi_init, dt [s], phi_new, por_change)")
    println("const _HK_I2ELVIS_REPLAY = (")
    for r in rows
        println("    (", join([repr(x) for x in r], ", "), "),")
    end
    println(")")
end
