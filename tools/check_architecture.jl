#!/usr/bin/env julia
# AST-based architecture ratchet verifying global bindings, argument budgets,
# RNG safety, coords guards, and function line spans against baseline.

const ROOT_DIR = normpath(joinpath(@__DIR__, ".."))
if Base.find_package("JSON") === nothing
    pushfirst!(LOAD_PATH, ROOT_DIR)
end

using JSON

const SRC_DIR = joinpath(ROOT_DIR, "src")
const BASELINE_PATH = joinpath(@__DIR__, "architecture_baseline.json")
const ALLOWLIST_PATH = joinpath(@__DIR__, "architecture_allowlist.txt")

const RNG_SERIAL_ALLOWLIST = Set([
    "setup_staggered_grid_properties",
    "setup_staggered_grid_properties_helpers",
    "setup_marker_properties",
    "setup_marker_properties_helpers",
    "define_markers!",
    "replenish_markers!",
    "sample_parameters",
])

# Mutable types forbidden in const bindings unless in allowlist
const FORBIDDEN_MUTABLE_CONSTRUCTORS = Set([
    :MersenneTwister, :TimerOutput, :zeros, :fill, :Array, :Dict, :Ref, :Vector, :Matrix
])

function parse_all_expressions(content::String)
    exprs = Tuple{Any,Int}[]
    pos = 1
    len = ncodeunits(content)
    line_num = 1
    while pos <= len
        # Skip whitespace and comments
        while pos <= len
            b = codeunit(content, pos)
            if b == UInt8(' ') || b == UInt8('\t') || b == UInt8('\r') || b == UInt8('\n')
                if b == UInt8('\n')
                    line_num += 1
                end
                pos += 1
            elseif b == UInt8('#')
                if pos < len && codeunit(content, pos + 1) == UInt8('=')
                    pos += 2
                    while pos < len && !(
                        codeunit(content, pos) == UInt8('=') &&
                        codeunit(content, pos + 1) == UInt8('#')
                    )
                        if codeunit(content, pos) == UInt8('\n')
                            line_num += 1
                        end
                        pos += 1
                    end
                    pos = min(len + 1, pos + 2)
                else
                    while pos <= len && codeunit(content, pos) != UInt8('\n')
                        pos += 1
                    end
                    if pos <= len && codeunit(content, pos) == UInt8('\n')
                        line_num += 1
                        pos += 1
                    end
                end
            else
                break
            end
        end
        pos > len && break
        start_line = line_num
        ex, next_pos = Meta.parse(content, pos)
        ex === nothing && break
        push!(exprs, (ex, start_line))
        for i in pos:(next_pos - 1)
            if codeunit(content, i) == UInt8('\n')
                line_num += 1
            end
        end
        pos = next_pos
    end
    return exprs
end

function flatten_toplevel(ex, line::Int)
    out = Tuple{Any,Int}[]
    if Meta.isexpr(ex, :toplevel) || (Meta.isexpr(ex, :block) && !isempty(ex.args))
        curr_line = line
        for a in ex.args
            if a isa LineNumberNode
                curr_line = a.line
            else
                append!(out, flatten_toplevel(a, curr_line))
            end
        end
    else
        push!(out, (ex, line))
    end
    return out
end

function is_function_assignment(lhs)
    cur = lhs
    while Meta.isexpr(cur, :where) || Meta.isexpr(cur, :(::))
        cur = cur.args[1]
    end
    return Meta.isexpr(cur, :call)
end

function read_allowlist(allowlist_path=ALLOWLIST_PATH)
    isfile(allowlist_path) || return Set{String}()
    allowed = Set{String}()
    for line in eachline(allowlist_path)
        s = strip(line)
        isempty(s) && continue
        startswith(s, "#") && continue
        push!(allowed, s)
    end
    return allowed
end

function get_src_files(src_dir=SRC_DIR)
    files = String[]
    for (root, _, fnames) in walkdir(src_dir)
        for fn in fnames
            endswith(fn, ".jl") || continue
            push!(files, joinpath(root, fn))
        end
    end
    return sort(files)
end

function extract_function_name(sig)
    cur = sig
    while Meta.isexpr(cur, :where) || Meta.isexpr(cur, :(::))
        cur = cur.args[1]
    end
    if cur isa Symbol
        return string(cur)
    elseif Meta.isexpr(cur, :(.))
        return string(cur.args[1], ".", cur.args[2].value)
    elseif Meta.isexpr(cur, :call)
        return extract_function_name(cur.args[1])
    end
    return "anonymous"
end

function count_positional_args(sig)
    cur = sig
    while Meta.isexpr(cur, :where) || Meta.isexpr(cur, :(::))
        cur = cur.args[1]
    end
    if Meta.isexpr(cur, :call)
        count = 0
        for i in 2:length(cur.args)
            arg = cur.args[i]
            # Keyword arguments are in a :parameters expr as first arg in call_expr.args
            if Meta.isexpr(arg, :parameters)
                continue
            end
            count += 1
        end
        return count
    end
    return 0
end

function check_expr_for_threading(ex)::Bool
    found = false
    function walk(e)
        if Meta.isexpr(e, :macrocall) && length(e.args) >= 1
            mname = string(e.args[1])
            if occursin("threads", mname) || occursin("spawn", mname)
                found = true
            end
        end
        if isa(e, Expr)
            for a in e.args
                walk(a)
            end
        end
    end
    walk(ex)
    return found
end

function check_expr_for_rng(ex)::Bool
    found = false
    function walk(e)
        if Meta.isexpr(e, :call) && length(e.args) >= 1
            fn = e.args[1]
            if fn in (:rand, :randn, :(Random.rand), :(Random.randn))
                found = true
            end
        end
        if isa(e, Expr)
            for a in e.args
                walk(a)
            end
        end
    end
    walk(ex)
    return found
end

function is_mutable_construction(rhs)::Bool
    if Meta.isexpr(rhs, :call) && length(rhs.args) >= 1
        fn = rhs.args[1]
        base_fn = fn isa Symbol ? fn : (Meta.isexpr(fn, :(.)) ? fn.args[2].value : nothing)
        if base_fn in FORBIDDEN_MUTABLE_CONSTRUCTORS
            return true
        end
    end
    return false
end

function analyze_codebase(src_dir=SRC_DIR)
    files = get_src_files(src_dir)

    global_bindings = Tuple{String,String,Int,Bool}[] # (name, file, line, is_mutable)
    function_arg_counts = Dict{String,Int}()
    function_line_spans = Dict{String,Int}()
    includes = String[]
    coords_guard_count = 0
    rng_violations = String[]

    for fpath in files
        rel_path = relpath(fpath, ROOT_DIR)
        content = read(fpath, String)
        lines = split(content, '\n')

        # 1. Coords guards in text (excluding comments)
        for (line_no, line) in enumerate(lines)
            stripped = strip(line)
            startswith(stripped, "#") && continue
            if occursin("coords === nothing", line) ||
                occursin("coords !== nothing", line) ||
                occursin("isnothing(coords)", line) ||
                occursin("coords == nothing", line)
                coords_guard_count += 1
            end
        end

        # 2. AST-based checks
        exprs = parse_all_expressions(content)
        for (top_ex, top_line) in exprs
            for (ex, start_line) in flatten_toplevel(top_ex, top_line)
                # Top-level includes
                if Meta.isexpr(ex, :call) && length(ex.args) >= 2 && ex.args[1] === :include
                    push!(includes, string(ex.args[2]))
                end

                # Check global / const bindings
                if Meta.isexpr(ex, :const)
                    # const binding
                    inner = ex.args[1]
                    if Meta.isexpr(inner, :(=))
                        lhs = inner.args[1]
                        rhs = inner.args[2]
                        cur = lhs
                        while Meta.isexpr(cur, :where) || Meta.isexpr(cur, :(::))
                            cur = cur.args[1]
                        end
                        vname = cur isa Symbol ? string(cur) : string(lhs)
                        push!(
                            global_bindings,
                            (vname, rel_path, start_line, is_mutable_construction(rhs)),
                        )
                    end
                elseif Meta.isexpr(ex, :(=))
                    lhs = ex.args[1]
                    rhs = ex.args[2]
                    if !is_function_assignment(lhs) &&
                        !Meta.isexpr(lhs, :ref) &&
                        !Meta.isexpr(lhs, :(.))
                        cur = lhs
                        while Meta.isexpr(cur, :where) || Meta.isexpr(cur, :(::))
                            cur = cur.args[1]
                        end
                        vname = cur isa Symbol ? string(cur) : string(lhs)
                        push!(
                            global_bindings,
                            (vname, rel_path, start_line, is_mutable_construction(rhs)),
                        )
                    end
                end

                # Walk all expressions inside ex to find functions and rand calls
                function walk_ast(e, parent_fn=nothing, enclosing_fn_expr=nothing)
                    if Meta.isexpr(e, :function) ||
                        (Meta.isexpr(e, :(=)) && is_function_assignment(e.args[1]))
                        sig = e.args[1]
                        body = e.args[2]
                        fn_name = extract_function_name(sig)
                        p_args = count_positional_args(sig)

                        fn_key = "$rel_path:$fn_name"
                        # Track max positional args per function
                        if !haskey(function_arg_counts, fn_name) ||
                            function_arg_counts[fn_name] < p_args
                            function_arg_counts[fn_name] = p_args
                        end

                        # Compute line span estimate
                        body_str = sprint(Base.show_unquoted, e)
                        span = count(c -> c == '\n', body_str) + 1
                        function_line_spans[fn_key] = max(
                            get(function_line_spans, fn_key, 0), span
                        )

                        # Walk body with parent function set
                        for a in e.args[2:end]
                            walk_ast(a, fn_name, e)
                        end
                        return nothing
                    end

                    # Check RNG call safety
                    if Meta.isexpr(e, :call) && length(e.args) >= 1
                        fn = e.args[1]
                        if fn in (:rand, :randn, :(Random.rand), :(Random.randn))
                            if parent_fn !== nothing
                                if !(parent_fn in RNG_SERIAL_ALLOWLIST)
                                    push!(
                                        rng_violations,
                                        "$rel_path: rand() called in unapproved function '$parent_fn'",
                                    )
                                elseif enclosing_fn_expr !== nothing &&
                                    check_expr_for_threading(enclosing_fn_expr)
                                    push!(
                                        rng_violations,
                                        "$rel_path: rand() allowlisted function '$parent_fn' uses threading",
                                    )
                                end
                            else
                                push!(
                                    rng_violations,
                                    "$rel_path:$start_line: Top-level rand() call detected",
                                )
                            end
                        end
                    end

                    if isa(e, Expr)
                        for a in e.args
                            walk_ast(a, parent_fn, enclosing_fn_expr)
                        end
                    end
                end
                walk_ast(ex)
            end
        end
    end

    sort!(includes)
    return (;
        global_bindings,
        function_arg_counts,
        function_line_spans,
        includes,
        coords_guard_count,
        rng_violations,
    )
end

function generate_baseline(src_dir=SRC_DIR, baseline_path=BASELINE_PATH)
    println("Analyzing codebase for architecture baseline...")
    analysis = analyze_codebase(src_dir)

    baseline_data = Dict{String,Any}(
        "coords_guard_count" => analysis.coords_guard_count,
        "includes" => analysis.includes,
        "function_max_positional_args" => analysis.function_arg_counts,
        "function_line_spans" => analysis.function_line_spans,
    )

    open(baseline_path, "w") do io
        JSON.print(io, baseline_data, 4)
        return println(io)
    end
    println("Architecture baseline saved to $baseline_path")
    println("  Coords guard count: $(analysis.coords_guard_count)")
    println("  Tracked functions: $(length(analysis.function_arg_counts))")
    return println("  Tracked includes: $(length(analysis.includes))")
end

function check_architecture(;
    src_dir=SRC_DIR,
    baseline_path=BASELINE_PATH,
    allowlist_path=ALLOWLIST_PATH,
    exit_on_failure=true,
)
    println("Running AST architecture ratchet checks...")
    analysis = analyze_codebase(src_dir)
    allowlist = read_allowlist(allowlist_path)

    errors = String[]

    # 1. Global bindings allowlist & mutability check
    for (name, file, line, is_mutable) in analysis.global_bindings
        if !(name in allowlist)
            push!(
                errors,
                "[GLOBAL_BINDING] $file:$line - Unapproved module binding '$name' (must be in allowlist)",
            )
        elseif is_mutable
            # If mutable, only 'to' is allowed as legacy exception
            if name != "to" &&
                name != "iparms_dict" &&
                name != "iparms" &&
                name != "SPECIATION_SPECIES" &&
                name != "SPECIES_AMU" &&
                name != "VALID_SECTIONS" &&
                name != "METAL_SILICATE_CAP_WARNING_COUNTER" &&
                name != "PICARD_WARNING_COUNTER" &&
                name != "SPECIES_MOLAR_MASS" &&
                name != "SPECIES_O_STOICH"
                push!(
                    errors,
                    "[MUTABLE_GLOBAL] $file:$line - Const binding '$name' instantiates a mutable object",
                )
            end
        end
    end

    # 2. RNG safety
    for viol in analysis.rng_violations
        push!(errors, "[RNG_THREAD_SAFETY] $viol")
    end

    # 3. Baseline comparisons
    if !isfile(baseline_path)
        push!(
            errors,
            "[BASELINE_MISSING] Architecture baseline file not found at $baseline_path. Run with --baseline.",
        )
    else
        baseline = JSON.parsefile(baseline_path)

        # Coords guards ratchet
        base_coords = get(baseline, "coords_guard_count", 0)
        if analysis.coords_guard_count > base_coords
            push!(
                errors,
                "[COORDS_GUARD_RATCHET] Coords guard count increased: $(analysis.coords_guard_count) > $base_coords",
            )
        end

        # Includes ratchet
        base_includes = Set(get(baseline, "includes", String[]))
        for inc in analysis.includes
            if !(inc in base_includes)
                push!(errors, "[INCLUDE_RATCHET] New unapproved include detected: '$inc'")
            end
        end

        # Positional arguments ratchet
        base_args = get(baseline, "function_max_positional_args", Dict{String,Any}())
        for (fname, count) in analysis.function_arg_counts
            if haskey(base_args, fname)
                bcount = base_args[fname]
                if count > bcount
                    push!(
                        errors,
                        "[ARG_BUDGET] Function '$fname' positional args increased: $count > $bcount",
                    )
                end
            else
                if count > 12
                    push!(
                        errors,
                        "[ARG_BUDGET] New function '$fname' exceeds max positional args (12): got $count",
                    )
                end
            end
        end

        # Function line spans ratchet
        base_spans = get(baseline, "function_line_spans", Dict{String,Any}())
        for (fkey, span) in analysis.function_line_spans
            if haskey(base_spans, fkey)
                bspan = base_spans[fkey]
                # Allow a tiny tolerance for formatting / comments, but ratchet growth
                if span > bspan + 25
                    push!(
                        errors,
                        "[LINE_BUDGET] Function '$fkey' line span grew significantly: $span > $bspan",
                    )
                end
            else
                if span > 400
                    push!(
                        errors,
                        "[LINE_BUDGET] New function '$fkey' exceeds line span budget (400): got $span",
                    )
                end
            end
        end
    end

    if !isempty(errors)
        println(
            stderr, "Architecture ratchet check FAILED with $(length(errors)) violations:"
        )
        for err in errors
            println(stderr, "  - $err")
        end
        if exit_on_failure
            exit(1)
        end
        return false, errors
    end

    println("Architecture ratchet check PASSED clean.")
    return true, errors
end

function main()
    mode = length(ARGS) >= 1 ? ARGS[1] : "--check"
    if mode == "--baseline"
        generate_baseline()
    elseif mode == "--check"
        check_architecture(; exit_on_failure=true)
    else
        println(stderr, "Unknown mode '$mode'. Use --baseline or --check.")
        exit(1)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
