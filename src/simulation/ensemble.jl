"""
Ensemble simulation sweep and parameter exploration framework.

Provides structured parameter sampling (grid, Latin Hypercube, random),
perturbed configuration generation, and batch execution with scalar telemetry
and catalog compilation.
"""

"""
Specification for an ensemble parameter sweep.

$(FIELDS)
"""
struct EnsembleSweepSpec
    base_config::SimulationConfig
    output_dir::String
    sampling_method::Symbol
    num_samples::Int
    parameters::Dict{String,Any}
    seed::Int

    function EnsembleSweepSpec(
        base_config::SimulationConfig,
        output_dir::AbstractString,
        sampling_method::Symbol,
        num_samples::Integer,
        parameters::AbstractDict,
        seed::Integer,
    )
        sampling_method in (:grid, :lhs, :random) || throw(
            ArgumentError(
                "sampling_method must be one of :grid, :lhs, :random, got :$sampling_method",
            ),
        )
        num_samples >= 1 ||
            throw(ArgumentError("num_samples must be >= 1, got $num_samples"))
        return new(
            base_config,
            String(output_dir),
            sampling_method,
            Int(num_samples),
            Dict{String,Any}(String(k) => v for (k, v) in parameters),
            Int(seed),
        )
    end
end

"""
    EnsembleSweepSpec(base_config; output_dir="ensemble_output", sampling_method=:lhs,
                      num_samples=10, parameters=Dict(), seed=42)

Construct an ensemble sweep specification.
"""
function EnsembleSweepSpec(
    base_config::SimulationConfig;
    output_dir::AbstractString="ensemble_output",
    sampling_method::Symbol=:lhs,
    num_samples::Integer=10,
    parameters::AbstractDict=Dict{String,Any}(),
    seed::Integer=42,
)
    return EnsembleSweepSpec(
        base_config, output_dir, sampling_method, num_samples, parameters, seed
    )
end

"""
    load_ensemble_config(source::AbstractString)::EnsembleSweepSpec

Load an ensemble sweep specification from a TOML file or string.
"""
function load_ensemble_config(source::AbstractString)::EnsembleSweepSpec
    parsed = if isfile(source)
        TOML.parsefile(source)
    else
        TOML.parse(source)
    end

    haskey(parsed, "base") || throw(
        ArgumentError("Missing required [base] section in ensemble configuration"),
    )
    base_sec = parsed["base"]
    haskey(base_sec, "config") || throw(
        ArgumentError("Missing required 'config' key in [base] section of ensemble configuration"),
    )

    base_path_raw = String(base_sec["config"])
    base_cfg_path = if isfile(source) && !isabspath(base_path_raw)
        cand1 = joinpath(dirname(source), base_path_raw)
        if isfile(cand1)
            cand1
        elseif isfile(base_path_raw)
            base_path_raw
        else
            cand_pkg = joinpath(dirname(dirname(@__DIR__)), base_path_raw)
            isfile(cand_pkg) ? cand_pkg : cand1
        end
    else
        if isfile(base_path_raw)
            base_path_raw
        else
            cand_pkg = joinpath(dirname(dirname(@__DIR__)), base_path_raw)
            isfile(cand_pkg) ? cand_pkg : base_path_raw
        end
    end

    isfile(base_cfg_path) || throw(
        ArgumentError("Base configuration file not found: '$base_path_raw' (resolved to '$base_cfg_path')"),
    )

    base_cfg = load_config(base_cfg_path)

    ens_sec = get(parsed, "ensemble", Dict{String,Any}())
    output_dir = get(ens_sec, "output_dir", "ensemble_output")
    method_str = get(ens_sec, "method", "lhs")
    method = Symbol(method_str)
    num_samples = get(ens_sec, "num_samples", 10)
    seed = get(ens_sec, "seed", 42)

    params = get(parsed, "parameters", Dict{String,Any}())

    return EnsembleSweepSpec(
        base_cfg;
        output_dir=output_dir,
        sampling_method=method,
        num_samples=num_samples,
        parameters=params,
        seed=seed,
    )
end

"""
    override_config(cfg, overrides)

Return a new `SimulationConfig` with nested field overrides applied.
Each key in `overrides` must follow the `"section.field"` naming convention.
"""
function override_config(
    cfg::SimulationConfig, overrides::AbstractDict{<:AbstractString,<:Any}
)
    d = config_to_dict(cfg)

    for (k, v) in overrides
        parts = split(String(k), ".")
        length(parts) == 2 ||
            throw(ArgumentError("Override key must be in format 'section.field', got '$k'"))
        sec, fld = parts
        haskey(d, sec) || throw(ArgumentError("Unknown section '$sec' in override '$k'"))
        haskey(d[sec], fld) ||
            throw(ArgumentError("Unknown field '$fld' in section '$sec'"))
        d[sec][fld] = v isa Symbol ? String(v) : v
    end

    io = IOBuffer()
    TOML.print(x -> x isa Symbol ? String(x) : nothing, io, d)
    seekstart(io)
    return parse_config_string(String(take!(io)))
end

"""
Determines whether a sweep parameter corresponds to an integer quantity.
"""
function _is_integer_parameter(cfg::SimulationConfig, key::AbstractString, spec_val::Any)::Bool
    if spec_val isa AbstractVector && length(spec_val) == 2 && all(x -> x isa Integer, spec_val)
        return true
    end
    parts = split(String(key), ".")
    if length(parts) == 2
        sec_sym = Symbol(parts[1])
        fld_sym = Symbol(parts[2])
        if hasfield(SimulationConfig, sec_sym)
            sec_val = getfield(cfg, sec_sym)
            T = typeof(sec_val)
            if hasfield(T, fld_sym)
                fld_type = fieldtype(T, fld_sym)
                return fld_type <: Integer
            end
        end
    end
    return false
end

"""
Constructs Cartesian product combinations of parameter lists.
"""
function _sample_grid_parameters(
    param_keys::Vector{String}, parameters::AbstractDict
)::Vector{Dict{String,Any}}
    sampled_rows = Vector{Dict{String,Any}}()
    value_lists = [parameters[k] for k in param_keys]
    for combo in Iterators.product(value_lists...)
        row = Dict{String,Any}()
        for (k, val) in zip(param_keys, combo)
            row[k] = val
        end
        push!(sampled_rows, row)
    end
    return sampled_rows
end

"""
Generates configurations and run identifiers from sampled parameter combinations.
"""
function _build_ensemble_results(
    spec::EnsembleSweepSpec, sampled_rows::Vector{Dict{String,Any}}
)::Vector{Tuple{String,Dict{String,Any},SimulationConfig}}
    results = Vector{Tuple{String,Dict{String,Any},SimulationConfig}}()
    for (idx, row) in enumerate(sampled_rows)
        run_id = string(spec.sampling_method, "_", lpad(idx, 4, '0'))
        run_out_dir = joinpath(spec.output_dir, run_id)
        overrides = copy(row)
        overrides["output.output_dir"] = run_out_dir
        if !haskey(overrides, "solver.seed")
            overrides["solver.seed"] = spec.seed + idx
        end
        run_cfg = override_config(spec.base_config, overrides)
        push!(results, (run_id, row, run_cfg))
    end
    return results
end

"""
Handles sampling when no parameters are swept.
"""
function _sample_empty_parameters(
    spec::EnsembleSweepSpec,
)::Vector{Tuple{String,Dict{String,Any},SimulationConfig}}
    N = spec.num_samples
    results = Vector{Tuple{String,Dict{String,Any},SimulationConfig}}()
    for i in 1:N
        run_id = string(spec.sampling_method, "_", lpad(i, 4, '0'))
        run_out_dir = joinpath(spec.output_dir, run_id)
        overrides = Dict{String,Any}(
            "output.output_dir" => run_out_dir, "solver.seed" => spec.seed + i
        )
        run_cfg = override_config(spec.base_config, overrides)
        push!(results, (run_id, Dict{String,Any}(), run_cfg))
    end
    return results
end

"""
    sample_parameters(spec)

Sample parameter combinations according to the sweep specification and return
an array of `(run_id, params_dict, config)` tuples.
"""
function sample_parameters(spec::EnsembleSweepSpec)
    param_keys = sort(collect(keys(spec.parameters)))
    D = length(param_keys)
    D == 0 && return _sample_empty_parameters(spec)

    if spec.sampling_method == :grid
        sampled_rows = _sample_grid_parameters(param_keys, spec.parameters)
        return _build_ensemble_results(spec, sampled_rows)
    end

    rng = Random.MersenneTwister(spec.seed)
    sampled_rows = Vector{Dict{String,Any}}()

    if spec.sampling_method == :lhs
        # Latin Hypercube Sampling across D dimensions
        N = spec.num_samples
        sampled_matrix = zeros(Float64, N, D)
        for (j, k) in enumerate(param_keys)
            spec_val = spec.parameters[k]
            lo, hi = if spec_val isa AbstractVector && length(spec_val) == 2
                Float64(spec_val[1]), Float64(spec_val[2])
            else
                throw(
                    ArgumentError(
                        "LHS sampling requires [min, max] range for parameter '$k', got $spec_val",
                    ),
                )
            end
            if _is_integer_parameter(spec.base_config, k, spec_val)
                lo_int = Int(round(lo))
                hi_int = Int(round(hi))
                span = Float64(hi_int - lo_int + 1)
                bin_vals = [lo_int - 0.5 + (i - 1 + rand(rng)) * span / N for i in 1:N]
                sampled_matrix[:, j] .= Random.shuffle(rng, bin_vals)
            else
                bin_vals = [lo + (i - 1 + rand(rng)) * (hi - lo) / N for i in 1:N]
                sampled_matrix[:, j] .= Random.shuffle(rng, bin_vals)
            end
        end
        for i in 1:N
            row = Dict{String,Any}()
            for (j, k) in enumerate(param_keys)
                spec_val = spec.parameters[k]
                if _is_integer_parameter(spec.base_config, k, spec_val)
                    lo_int = Int(round(spec_val[1]))
                    hi_int = Int(round(spec_val[2]))
                    row[k] = clamp(round(Int, sampled_matrix[i, j]), lo_int, hi_int)
                else
                    row[k] = sampled_matrix[i, j]
                end
            end
            push!(sampled_rows, row)
        end
    elseif spec.sampling_method == :random
        # Independent uniform random sampling
        N = spec.num_samples
        for i in 1:N
            row = Dict{String,Any}()
            for k in param_keys
                spec_val = spec.parameters[k]
                lo, hi = if spec_val isa AbstractVector && length(spec_val) == 2
                    Float64(spec_val[1]), Float64(spec_val[2])
                else
                    throw(
                        ArgumentError(
                            "Random sampling requires [min, max] range for parameter '$k', got $spec_val",
                        ),
                    )
                end
                if _is_integer_parameter(spec.base_config, k, spec_val)
                    lo_int = Int(round(lo))
                    hi_int = Int(round(hi))
                    span = Float64(hi_int - lo_int + 1)
                    val = clamp(round(Int, lo_int - 0.5 + rand(rng) * span), lo_int, hi_int)
                    row[k] = val
                else
                    row[k] = lo + rand(rng) * (hi - lo)
                end
            end
            push!(sampled_rows, row)
        end
    else
        throw(
            ArgumentError(
                "Unknown sampling method '$(spec.sampling_method)'. Supported methods: :grid, :lhs, :random",
            ),
        )
    end

    return _build_ensemble_results(spec, sampled_rows)
end

"""
Escapes string values for standard CSV output complying with RFC 4180.
"""
function _escape_csv(val::Any)::String
    s = string(val)
    if occursin(',', s) || occursin('"', s) || occursin('\n', s) || occursin('\r', s)
        return "\"" * replace(s, "\"" => "\"\"") * "\""
    end
    return s
end

"""
    save_ensemble_catalog(catalog, filepath)

Save ensemble catalog records to a CSV file.
"""
function save_ensemble_catalog(
    catalog::AbstractVector{<:AbstractDict{String,Any}}, filepath::AbstractString
)
    isempty(catalog) && return nothing
    mkpath(dirname(filepath))
    key_set = Set{String}()
    for row in catalog
        for k in keys(row)
            push!(key_set, String(k))
        end
    end
    all_keys = sort(collect(key_set))
    # Ensure run_id and status are first
    priority = ["run_id", "status", "walltime_s"]
    ordered_keys = filter(k -> k in all_keys, priority)
    for k in all_keys
        k in ordered_keys || push!(ordered_keys, k)
    end

    open(filepath, "w") do io
        println(io, join([_escape_csv(k) for k in ordered_keys], ","))
        for row in catalog
            vals = [_escape_csv(get(row, k, "")) for k in ordered_keys]
            println(io, join(vals, ","))
        end
    end
    return nothing
end

"""
Execute one ensemble member simulation and record its status.
"""
function _execute_ensemble_member(
    run_id::String, params::Dict{String,Any}, cfg::SimulationConfig, verbose::Bool
)
    if verbose
        @info "Executing ensemble member: $run_id"
    end
    t_start = time()
    status = :success
    err_msg = ""
    try
        simulation_loop(cfg)
    catch e
        status = :failed
        err_msg = sprint(showerror, e)
        @warn "Ensemble member $run_id failed" exception = (e, catch_backtrace())
    end
    elapsed = time() - t_start

    record = Dict{String,Any}(
        "run_id" => run_id,
        "status" => String(status),
        "walltime_s" => round(elapsed; digits=3),
        "error_msg" => err_msg,
    )
    for (k, v) in params
        record[k] = v
    end
    return record
end

"""
    run_ensemble(spec; max_workers=1, verbose=true)

Execute all simulations specified by `spec`, streaming scalar telemetry and compiling
an ensemble catalog file `catalog.csv` in `spec.output_dir`.

When `max_workers > 1`, worker processes must be initialized beforehand (for example via
`Distributed.addprocs`) and `Erebus` loaded on all workers (`@everywhere using Erebus`).
"""
function run_ensemble(spec::EnsembleSweepSpec; max_workers::Integer=1, verbose::Bool=true)
    max_workers >= 1 || throw(ArgumentError("max_workers must be >= 1, got $max_workers"))
    max_workers <= Distributed.nprocs() || throw(
        ArgumentError(
            "max_workers ($max_workers) exceeds available Distributed processes ($(Distributed.nprocs()))",
        ),
    )

    mkpath(spec.output_dir)
    runs = sample_parameters(spec)
    N_runs = length(runs)
    catalog = Vector{Dict{String,Any}}(undef, N_runs)

    if verbose
        @info "Starting ensemble sweep" method = spec.sampling_method total_runs = N_runs output_dir =
            spec.output_dir max_workers = max_workers
    end

    if max_workers > 1
        workers_list = Distributed.workers()
        available_workers = if length(workers_list) >= max_workers
            workers_list[1:max_workers]
        else
            Distributed.procs()[1:max_workers]
        end
        pool = Distributed.WorkerPool(available_workers)
        results = Distributed.pmap(
            run_entry -> begin
                run_id, params, cfg = run_entry
                return _execute_ensemble_member(run_id, params, cfg, verbose)
            end,
            pool,
            runs;
            batch_size=1,
        )
        for idx in 1:N_runs
            catalog[idx] = results[idx]
        end
    else
        for idx in 1:N_runs
            run_id, params, cfg = runs[idx]
            catalog[idx] = _execute_ensemble_member(run_id, params, cfg, verbose)
        end
    end

    catalog_path = joinpath(spec.output_dir, "catalog.csv")
    save_ensemble_catalog(catalog, catalog_path)

    if verbose
        @info "Ensemble sweep completed" total = N_runs catalog = catalog_path
    end

    return catalog
end
