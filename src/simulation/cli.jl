"""
Parse command line arguments and feed them to the main function.

$(SIGNATURES)

# Details:
    
    - nothing

# Returns

    - parsed_args: parsed command line arguments
"""
function parse_commandline(; args::Vector{String}=ARGS)
    s = ArgParseSettings()
    @add_arg_table! s begin
        "config_or_output"
        help = "path to TOML configuration file (.toml) or output directory"
        default = "output"
        "--output_path", "-o"
        help = "output path for simulation data (overrides config output_dir if provided)"
        default = ""
        "--restart", "-r"
        help = "path to JLD2 checkpoint file to resume from"
        default = ""
        "--show_timer"
        help = "show timing results?"
        action = :store_true
        "--force-restart-config"
        help = "override configuration parameters recorded in restart checkpoint"
        action = :store_true
    end
    return parse_args(args, s)
end

"""
    rebuild_cli_restart_config(cfg::SimulationConfig, restart_path::AbstractString)

Return a new `SimulationConfig` identical to `cfg` except with `output.restart_from` set to `restart_path`.
Preserves all other configuration sections intact.
"""
function rebuild_cli_restart_config(cfg::SimulationConfig, restart_path::AbstractString)
    new_output = OutputConfig(;
        output_dir=cfg.output.output_dir,
        savematstep=cfg.output.savematstep,
        visstep=cfg.output.visstep,
        restart_from=String(restart_path),
        mode=cfg.output.mode,
        telemetrystep=cfg.output.telemetrystep,
        telemetry_file=cfg.output.telemetry_file,
        save_final=cfg.output.save_final,
    )
    fields = Pair{Symbol,Any}[]
    for fn in fieldnames(SimulationConfig)
        if fn === :output
            push!(fields, fn => new_output)
        else
            push!(fields, fn => getfield(cfg, fn))
        end
    end
    return SimulationConfig(; fields...)
end

"""
Runs the simulation with the given parameters.

$(SIGNATURES)

# Details

    - nothing

# Returns

    - state: simulation state at completion
"""
function run_simulation(
    config_or_output::AbstractString=""; force_restart_config::Bool=false
)
    if isempty(config_or_output)
        parsed_args = parse_commandline()
        target = parsed_args["config_or_output"]
        cli_output = parsed_args["output_path"]
        cli_restart = parsed_args["restart"]
        show_timer = parsed_args["show_timer"]
        force_restart = parsed_args["force-restart-config"]
    else
        target = config_or_output
        cli_output = ""
        cli_restart = ""
        show_timer = false
        force_restart = force_restart_config
    end

    cfg = if endswith(target, ".toml")
        isfile(target) || throw(ArgumentError("Configuration file does not exist: $target"))
        load_config(target)
    elseif !isempty(target)
        SimulationConfig(; output=OutputConfig(; output_dir=target))
    else
        default_config()
    end

    if !isempty(cli_restart)
        cfg = rebuild_cli_restart_config(cfg, cli_restart)
    end

    actual_output = isempty(cli_output) ? cfg.output.output_dir : cli_output
    actual_output = endswith(actual_output, "/") ? actual_output : actual_output * "/"
    mkpath(actual_output)

    io = open(actual_output * "Erebus_run.log", "w+")
    logger = SimpleLogger(io)
    old_logger = global_logger(logger)
    try
        @info "=========== Erebus simulation run ==========="
        @info "system information: Apple=$(Sys.isapple()) Linux=$(Sys.islinux()) Win=$(Sys.iswindows())" Sys.cpu_info()
        @info "writing results to $actual_output"
        t1 = now()
        @info "start time = $t1"
        state = simulation_loop(
            cfg;
            output_path=actual_output,
            restart_from=cfg.output.restart_from,
            force_restart_config=force_restart,
        )
        t2 = now()
        @info "end time = $t2"
        @info "total run time = $(Dates.canonicalize(
            Dates.CompoundPeriod(t2-t1)))"
        if show_timer && state isa SimulationState
            show(state.timer)
        end
        return state
    finally
        close(io)
        global_logger(old_logger)
    end
end

"""
    run_simulation(cfg::SimulationConfig; restart_from::AbstractString = "", output_path::AbstractString = "", force_restart_config::Bool = false)

Execute a simulation using a pre-loaded `SimulationConfig` object.
"""
function run_simulation(
    cfg::SimulationConfig;
    restart_from::AbstractString="",
    output_path::AbstractString="",
    force_restart_config::Bool=false,
)
    actual_restart = isempty(restart_from) ? cfg.output.restart_from : restart_from
    actual_output = isempty(output_path) ? cfg.output.output_dir : output_path
    return simulation_loop(
        cfg;
        output_path=actual_output,
        restart_from=actual_restart,
        force_restart_config=force_restart_config,
    )
end
