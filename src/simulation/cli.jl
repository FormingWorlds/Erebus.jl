"""
Parse command line arguments and feed them to the main function.

$(SIGNATURES)

# Details:
    
    - nothing

# Returns

    - parsed_args: parsed command line arguments
"""
function parse_commandline()
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
        arg_type = Bool
        default = false
    end
    return parse_args(s)
end

"""
Runs the simulation with the given parameters.

$(SIGNATURES)

# Details

    - nothing

# Returns

    - nothing 
"""
function run_simulation(config_or_output::AbstractString="")
    if isempty(config_or_output)
        parsed_args = parse_commandline()
        target = parsed_args["config_or_output"]
        cli_output = parsed_args["output_path"]
        cli_restart = parsed_args["restart"]
        show_timer = parsed_args["show_timer"]
    else
        target = config_or_output
        cli_output = ""
        cli_restart = ""
        show_timer = false
    end

    cfg = if endswith(target, ".toml")
        isfile(target) || throw(ArgumentError("Configuration file does not exist: $target"))
        load_config(target)
    elseif !isempty(target)
        def = default_config()
        SimulationConfig(;
            grid=def.grid,
            geometry=def.geometry,
            time=def.time,
            solver=def.solver,
            poroelasticity=def.poroelasticity,
            thermodynamics=def.thermodynamics,
            materials=def.materials,
            output=OutputConfig(; output_dir=target),
        )
    else
        default_config()
    end

    if !isempty(cli_restart)
        cfg = SimulationConfig(;
            grid=cfg.grid,
            geometry=cfg.geometry,
            time=cfg.time,
            solver=cfg.solver,
            poroelasticity=cfg.poroelasticity,
            thermodynamics=cfg.thermodynamics,
            materials=cfg.materials,
            output=OutputConfig(;
                output_dir=cfg.output.output_dir,
                savematstep=cfg.output.savematstep,
                visstep=cfg.output.visstep,
                restart_from=cli_restart,
            ),
        )
    end

    actual_output = isempty(cli_output) ? cfg.output.output_dir : cli_output
    actual_output = endswith(actual_output, "/") ? actual_output : actual_output * "/"
    mkpath(actual_output)

    io = open(actual_output * "Erebus_run.log", "w+")
    logger = SimpleLogger(io)
    global_logger(logger)
    if show_timer
        reset_timer!(to)
    end
    @info "=========== Erebus simulation run ==========="
    @info "system information: Apple=$(Sys.isapple()) Linux=$(Sys.islinux()) Win=$(Sys.iswindows())" Sys.cpu_info()
    @info "writing results to $actual_output"
    t1 = now()
    @info "start time = $t1"
    simulation_loop(cfg; output_path=actual_output, restart_from=cfg.output.restart_from)
    t2 = now()
    @info "end time = $t2"
    @info "total run time = $(Dates.canonicalize(
        Dates.CompoundPeriod(t2-t1)))"
    if show_timer
        show(to)
    end
    return close(io)
end

"""
    run_simulation(cfg::SimulationConfig; restart_from::AbstractString = "", output_path::AbstractString = "")

Execute a simulation using a pre-loaded `SimulationConfig` object.
"""
function run_simulation(
    cfg::SimulationConfig; restart_from::AbstractString="", output_path::AbstractString=""
)
    actual_restart = isempty(restart_from) ? cfg.output.restart_from : restart_from
    actual_output = isempty(output_path) ? cfg.output.output_dir : output_path
    return simulation_loop(cfg; output_path=actual_output, restart_from=actual_restart)
end
