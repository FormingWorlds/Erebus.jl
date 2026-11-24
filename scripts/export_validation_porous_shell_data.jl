#!/usr/bin/env julia
# Write the data of the porous shell validation figures as CSV files into a folder.
# Usage: julia --project=. scripts/export_validation_porous_shell_data.jl <output_dir> <suite_dir> <table_dir>
# The source folders hold the runs of scripts/run_validation_porous_shell_suite.jl and
# scripts/run_validation_porous_shell_table.jl; this script checks and copies them.

isdefined(@__MODULE__, :PorousShellValidation) ||
    include(joinpath(@__DIR__, "..", "benchmarks", "validation_porous_shell.jl"))
isdefined(@__MODULE__, :PorousShellSuiteExport) ||
    include(joinpath(@__DIR__, "..", "tools", "export_validation_porous_shell_suite.jl"))

const PROFILE_HEADER = ["r_m", "tmax_K", "tmax_ref_K", "phi", "phi_ref"]
const TABLE_KEY_COLUMNS = ("t_form_myr", "phi_init")

"""
    read_rows(path)

Read the comma-separated file `path` into a vector of rows, each a vector of field strings.
Blank lines are skipped.

# Arguments
- `path::String`: CSV file.

# Returns
- `Vector{Vector{String}}`: the rows in file order, header included.
"""
function read_rows(path::String)
    return [String.(split(strip(l), ",")) for l in eachline(path) if !isempty(strip(l))]
end

"""
    write_rows(path, header, rows)

Write `header` and `rows` to the CSV file `path`.

# Arguments
- `path::String`: output CSV file.
- `header::Vector{String}`: column names.
- `rows::Vector{Vector{String}}`: field strings, each row as long as `header`.

# Returns
- `String`: `path`.
"""
function write_rows(path::String, header::Vector{String}, rows::Vector{Vector{String}})
    open(path, "w") do io
        println(io, join(header, ","))
        for row in rows
            println(io, join(row, ","))
        end
    end
    return path
end

"""
    check_header(path, found, expected)

Check that the header `found` of the file `path` equals `expected`.

# Arguments
- `path::String`: file the header comes from, for the error message.
- `found::Vector{String}`: header read from the file.
- `expected::Vector{String}`: required column names.

# Returns
- `Nothing`.

# Throws
- `ArgumentError`: if the headers differ.
"""
function check_header(path::String, found::Vector{String}, expected::Vector{String})
    found == expected ||
        throw(ArgumentError("$path has columns $(found), expected $(expected)"))
    return nothing
end

"""
    read_suite(suite_dir)

Read the full-evolution suite runs through the checks of `read_suite_results` of
`tools/export_validation_porous_shell_suite.jl`: every case of `SUITE_CASES` present, no
other tag, the case parameters of each row, the pinned 1-D depth and the energy residual
limit. Each case must also have its radial profile `prof_<tag>.csv` with at least one row.

# Arguments
- `suite_dir::String`: folder with `header.csv`, `result_<tag>.csv` and `prof_<tag>.csv`.

# Returns
- `NamedTuple`: `(; tags, rows, profiles)`, the tags in the order of `SUITE_CASES`, the
  result row of each tag as field strings and the profile rows of each tag without the
  header.

# Throws
- `ArgumentError`: from `read_suite_results`, or if a profile is missing, has other
  columns than `PROFILE_HEADER` or holds no row.
"""
function read_suite(suite_dir::String)
    results = PorousShellSuiteExport.read_suite_results(suite_dir)
    tags = [String(r.tag) for r in results]
    rows = [only(read_rows(joinpath(suite_dir, "result_$(tag).csv"))) for tag in tags]
    profiles = Vector{Vector{String}}[]
    for tag in tags
        path = joinpath(suite_dir, "prof_$(tag).csv")
        isfile(path) || throw(ArgumentError("missing $path"))
        content = read_rows(path)
        check_header(path, isempty(content) ? String[] : content[1], PROFILE_HEADER)
        length(content) >= 2 || throw(ArgumentError("$path holds no profile row"))
        push!(profiles, content[2:end])
    end
    return (; tags, rows, profiles)
end

"""
    read_table(table_dir)

Read the 1-D reference table files `g<n>.csv` of `table_dir` into one table, sorted by
`phi_init` and then `t_form_myr`.

# Arguments
- `table_dir::String`: folder with the files `g1.csv`, `g2.csv`, ... of equal header.

# Returns
- `NamedTuple`: `(; header, rows)`, the column names and the sorted rows.

# Throws
- `ArgumentError`: if the folder holds no table file, the headers differ, a row length
  differs from the header, or a (`t_form_myr`, `phi_init`) pair appears twice.
"""
function read_table(table_dir::String)
    files = sort(filter(f -> occursin(r"^g\d+\.csv$", f), readdir(table_dir)))
    isempty(files) && throw(ArgumentError("$table_dir holds no g<n>.csv table file"))
    header = read_rows(joinpath(table_dir, files[1]))[1]
    rows = Vector{String}[]
    for f in files
        path = joinpath(table_dir, f)
        content = read_rows(path)
        check_header(path, content[1], header)
        for row in content[2:end]
            length(row) == length(header) ||
                throw(ArgumentError("$path has a row of $(length(row)) fields"))
            push!(rows, row)
        end
    end
    keycols = [findfirst(==(c), header) for c in TABLE_KEY_COLUMNS]
    any(isnothing, keycols) &&
        throw(ArgumentError("the table header lacks $(TABLE_KEY_COLUMNS)"))
    key(row) = (parse(Float64, row[keycols[2]]), parse(Float64, row[keycols[1]]))
    keys_seen = key.(rows)
    allunique(keys_seen) ||
        throw(ArgumentError("a (t_form_myr, phi_init) pair appears twice in $table_dir"))
    return (; header, rows=rows[sortperm(keys_seen)])
end

"""
    export_all(outdir; suite_dir, table_dir)

Write the suite summary `suite_summary.csv`, the suite profiles `profile_<tag>.csv` and
the 1-D table `table_1d.csv` into `outdir`. Both source folders are read and checked
before `outdir` is created, so a refused folder writes nothing.

# Arguments
- `outdir::String`: output folder, created if missing.
- `suite_dir::String`: folder of the full-evolution suite runs.
- `table_dir::String`: folder of the 1-D reference table files.

# Returns
- `NamedTuple`: `(; tags, n_table)`, the suite tags written and the number of table rows.

# Throws
- `ArgumentError`: from `read_suite` and `read_table`.
"""
function export_all(outdir::String; suite_dir::String, table_dir::String)
    suite = read_suite(suite_dir)
    table = read_table(table_dir)
    mkpath(outdir)
    for (tag, profile) in zip(suite.tags, suite.profiles)
        write_rows(joinpath(outdir, "profile_$(tag).csv"), PROFILE_HEADER, profile)
    end
    header = collect(String, PorousShellValidation.SUITE_RESULT_HEADER)
    write_rows(joinpath(outdir, "suite_summary.csv"), header, suite.rows)
    write_rows(joinpath(outdir, "table_1d.csv"), table.header, table.rows)
    return (; tags=suite.tags, n_table=length(table.rows))
end

"""
    data_export_main(args)

Run the command line `<output_dir> <suite_dir> <table_dir>`. A refused source folder, for
example a suite folder without every case, prints its error.

# Arguments
- `args::Vector{String}`: command-line arguments.

# Returns
- `Int`: exit code, `0` on success, `1` for a refused source folder and `2` for a wrong
  command line.
"""
function data_export_main(args::Vector{String})::Int
    if length(args) != 3 || any(startswith("-"), args)
        println(
            stderr,
            "Usage: julia --project=. scripts/export_validation_porous_shell_data.jl <output_dir> <suite_dir> <table_dir>",
        )
        return 2
    end
    out = try
        export_all(args[1]; suite_dir=args[2], table_dir=args[3])
    catch err
        err isa ArgumentError || rethrow()
        println(stderr, sprint(showerror, err))
        return 1
    end
    println("suite rows: $(length(out.tags)) ($(join(out.tags, ", ")))")
    println("1-D table rows: $(out.n_table)")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(data_export_main(ARGS))
end
