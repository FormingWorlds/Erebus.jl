#!/usr/bin/env julia
# Write the data of the porous shell validation figures as CSV files into a folder.
# Usage: julia scripts/export_validation_porous_shell_data.jl <output_dir> <suite_dir> <table_dir>
# The source folders hold the runs of scripts/run_validation_porous_shell_suite.jl and
# scripts/run_validation_porous_shell_table.jl; this script copies and checks them.

const SUITE_TAGS = (
    "s65_25_025",
    "s129_25_025",
    "s129_225_025",
    "s129_175_010",
    "s129_25_050",
    "s129_25_025_k",
)
const SUITE_HEADER = [
    "tag",
    "Nx",
    "t_form_myr",
    "phi_init",
    "n_steps",
    "k_air",
    "dx",
    "d_shell",
    "q25",
    "q75",
    "d_ref",
    "err_cells",
    "tmax_rel_l2",
    "phi_rms",
    "T_peak",
    "T_peak_ref",
    "energy_res",
    "wall_s",
]
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
    export_suite(suite_dir, outdir)

Copy the full-evolution suite runs: one summary CSV `suite_summary.csv` with a row per
finished tag of `SUITE_TAGS`, and the radial profile of each such tag as
`profile_<tag>.csv`. A tag counts as finished when both `result_<tag>.csv` and
`prof_<tag>.csv` exist in `suite_dir`.

# Arguments
- `suite_dir::String`: folder with `header.csv`, `result_<tag>.csv` and `prof_<tag>.csv`.
- `outdir::String`: output folder.

# Returns
- `Vector{String}`: the tags written, in the order of `SUITE_TAGS`.

# Throws
- `ArgumentError`: if a header differs from the expected columns, a result file does not
  hold exactly one row of the header length, or its tag differs from its file name.
"""
function export_suite(suite_dir::String, outdir::String)
    header_path = joinpath(suite_dir, "header.csv")
    check_header(header_path, read_rows(header_path)[1], SUITE_HEADER)
    rows = Vector{String}[]
    written = String[]
    for tag in SUITE_TAGS
        result_path = joinpath(suite_dir, "result_$(tag).csv")
        profile_path = joinpath(suite_dir, "prof_$(tag).csv")
        (isfile(result_path) && isfile(profile_path)) || continue
        result = read_rows(result_path)
        (length(result) == 1 && length(result[1]) == length(SUITE_HEADER)) || throw(
            ArgumentError("$result_path must hold 1 row of $(length(SUITE_HEADER)) fields"),
        )
        result[1][1] == tag ||
            throw(ArgumentError("$result_path holds tag $(result[1][1]), not $tag"))
        profile = read_rows(profile_path)
        check_header(profile_path, profile[1], PROFILE_HEADER)
        write_rows(joinpath(outdir, "profile_$(tag).csv"), PROFILE_HEADER, profile[2:end])
        push!(rows, result[1])
        push!(written, tag)
    end
    write_rows(joinpath(outdir, "suite_summary.csv"), SUITE_HEADER, rows)
    return written
end

"""
    export_table(table_dir, outdir)

Merge the 1-D reference table files `g<n>.csv` of `table_dir` into `table_1d.csv`, sorted by
`phi_init` and then `t_form_myr`.

# Arguments
- `table_dir::String`: folder with the files `g1.csv`, `g2.csv`, ... of equal header.
- `outdir::String`: output folder.

# Returns
- `Int`: the number of table rows written.

# Throws
- `ArgumentError`: if the folder holds no table file, the headers differ, a row length
  differs from the header, or a (`t_form_myr`, `phi_init`) pair appears twice.
"""
function export_table(table_dir::String, outdir::String)
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
    write_rows(joinpath(outdir, "table_1d.csv"), header, rows[sortperm(keys_seen)])
    return length(rows)
end

"""
    export_all(outdir; suite_dir, table_dir)

Write the suite summary, the suite profiles and the 1-D table into `outdir`.

# Arguments
- `outdir::String`: output folder, created if missing.
- `suite_dir::String`: folder of the full-evolution suite runs.
- `table_dir::String`: folder of the 1-D reference table files.

# Returns
- `NamedTuple`: `(; tags, n_table)`, the suite tags written and the number of table rows.

# Throws
- `ArgumentError`: from `export_suite` and `export_table`.
"""
function export_all(outdir::String; suite_dir::String, table_dir::String)
    mkpath(outdir)
    tags = export_suite(suite_dir, outdir)
    n_table = export_table(table_dir, outdir)
    return (; tags, n_table)
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 3 || error(
        "Usage: export_validation_porous_shell_data.jl <output_dir> <suite_dir> <table_dir>",
    )
    out = export_all(ARGS[1]; suite_dir=ARGS[2], table_dir=ARGS[3])
    missing_tags = [t for t in SUITE_TAGS if !(t in out.tags)]
    println(
        "suite rows: $(length(out.tags)) of $(length(SUITE_TAGS)) ($(join(out.tags, ", ")))"
    )
    isempty(missing_tags) || println("not finished: $(join(missing_tags, ", "))")
    println("1-D table rows: $(out.n_table)")
end
