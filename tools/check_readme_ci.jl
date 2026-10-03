#!/usr/bin/env julia
# tools/check_readme_ci.jl
# Verifies that README.md accurately describes the CI workflow defined in .github/workflows/CI.yml,
# and verifies that all README capability, architecture, and restart claims map to real test files
# or CI jobs as defined in tools/readme_claim_map.toml.

using TOML

const ROOT_DIR = normpath(joinpath(@__DIR__, ".."))
const CI_YML_PATH = joinpath(ROOT_DIR, ".github", "workflows", "CI.yml")
const README_PATH = joinpath(ROOT_DIR, "README.md")
const CLAIM_MAP_PATH = joinpath(@__DIR__, "readme_claim_map.toml")

function parse_ci_workflow(ci_path::String)
    isfile(ci_path) || error("CI workflow file not found: $ci_path")
    content = read(ci_path, String)
    lines = split(content, "\n")

    jobs = String[]
    in_jobs = false
    in_version = false
    test_versions = String[]
    test_groups = String[]

    for (i, line) in enumerate(lines)
        stripped = rstrip(line)
        if startswith(stripped, "jobs:")
            in_jobs = true
            continue
        end

        if in_jobs
            # Top-level job names under jobs: are indented by 2 spaces
            m = match(r"^  ([a-zA-Z0-9_\-]+):\s*$", stripped)
            if m !== nothing
                push!(jobs, m.captures[1])
                in_version = false
            end

            if occursin(r"version:\s*$", stripped)
                in_version = true
                continue
            elseif in_version &&
                startswith(stripped, "      ") &&
                !startswith(stripped, "        ") &&
                !startswith(stripped, "      -")
                in_version = false
            end

            # Check test versions under version:
            if in_version
                m_ver = match(r"-\s*'([0-9]+\.[0-9]+)'", stripped)
                if m_ver !== nothing && !(m_ver.captures[1] in test_versions)
                    push!(test_versions, m_ver.captures[1])
                end
            end

            # Check test group
            m_grp = match(r"EREBUS_TEST_GROUP:\s*([a-zA-Z0-9_\-]+)", stripped)
            if m_grp !== nothing && !(m_grp.captures[1] in test_groups)
                push!(test_groups, m_grp.captures[1])
            end
        end
    end

    return (jobs=jobs, versions=test_versions, groups=test_groups)
end

function generate_readme_ci_sentence(ci_info)
    groups_str = join(ci_info.groups, ", ")
    versions_str = join(ci_info.versions, " and ")

    return "3. Open a pull request against `main`. Continuous Integration verifies the test suite (group `$groups_str`) on Julia $versions_str (`ubuntu-latest`), quick simulation on `macos-latest` (`macos-test`), coverage, documentation (`docs`), code style formatting (`format`), and dedicated checks for architecture ratchet (`check_architecture`), performance budget (`check_budget`), bitwise determinism (`check_determinism`), and test quality standards (`test-quality`)."
end

function check_readme_ci(; rewrite::Bool=false)
    ci_info = parse_ci_workflow(CI_YML_PATH)
    expected_sentence = generate_readme_ci_sentence(ci_info)

    readme_content = read(README_PATH, String)
    lines = split(readme_content, "\n")

    # Locate the CI line in Contributing and Code Style
    ci_line_idx = 0
    for (i, line) in enumerate(lines)
        if occursin(r"^3\.\s+Open a pull request against `main`", line)
            ci_line_idx = i
            break
        end
    end

    ci_line_idx > 0 || error("Could not find pull request item 3 in README.md")

    current_line = lines[ci_line_idx]

    if rewrite
        if current_line != expected_sentence
            println("Updating README.md line $(ci_line_idx) to match CI workflow...")
            lines[ci_line_idx] = expected_sentence
            write(README_PATH, join(lines, "\n"))
            println("README.md updated.")
            current_line = expected_sentence
        else
            println("README.md line $(ci_line_idx) is already up to date.")
        end
    end

    # Verification assertions
    all_ok = true

    # Assert exact sentence match
    if current_line != expected_sentence
        println(
            stderr,
            "FAIL: README.md line $(ci_line_idx) does not match CI workflow definition.",
        )
        println(stderr, "  Expected: ", expected_sentence)
        println(stderr, "  Got:      ", current_line)
        all_ok = false
    end

    # 3. Assert required jobs are mentioned
    required_jobs = ["check_architecture", "check_budget", "check_determinism"]
    for j in required_jobs
        if !(j in ci_info.jobs)
            println(stderr, "FAIL: Required job '$j' missing from CI.yml")
            all_ok = false
        end
        if !occursin(j, current_line) && !rewrite
            println(
                stderr,
                "FAIL: Required job '$j' is not mentioned in README.md line $(ci_line_idx)",
            )
            all_ok = false
        end
    end

    # 4. Verify claim mappings if claim map exists
    if isfile(CLAIM_MAP_PATH)
        claim_map = TOML.parsefile(CLAIM_MAP_PATH)

        # Verify capabilities
        if haskey(claim_map, "capabilities")
            for (claim, test_target) in claim_map["capabilities"]
                if !occursin(claim, readme_content)
                    println(stderr, "FAIL: Capability claim not found in README: '$claim'")
                    all_ok = false
                end
                target_path = joinpath(ROOT_DIR, test_target)
                if !isfile(target_path)
                    println(
                        stderr,
                        "FAIL: Test file mapped to capability claim does not exist: '$test_target'",
                    )
                    all_ok = false
                end
            end
        end

        # Verify architecture
        if haskey(claim_map, "architecture")
            for (mod_name, test_target) in claim_map["architecture"]
                if !occursin(mod_name, readme_content)
                    println(
                        stderr, "FAIL: Architecture module '$mod_name' not found in README"
                    )
                    all_ok = false
                end
                target_path = joinpath(ROOT_DIR, test_target)
                if !isfile(target_path)
                    println(
                        stderr,
                        "FAIL: Test file mapped to architecture module does not exist: '$test_target'",
                    )
                    all_ok = false
                end
            end
        end

        # Verify restart
        if haskey(claim_map, "restart")
            for (claim, test_target) in claim_map["restart"]
                if !occursin(claim, readme_content)
                    println(stderr, "FAIL: Restart claim '$claim' not found in README")
                    all_ok = false
                end
                target_path = joinpath(ROOT_DIR, test_target)
                if !isfile(target_path)
                    println(
                        stderr,
                        "FAIL: Test file mapped to restart claim does not exist: '$test_target'",
                    )
                    all_ok = false
                end
            end
        end
    end

    return all_ok
end

function main()
    rewrite = "--rewrite" in ARGS || "--fix" in ARGS
    ok = check_readme_ci(; rewrite=rewrite)

    if ok
        println("All README CI description and claim checks passed.")
        exit(0)
    else
        println(stderr, "README CI check failed.")
        exit(1)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
