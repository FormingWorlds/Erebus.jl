# test/test_tooling.jl
# Self-tests for developer tooling and linter infrastructure.

using Test
using TOML

include(joinpath(@__DIR__, "..", "tools", "check_doc_numbers.jl"))
include(joinpath(@__DIR__, "..", "tools", "check_test_quality.jl"))

@testset "Tooling Infrastructure Verification" begin
    @testset "check_doc_numbers validation and rejection self-tests" begin
        # Valid empty map
        mktempdir() do tmpdir
            empty_map = joinpath(tmpdir, "empty_map.toml")
            write(empty_map, "# empty map\n")
            valid, errors = check_doc_numbers(empty_map)
            @test valid == true
            @test length(errors) == 0

            # Missing file error contract
            missing_map = joinpath(tmpdir, "nonexistent.toml")
            valid_missing, errors_missing = check_doc_numbers(missing_map)
            @test valid_missing == false
            @test length(errors_missing) == 1

            # Planted invalid entry with non-existent number
            planted_map = joinpath(tmpdir, "planted_map.toml")
            # Point to an existing doc page and test file, but with an impossible number
            doc_file = joinpath(tmpdir, "doc.md")
            test_file = joinpath(tmpdir, "test.jl")
            write(doc_file, "This page documents speed = 999999.987654321 m/s.\n")
            write(test_file, "@test speed ≈ 123456.0\n")

            # We create the planted map pointing to these files relative to root dirs
            open(planted_map, "w") do io
                println(io, "[[entry]]")
                println(io, "page = \"validation/core_geochemistry.md\"")
                println(io, "number = \"999999999999999.999999999\"")
                println(io, "test_name = \"test_config.jl\"")
            end

            valid_planted, errors_planted = check_doc_numbers(planted_map)
            @test valid_planted == false
            @test length(errors_planted) == 2
            @test any(occursin("Documentation page", err) for err in errors_planted)
            @test any(occursin("Test file", err) for err in errors_planted)
        end
    end

    @testset "check_test_quality AST linter detection self-tests" begin
        # 1. Float equality detection on literal and typed float operands
        expr_float_literal = Meta.parse("@test x == 1.0")
        violations1 = Violation[]
        check_float_equality(expr_float_literal, "planted.jl", 10, violations1)
        @test length(violations1) == 1
        @test violations1[1].rule === :float_equality

        expr_typed_float = Meta.parse("@test x::Float64 == y")
        violations2 = Violation[]
        check_float_equality(expr_typed_float, "planted.jl", 20, violations2)
        @test length(violations2) == 1
        @test violations2[1].rule === :float_equality

        expr_valid_approx = Meta.parse("@test x ≈ 1.0 rtol=1e-6")
        violations_valid = Violation[]
        check_float_equality(expr_valid_approx, "planted.jl", 30, violations_valid)
        @test length(violations_valid) == 0

        # 2. Weak assert detection on !== nothing and numeric thresholds
        expr_nothing = Meta.parse("@test res !== nothing")
        violations3 = Violation[]
        check_weak_asserts(expr_nothing, "planted.jl", 40, violations3)
        @test length(violations3) == 1
        @test violations3[1].rule === :weak_assert

        expr_positivity_zero = Meta.parse("@test temperature > 0.0")
        violations4 = Violation[]
        check_weak_asserts(expr_positivity_zero, "planted.jl", 50, violations4)
        @test length(violations4) == 1
        @test violations4[1].rule === :weak_assert

        expr_threshold = Meta.parse("@test flux > 1e-5")
        violations5 = Violation[]
        check_weak_asserts(expr_threshold, "planted.jl", 60, violations5)
        @test length(violations5) == 1
        @test violations5[1].rule === :weak_assert

        # Valid concrete type assert produces zero violations
        expr_isa = Meta.parse("@test res isa NamedTuple")
        violations_isa = Violation[]
        check_weak_asserts(expr_isa, "planted.jl", 70, violations_isa)
        @test length(violations_isa) == 0

        # 3. Testset assertion count verification
        expr_single_assert_testset = Meta.parse("@testset \"Single\" begin @test 1 == 1 end")
        violations_testset = Violation[]
        check_testsets(expr_single_assert_testset, "planted.jl", 80, violations_testset)
        @test length(violations_testset) == 1
        @test violations_testset[1].rule === :min_asserts
    end
end
