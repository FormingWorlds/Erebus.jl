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

            # Valid passing map using custom docs_dir and test_dir
            mock_doc = joinpath(tmpdir, "mock_page.md")
            mock_test = joinpath(tmpdir, "mock_test.jl")
            write(mock_doc, "Documented speed = 299792458 m/s.\n")
            write(mock_test, "@test speed == 299792458\n")
            valid_map = joinpath(tmpdir, "valid_map.toml")
            open(valid_map, "w") do io
                println(io, "[[entry]]")
                println(io, "page = \"mock_page.md\"")
                println(io, "number = \"299792458\"")
                println(io, "test_name = \"mock_test.jl\"")
            end
            valid_res, errors_res = check_doc_numbers(
                valid_map; docs_dir=tmpdir, test_dir=tmpdir
            )
            @test valid_res == true
            @test isempty(errors_res)

            # Planted invalid entry with non-existent number in existing files
            planted_map_num = joinpath(tmpdir, "planted_num_map.toml")
            open(planted_map_num, "w") do io
                println(io, "[[entry]]")
                println(io, "page = \"mock_page.md\"")
                println(io, "number = \"999999999999\"")
                println(io, "test_name = \"mock_test.jl\"")
            end
            valid_p_num, errors_p_num = check_doc_numbers(
                planted_map_num; docs_dir=tmpdir, test_dir=tmpdir
            )
            @test valid_p_num == false
            @test length(errors_p_num) == 2
            @test any(occursin("does not contain number", err) for err in errors_p_num)

            # Planted invalid entry with missing documentation and test files
            missing_file_map = joinpath(tmpdir, "missing_file_map.toml")
            open(missing_file_map, "w") do io
                println(io, "[[entry]]")
                println(io, "page = \"nonexistent_doc.md\"")
                println(io, "number = \"299792458\"")
                println(io, "test_name = \"nonexistent_test.jl\"")
            end
            valid_mf, errors_mf = check_doc_numbers(
                missing_file_map; docs_dir=tmpdir, test_dir=tmpdir
            )
            @test valid_mf == false
            @test length(errors_mf) == 2
            @test any(occursin("Documentation page not found", err) for err in errors_mf)
            @test any(occursin("Test file not found", err) for err in errors_mf)
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

        # Valid error bound and memory allocation checks do NOT produce weak_assert violations
        expr_err_bound = Meta.parse("@test rel_diff < 1e-12")
        violations_err = Violation[]
        check_weak_asserts(expr_err_bound, "planted.jl", 65, violations_err)
        @test length(violations_err) == 0

        expr_alloc_bound = Meta.parse("@test alloc_ldiv < 1024")
        violations_alloc = Violation[]
        check_weak_asserts(expr_alloc_bound, "planted.jl", 66, violations_alloc)
        @test length(violations_alloc) == 0

        # Negativity / zero upper bound assertions
        expr_neg_zero = Meta.parse("@test x < 0")
        violations_neg = Violation[]
        check_weak_asserts(expr_neg_zero, "planted.jl", 67, violations_neg)
        @test length(violations_neg) == 1
        @test violations_neg[1].rule === :weak_assert

        expr_neg_zero_f = Meta.parse("@test x <= 0.0")
        violations_neg_f = Violation[]
        check_weak_asserts(expr_neg_zero_f, "planted.jl", 68, violations_neg_f)
        @test length(violations_neg_f) == 1
        @test violations_neg_f[1].rule === :weak_assert

        # Reversed operand checks
        expr_rev_pos = Meta.parse("@test 0 < speed")
        violations_rev_pos = Violation[]
        check_weak_asserts(expr_rev_pos, "planted.jl", 69, violations_rev_pos)
        @test length(violations_rev_pos) == 1
        @test violations_rev_pos[1].rule === :weak_assert

        expr_rev_type = Meta.parse("@test Float64 === typeof(val)")
        violations_rev_type = Violation[]
        check_weak_asserts(expr_rev_type, "planted.jl", 70, violations_rev_type)
        @test length(violations_rev_type) == 1
        @test violations_rev_type[1].rule === :weak_assert

        # Macro kwargs support (e.g. broken=true)
        expr_broken_float = Meta.parse("@test broken=true x == 1.0")
        violations_broken = Violation[]
        check_float_equality(expr_broken_float, "planted.jl", 71, violations_broken)
        @test length(violations_broken) == 1
        @test violations_broken[1].rule === :float_equality

        # @test_broken macro support
        expr_tb_float = Meta.parse("@test_broken x == 1.0")
        violations_tb = Violation[]
        check_float_equality(expr_tb_float, "planted.jl", 72, violations_tb)
        @test length(violations_tb) == 1
        @test violations_tb[1].rule === :float_equality

        # Valid concrete type assert produces zero violations
        expr_isa = Meta.parse("@test res isa NamedTuple")
        violations_isa = Violation[]
        check_weak_asserts(expr_isa, "planted.jl", 70, violations_isa)
        @test length(violations_isa) == 0

        # 3. Testset assertion count verification
        expr_single_assert_testset = Meta.parse(
            "@testset \"Single\" begin @test 1 == 1 end"
        )
        violations_testset = Violation[]
        check_testsets(expr_single_assert_testset, "planted.jl", 80, violations_testset)
        @test length(violations_testset) == 1
        @test violations_testset[1].rule === :min_asserts

        # Parameterized testset with single assertion triggers :min_asserts
        expr_for_testset = Meta.parse("@testset \"Loop\" for i in 1:3 @test i == i end")
        violations_for_testset = Violation[]
        check_testsets(expr_for_testset, "planted.jl", 85, violations_for_testset)
        @test length(violations_for_testset) == 1
        @test violations_for_testset[1].rule === :min_asserts
    end

    @testset "GoldenHelpers bitwise comparator self-tests" begin
        include(joinpath(@__DIR__, "golden_helpers.jl"))
        using .GoldenHelpers: compare_golden

        # Identical NamedTuples match
        nt1 = (; a=1.0, b=[2.0, NaN], c="test")
        nt2 = (; a=1.0, b=[2.0, NaN], c="test")
        @test compare_golden(nt1, nt2) == true

        # NamedTuple vs Dict with same string/symbol keys matches
        dict_rep = Dict("a" => 1.0, "b" => [2.0, NaN], "c" => "test")
        @test compare_golden(nt1, dict_rep) == true
        @test compare_golden(dict_rep, nt1) == true

        # Key set mismatch throws
        dict_mismatch = Dict("a" => 1.0, "b" => [2.0, NaN])
        @test_throws ErrorException compare_golden(nt1, dict_mismatch)

        # Array size mismatch throws
        nt_size_mismatch = (; a=1.0, b=[2.0], c="test")
        @test_throws ErrorException compare_golden(nt1, nt_size_mismatch)

        # Array eltype mismatch throws
        nt_eltype_mismatch = (; a=1.0, b=Float32[2.0, NaN], c="test")
        @test_throws ErrorException compare_golden(nt1, nt_eltype_mismatch)

        # Signed zero (+0.0 vs -0.0) throws
        nt_pos_zero = (; z=0.0)
        nt_neg_zero = (; z=-0.0)
        @test_throws ErrorException compare_golden(nt_pos_zero, nt_neg_zero)
    end
end
