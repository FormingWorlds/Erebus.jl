content = read("tools/check_test_quality.jl", String)
content = replace(
    content,
    "return assert_count, has_sub_testsets\\nend\\n\\nfunction check_testsets" => "return assert_count, has_sub_testsets, throws_count\\nend\\n\\nfunction check_testsets",
)
write("tools/check_test_quality.jl", content)
