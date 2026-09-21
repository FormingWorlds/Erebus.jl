content = read("tools/check_test_quality.jl", String)
old = """
                    )
                end
                # Check for bare positivity: @test x > 0, @test x >= 0, @test x < 0, etc.
                if (op === :(>) || op === :(>=) || op === :(<) || op === :(<=)) &&
"""
new = """
                    )
                # Check for bare positivity: @test x > 0, @test x >= 0, @test x < 0, etc.
                elseif (op === :(>) || op === :(>=) || op === :(<) || op === :(<=)) &&
"""
content = replace(content, old => new)
write("tools/check_test_quality.jl", content)
