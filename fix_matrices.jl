for file in [".github/workflows/CI.yml", ".github/workflows/nightly.yml"]
    if isfile(file)
        content = read(file, String)
        content = replace(
            content, "          - '1.13'" => "          - '1.12'\n          - '1.13'"
        )
        write(file, content)
    end
end

content = read(".github/workflows/Coverage.yml", String)
content = replace(content, "version: '1.13'" => "version: '1.12'")
write(".github/workflows/Coverage.yml", content)
