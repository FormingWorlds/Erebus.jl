for file in [
    ".github/workflows/CI.yml",
    ".github/workflows/nightly.yml",
    ".github/workflows/Coverage.yml",
]
    if isfile(file)
        content = read(file, String)

        # Replace matrix versions
        content = replace(
            content, "          - '1.10'\\n          - '1.11'" => "          - '1.13'"
        )

        # Replace standalone '1.10'
        content = replace(content, "version: '1.10'" => "version: '1.13'")

        write(file, content)
    end
end
