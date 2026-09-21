for file in [
    ".github/workflows/CI.yml",
    ".github/workflows/nightly.yml",
    ".github/workflows/Coverage.yml",
]
    if isfile(file)
        content = read(file, String)

        # We can just use regex to replace the versions list
        content = replace(
            content,
            r"version:\s*\n\s+- '1\.10'\s*\n\s+- '1\.11'" => "version:\n          - '1.13'",
        )

        write(file, content)
    end
end
