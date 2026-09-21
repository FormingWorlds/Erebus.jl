content = read("Project.toml", String)
content = replace(content, "julia = \"1\"" => "julia = \"1.13\"")
write("Project.toml", content)
