using TOML

input = TOML.parsefile("input.toml")
print(input["Grid"])
# input["BoundaryCondition"][1]["include"]