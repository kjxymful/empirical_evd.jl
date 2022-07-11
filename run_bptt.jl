#using Revise
using BPTT # bppt-julia at https://gitlab.zi.local/Florian.Hess/bptt-julia.git
ENV["GKSwstype"] = "nul"

main_routine(parse_commandline())
