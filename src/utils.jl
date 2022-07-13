using BSON: @load
import JLD2: load as jld2load
using BPTT

function load_model(exp::String, name::String, run::Int, epoch::Int; eval=true)::BPTT.AbstractPLRNN
    load_path = joinpath(["Results", exp, name, format_run_ID(run)])
    if eval
        load_path = "../" * load_path
    end
    @load joinpath(load_path, "checkpoints", "model_$epoch.bson") m
    return m
end

function format_run_ID(run::Int)::String
    run_str = string(run)
    n_digits = length(run_str)
    @assert n_digits < 4
    add_zeros = 3 - n_digits
    return repeat("0", add_zeros) * run_str
end

function load_data(path::String, name::String; device=cpu)
    data_dict = jld2load(path)
    X = Float32.(data_dict["time series"]) |> device
    if "par" in keys(data_dict)
        par = data_dict["par"]
    else
        par = nothing
    end

    return BPTT.Dataset(X, name), par
end

function load_result_path(exp::String, epoch, run; eval=true)
    run = format_run_ID(run)
    path = "Results/" * exp * "/evd/$run/checkpoints/" * "model_$epoch.bson"
    if eval
        path = "../" * path
    end
    return path
end