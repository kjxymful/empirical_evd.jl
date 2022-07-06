using JLD2: load
using BPTT

function load_data(path::String, name::String; device = cpu)
    data_dict = load(path)
    X = Float32.(data_dict["time series"]) |> device
    if "par" in keys(data_dict)
        par = data_dict["par"]
    else
        par = nothing
    end

    return BPTT.Dataset(X, name), par
end