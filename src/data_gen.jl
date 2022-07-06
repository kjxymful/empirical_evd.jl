using DynamicalSystems
using JLD2: save
using BPTT: Dataset as dataset

function create_series(exp::String, μ;u0=[], Δt=0.01, num_t=100000, start_up=200)
    T_end = num_t*Δt
    if exp == "lorenz"
        if isempty(u0)
            ds = Systems.lorenz(; ρ=μ)
        else
            ds = Systems.lorenz(u0; ρ=μ)
        end
    elseif exp == "bursting_neuron"
        if isempty(u0)
            ds = Systems.hodgkinhuxley(;I = μ)
        else
            ds = Systems.hodgkinhuxley(u0; I = μ)
        end
    else
        throw(ArgumentError("$exp not a valid experiment"))
    end

    ts = trajectory(ds, T_end; Δt=Δt, Ttr=start_up)
    return Matrix(ts)
end

function save_series(ts::AbstractMatrix,save_rep::String, save_file::String; opt_info=[])
    save_dict = Dict{String,Any}("time series"=>ts)
    save_path = save_rep*save_file
    if length(opt_info)>0
        for info in opt_info
            merge!(save_dict, Dict(info[1]=>info[2]))
        end
    end
    mkpath(save_rep)
    save(save_path, save_dict)
    return save_path
end

function gen_path(exp::String, μ; eval=false)::Tuple{String,String}
    save_rep = "data/$exp"
    save_file = "/$(exp)_$μ.jld2"
    if eval
        save_rep = "../"*save_rep
    end
    return save_rep, save_file
end

function gen_bif_pars(exp::String)::Vector
    if exp == "lorenz"
        μs = [22 + i for i in 0:7]
    elseif exp == "bursting_neuron"
        μs = [-5 + i for i in 0:17]
    else
        throw(ArgumentError("no data generation possible for experiment: $exp"))
    end
    return μs
end