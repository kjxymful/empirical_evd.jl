include("src/empirical_evd.jl")
using .empirical_evd
using BPTT # bppt-julia at https://gitlab.zi.local/Florian.Hess/bptt-julia.git
using JLD2
using Flux
using Base.Threads: @threads
ENV["GKSwstype"] = "nul"
# using Debugge
BPTT.load_dataset(path::String, name::String; device=cpu) = load_data(path, name;device)

function main()
    # parse args
    args = parse_commandline()

    # experiment
    experiment = args["experiment"]
    println("EXP $experiment")

    if experiment == "lorenz"
        args["latent_dim"] = 15
        args["teacher_forcing_interval"] = 10
    elseif experiment == "bursting_neuron"
        if args["latent_dim"] ==10
            args["latent_dim"] = 26
        end
        args["teacher_forcing_interval"] = 10
    end

    # data
    if !isdir("data/$experiment")
        println("Generating Data")
        μs = gen_bif_pars(experiment)
        for μ in μs
            tseries = create_series(experiment, μ)
            save_series(tseries,
                gen_path(experiment, μ)...
                ;
                opt_info=[("par", μ)])
        end
    end

    dir_path, _ = gen_path(experiment, 0)

    if args["path_to_data"] == "example_data/lorenz.npy"
        data_files = readdir(dir_path)

        for exp_file in data_files
            println("Fitting PlRNN, $exp_file on Thread $(Threads.threadid())")
            path_to_data = dir_path * "/" * exp_file
            args["path_to_data"] = path_to_data

            bptt_routine(args)
        end
    else
        println("Fitting PlRNN, $(args["path_to_data"]) on Thread $(Threads.threadid())")
        args["path_to_data"] = dir_path*"/"*args["path_to_data"]# in my case only give filename
        
        bptt_routine(args)
    end
end

main()