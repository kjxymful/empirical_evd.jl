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
        args["latent_dim"] = 26
        args["teacher_forcing_interval"] = 5
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

            n_threads = Threads.nthreads()
            println("Running on $n_threads Thread(s)")

            # get computing device
            device = get_device(args)

            # data
            D = load_dataset(path_to_data, "DATA", device=device)

            # model
            plrnn = initialize_model(args, D) |> device

            # optimizer
            opt = initialize_optimizer(args)

            # create directories
            save_path = create_folder_structure(args["experiment"], args["name"], args["run"])

            # store hypers
            store_hypers(args, save_path)

            train_!(plrnn, D, opt, args, save_path)
        end
    else
        args["path_to_data"] = dir_path*"/"*args["path_to_data"]# in my case only give filename
        main_routine(args)
    end
end

main()