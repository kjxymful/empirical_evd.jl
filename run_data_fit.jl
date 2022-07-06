include("src/empirical_evd.jl")
using .empirical_evd
using BPTT # bppt-julia at https://gitlab.zi.local/Florian.Hess/bptt-julia.git
using JLD2
using Flux
ENV["GKSwstype"] = "nul"

# using Debugger

const gen_data = true

function main()
    args = parse_commandline()
    # parse args
    println("Parsed args:")
    for (arg, val) in args
        println("  $arg  --->  $val")
    end
    # experiment
    experiment = args["experiment"]
    println("EXP $experiment")
    # run number
    run = args["run"]
    println("RUN $run")

    # num threads
    n_threads = Threads.nthreads()
    println("Running on $n_threads Thread(s)")

    # data
    if gen_data
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
    data_files = readdir(dir_path)
    
    for exp_file in data_files
        println("Fitting PlRNN on $exp_file")
        path_to_data = dir_path *"/"*exp_file
        D, par = load_data(path_to_data, "DATA", device=cpu)

        # model
        plrnn = initialize_model(args, D)

        # optimizer
        opt = initialize_optimizer(args)

        # create directories for saving ProgressMeter    # saving stuff here
        save_path = create_folder_structure(experiment*"_$par", run)
        train_!(plrnn, D, opt, args, save_path)
    end
end

# main()
