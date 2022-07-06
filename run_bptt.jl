#using Revise
using BPTT # bppt-julia at https://gitlab.zi.local/Florian.Hess/bptt-julia.git
using NPZ
using Flux
ENV["GKSwstype"] = "nul"

function main(args)
    # parse args
    println("Parsed args:")
    for (arg, val) in args
        println("  $arg  --->  $val")
    end

    # run number
    run = args["run"]
    println("RUN $run")

    # num threads
    n_threads = Threads.nthreads()
    println("Running on $n_threads Thread(s)")

    # data
    D = load_dataset(args["path_to_data"], "DATA", device=cpu)

    # model
    plrnn = initialize_model(args, D)

    # optimizer
    opt = initialize_optimizer(args)

    # create directories for saving ProgressMeter    # saving stuff here
    save_path = create_folder_structure(args["experiment"], run)

    train_!(plrnn, D, opt, args, save_path)
end

main(parse_commandline())
