using BPTT
"""
    main_routine(args)

Function executed by every worker process.
"""
function bptt_routine(args::Dict{String,Any})
    # num threads
    n_threads = Threads.nthreads()
    println("Running on $n_threads Thread(s)")

    # get computing device
    device = get_device(args)

    # data
    D, par = load_dataset(args["path_to_data"], "DATA", device=device)

    # model
    plrnn = initialize_model(args, D) |> device

    # optimizer
    opt = initialize_optimizer(args)

    # create directories
    save_path = create_folder_structure(args["experiment"]*"_$par", args["name"], args["run"])

    # store hypers
    store_hypers(args, save_path)

    train_!(plrnn, D, opt, args, save_path)
end
