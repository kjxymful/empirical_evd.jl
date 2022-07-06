module empirical_evd

# using Debugger

include("utils.jl")
export load_model

include("data_gen.jl")
export create_series, gen_path, save_series, gen_bif_pars

include("dataloader.jl")
export load_data

end # module
