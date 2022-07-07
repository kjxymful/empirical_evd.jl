module empirical_evd

# using Debugger

include("utils.jl")
export load_model, load_data

include("data_gen.jl")
export create_series, gen_path, save_series, gen_bif_pars

include("ds_models.jl")
export bursting_neuron

end # module
