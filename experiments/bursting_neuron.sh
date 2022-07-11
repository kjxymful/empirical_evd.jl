echo "starting exp bursting_neuron"
julia -t 3 --project=.. ../run_data_fit.jl --experiment bursting_neuron --name evd --path_to_data bursting_neuron_20.jld2 &
julia -t 3 --project=.. ../run_data_fit.jl --experiment bursting_neuron --name evd --path_to_data bursting_neuron_22.jld2 $
julia -t 3 --project=.. ../run_data_fit.jl --experiment bursting_neuron --name evd --path_to_data bursting_neuron_24.jld2 $
wait
julia -t 3 --project=.. ../run_data_fit.jl --experiment bursting_neuron --name evd --path_to_data bursting_neuron_25.jld2 $
julia -t 3 --project=.. ../run_data_fit.jl --experiment bursting_neuron --name evd --path_to_data bursting_neuron_28.jld2 $
julia -t 3 --project=.. ../run_data_fit.jl --experiment bursting_neuron --name evd --path_to_data bursting_neuron_34.jld2 $
wait
echo "done"