echo "starting exp bursting_neuron"
julia -t 3 --project=.. ../run_data_fit.jl --experiment bursting_neuron --name evd -d bursting_neuron_2.0.jld2 &
julia -t 3 --project=.. ../run_data_fit.jl --experiment bursting_neuron --name evd -d bursting_neuron_3.0.jld2 &
julia -t 3 --project=.. ../run_data_fit.jl --experiment bursting_neuron --name evd -d bursting_neuron_5.0.jld2 &
julia -t 3 --project=.. ../run_data_fit.jl --experiment bursting_neuron --name evd -d bursting_neuron_7.0.jld2 &
julia -t 3 --project=.. ../run_data_fit.jl --experiment bursting_neuron --name evd -d bursting_neuron_9.0.jld2 &
julia -t 3 --project=.. ../run_data_fit.jl --experiment bursting_neuron --name evd -d bursting_neuron_10.2.jld2 &
wait
echo "done"