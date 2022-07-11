echo "starting exp lorenz"
julia -t 3 --project=.. ../run_data_fit.jl --experiment lorenz --name evd --path_to_data lorenz_20.jld2 &
julia -t 3 --project=.. ../run_data_fit.jl --experiment lorenz --name evd --path_to_data lorenz_22.jld2 $
julia -t 3 --project=.. ../run_data_fit.jl --experiment lorenz --name evd --path_to_data lorenz_24.jld2 $
wait
julia -t 3 --project=.. ../run_data_fit.jl --experiment lorenz --name evd --path_to_data lorenz_25.jld2 $
julia -t 3 --project=.. ../run_data_fit.jl --experiment lorenz --name evd --path_to_data lorenz_28.jld2 $
julia -t 3 --project=.. ../run_data_fit.jl --experiment lorenz --name evd --path_to_data lorenz_34.jld2 $
wait
echo "done"