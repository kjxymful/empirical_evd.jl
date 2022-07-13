echo "starting exp lorenz"
julia -t 3 --project=.. ../run_data_fit.jl --experiment lorenz --name evd -d lorenz_20.jld2 -e 1500 &
julia -t 3 --project=.. ../run_data_fit.jl --experiment lorenz --name evd -d lorenz_22.jld2 -e 1500 &
julia -t 3 --project=.. ../run_data_fit.jl --experiment lorenz --name evd -d lorenz_24.jld2 -e 1500 &
wait
julia -t 3 --project=.. ../run_data_fit.jl --experiment lorenz --name evd -d lorenz_25.jld2 -e 1500 &
julia -t 3 --project=.. ../run_data_fit.jl --experiment lorenz --name evd -d lorenz_28.jld2 -e 1500 &
julia -t 3 --project=.. ../run_data_fit.jl --experiment lorenz --name evd -d lorenz_34.jld2 -e 1500 &
wait
echo "done"