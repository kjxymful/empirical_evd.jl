include("../src/empirical_evd.jl")
using .empirical_evd
using Plots
using BPTT
using LinearAlgebra
using Measures
using SplitApplyCombine
using StatsBase

function bursting_all()
    bins = 20
    system = "bursting_neuron"
    dir = "../Figures/$(system)_all"
    function remove!(a, item)
        deleteat!(a, findall(x -> x == item, a))
    end

    # for making subplots for each experiment
    function make_subplots(exp, plots, type)
        if exp == "lorenz"
            μ = "ρ"
            title = ["\n20" "ρ as bifurcation parameter\n22" "\n24" "\n25" "\n28" "\n34"]
        else
            μ = "\$g_{nmda}\$"
            title = ["\n2" "$μ as bifurcation parameter\n3" "\n5" "\n9" "\n10" "\n10.2"]
        end
        p = plot(plots..., layout=length(plots),
            title=title, titlefont=font(11),
            plot_title="$type of evd with \n $exp plrnn", plot_titlevspan=0.15
        )
        type_str = replace(type, " " => "_") # replace the spaces for savefig
        savefig(p, "../Figures/$(exp)_all/$(exp)_$(type_str)_evd.png")
    end


    experiments = readdir("../Results")
    remove!(experiments, "default")
    deleteat!(experiments, findall(x -> !occursin(system, x), experiments))

    # safe the latent dim of the plrnn(from the first PLRNN)
    model_path = load_result_path(experiments[1], 5000, 1)
    m = BPTT.load_model(model_path)
    dim = length(m.A)
    n_quadrants = 2^dim
    n_iter = 10000
    n_parts =n_quadrants÷n_iter

    # loop over all possible matrices D
    all_Ds = ([0, 1] for i in 1:dim)
    all_Ds = Iterators.product(all_Ds...)
    @assert length(all_Ds) == n_quadrants

    # split the iterator in 100 parts for memory reasons
    partition_Ds = Iterators.partition(all_Ds, n_parts)

    # loop over part of the possible values
for (d_part, Ds) in enumerate(partition_Ds)
    if d_part == 10
        break
    end
    @show d_part
    global exp_λs = nothing
    global exp_λs = Vector{Vector{Vector{ComplexF32}}}()
    # load all the experiments
    for exp in experiments
        model_path = load_result_path(exp, 5000, 1)
        m = BPTT.load_model(model_path)

        # copy all images the first time this is run
        if d_part == 1
            mkpath(dir)
            dst_path = "$dir/$(exp)_traj.png"
            # cp(load_result_path(exp, 5000, 1; img=true), dst_path, force=true)
        end

        global all_λs = Vector{Vector{ComplexF32}}()

        for D in Ds
            Dt = diagm(collect(D))
            # compose single matrix A+WD
            AWD = diagm(m.A) + m.W * Dt

            # get evd of AWD
            EVD = eigen(AWD)
            λs = EVD.values
            # add the same for Vs by doing the same with labda->v
            push!(all_λs, complex(λs))
        end
        push!(exp_λs, all_λs)
    end
    all_λs = nothing
    GC.gc()

    # convert to array
    arr_l = combinedims(exp_λs)
    exp_λs = nothing
    arr_l = combinedims(arr_l)
    arr_l = permutedims(arr_l, (3, 2, 1))
    # sort all the arrays
    res_l = zeros(size(arr_l))
    ims_l = zeros(size(arr_l))
    abs_l = zeros(size(arr_l))

    for ex in 1:length(experiments)
        v = nothing
        for i in 1:n_parts
            v = sortperm(abs.(arr_l[ex, i, :]))
            res_l[ex, i, :] = real(arr_l[ex, i, v])
            ims_l[ex, i, :] = imag(arr_l[ex, i, v])
            abs_l[ex, i, :] = abs.(arr_l[ex, i, v])
        end
    end
    arr_l = nothing

    if d_part == 1
        # histograms of single experiments of single experiments
        for exp_id in 1:length(experiments)
            h_abs = fit(Histogram, collect(Iterators.flatten(abs_l)), nbins=bins)
            h_ims = fit(Histogram, collect(Iterators.flatten(ims_l)), nbins=bins)
            h_res = fit(Histogram, collect(Iterators.flatten(res_l)), nbins=bins)

            if exp_id == 1
                global edges_abs = [h_abs.edges[1]]
                global edges_ims = [h_ims.edges[1]]
                global edges_res = [h_res.edges[1]]
                global counts_abs = [h_abs.weights]
                global counts_ims = [h_ims.weights]
                global counts_res = [h_res.weights]
            else
                push!(edges_abs, h_abs.edges[1])
                push!(edges_ims, h_ims.edges[1])
                push!(edges_res, h_res.edges[1])

                push!(counts_abs, h_abs.weights)
                push!(counts_ims, h_ims.weights)
                push!(counts_res, h_res.weights)
            end
        end

        # plot one detailed plot of the distribution the first time this is run
        _, index = findmax(abs_l[1, :, :])
        # abs values together for one case
        c = [:orange, :red, :darkred, :purple, :black]
        label = ["cycle", "1 loop", "4 loops", "almost bursting", "bursting"]
        for i in 1:length(experiments)
            try
                plot!(abs_plot, abs_l[i, index[1], :], 1:dim, c=c[i], xlabel="Ev value [b.E.]",
                    ylabel="Ev No", label=label[i])
            catch
                global abs_plot = plot(abs_l[i, index[1], :], 1:dim,
                    c=c[i],
                    label=label[i],
                    xlabel="Ev value [b.E.]",
                    ylabel="Ev No",
                    plot_title="Eigenvalue evolution of quadrant $(index[1])", plot_titlevspan=0.05)
            end
        end
        savefig(abs_plot, "../Figures/$(system)_all/change_in_abs_vals.png")
    else
        # histograms of single experiments of single experiments
        for exp_id in 1:length(experiments)
            h_abs = fit(Histogram, collect(Iterators.flatten(abs_l)), edges_abs[exp_id])
            h_ims = fit(Histogram, collect(Iterators.flatten(ims_l)), edges_ims[exp_id])
            h_res = fit(Histogram, collect(Iterators.flatten(res_l)), edges_res[exp_id])
            counts_abs[exp_id] += h_abs.weights
            counts_ims[exp_id] += h_ims.weights
            counts_res[exp_id] += h_res.weights
        end
    end

    change_abs = Vector{Matrix{Float64}}()
    # look at change of the parameters
    for id in 2:length(experiments) # for all experiment
        push!(change_abs, abs.(abs_l[id, :, :] - abs_l[id-1, :, :]))
    end

    change_abs_p = Vector{Matrix{Float64}}()
    # percentual change
    for id in 2:length(experiments) # for all experiments
        push!(change_abs_p, abs.(abs_l[id, :, :] ./ abs_l[id-1, :, :]))
    end

    arr_change_abs = combinedims(change_abs, 1)
    arr_change_abs_p = combinedims(change_abs_p, 1)

    change_abs = nothing
    change_abs_p = nothing
    idcs = findall(x -> x >= 2, arr_change_abs_p)
    for id in idcs
        arr_change_abs_p[id] = 0
    end

    if d_part == 1
        # plot the change of some of the imoprtant eigenvalues
        _, idx = findmax(sum(arr_change_abs, dims=1))
        _, idx_p = findmax(sum(arr_change_abs_p, dims=1))

        # specific plot, to show change of the most changing evs, also just the firs time this is run
        a = plot(abs_l[:, idx[2], idx[3]], xlabel="experiment", ylabel="b.E.")
        b = plot(abs_l[:, idx_p[2], idx_p[3]], xlabel="experiment", ylabel="%")

        d = plot(a, b, layout=2,
            size=(1000, 500),
            legend=nothing,
            margin=10mm,
            title=["absolute change" "percentual change"],
            plot_title="detailed changes of absolute values", plot_titlevspan=0.05)
        savefig(d, "../Figures/bursting_neuron_all/ev_t_detail.png")
    end

    if d_part == 1
        h_abs_all = fit(Histogram, collect(Iterators.flatten(arr_change_abs)), nbins=bins)
        h_abs_p = fit(Histogram, collect(Iterators.flatten(arr_change_abs_p)), nbins=bins)

        global counts_abs_p = h_abs_p.weights
        global counts_abs_all = h_abs_all.weights
        global edges_abs_p = h_abs_p.edges[1]
        global edges_abs_all = h_abs_all.edges[1]
    else
        h_abs_all = fit(Histogram, collect(Iterators.flatten(arr_change_abs)), edges_abs_all)
        h_abs_p = fit(Histogram, collect(Iterators.flatten(arr_change_abs_p)), edges_abs_p)

        counts_abs_p += h_abs_p.weights
        counts_abs_all += h_abs_all.weights
    end
    arr_change_abs = nothing
    arr_change_abs_p = nothing
end

    HIST_abs_p = Histogram(edges_abs_p, counts_abs_p)
    HIST_abs_all = Histogram(edges_abs_all, counts_abs_all)

    HISTS_abs = [Histogram(edges_abs[i], counts_abs[i]) for i in 1:length(experiments)]
    HISTS_ims = [Histogram(edges_ims[i], counts_ims[i]) for i in 1:length(experiments)]
    HISTS_res = [Histogram(edges_res[i], counts_res[i]) for i in 1:length(experiments)]

    # histogram of all abs values
    hists_abs_val = (plot(HISTS_abs[i], yaxis=:log, legend=nothing, xlabel="Ev value [b.E.]", ylabel="Hits") for i in 1:length(experiments))
    # im and real part histograms
    hists_ims = (plot(HISTS_ims[i], yaxis=:log, xlabel="Ev value [b.E.]") for i in 1:length(experiments))
    hists_res = (plot(HISTS_res[i], yaxis=:log, xlabel="Ev value [b.E.]") for i in 1:length(experiments))

    make_subplots("bursting_neuron", hists_abs_val, "histogram of absolute values")
    make_subplots("bursting_neuron", hists_ims, "histogram of imaginary part")
    make_subplots("bursting_neuron", hists_res, "historam of real parts")

    # hist of overall change in all eigenvalues
    hist_abs_p = plot(HIST_abs_p,
        title="percentual change of all eigenvalues",
        yaxis=:log,
        legend=nothing)

    # hist of overall change in all eigenvalues absou change
    hist_abs = plot(HIST_abs_all,
        yaxis=:log,
        title="absolute change of all eigenvalues",
        legend=nothing)

    savefig(hist_abs_p, "../Figures/bursting_neuron_all/overall_change_abs_val_p.png")
    savefig(hist_abs, "../Figures/bursting_neuron_all/overall_change_abs_val.png")
end
