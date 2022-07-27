include("../src/empirical_evd.jl")
using .empirical_evd
using Plots
using BPTT
using LinearAlgebra
using Measures
using SplitApplyCombine
using StatsBase

function main()
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

    # loop over all possible matrices D
    all_Ds = ([0, 1] for i in 1:dim)
    all_Ds = Iterators.product(all_Ds...)
    @assert length(all_Ds) == n_quadrants

    # split the iterator in 100 parts for memory reasons
    partition_Ds = Iterators.partition(all_Ds, n_quadrants÷100)

    # loop over part of the possible values
    for (d_part, Ds) in enumerate(partition_Ds)
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
    
            global all_λs = nothing
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
        # convert to array
        abs_l = combinedims(exp_λs)
        abs_l = combinedims(abs_l)
        exp_λs = splitdims(abs_l)
        abs_l = permutedims(abs_l, (3, 2, 1))
    
        # sort all the arrays
        res_l = zeros(6, 2^dim, dim)
        ims_l = zeros(6, 2^dim, dim)
    
        for (ex, all_lambdas) in enumerate(exp_λs)
            v = nothing
            for i in 1:size(all_lambdas)[2]
                v = sortperm(abs_ls[:, i])
                abs_l[ex, i, :] = abs.(all_lambdas)[v, i]
                res_l[ex, i, :] = real(all_lambdas)[v, i]
                ims_l[ex, i, :] = imag(all_lambdas)[v, i]
            end
        end
        exp_λs = nothing
    
        if d_part == 1
            # histograms of single experiments of single experiments
            for exp_id in 1:length(experiments)
                h_abs = fit(Histogram, abs_l, nbins=bins)
                h_ims = fit(Histogram, ims_l, nbins=bins)
                h_res = fit(Histogram, res_l, nbins=bins)
    
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
            savefig(abs_plot, "../Figures/$(exp_str)_all/change_in_abs_vals.png")
        else
            # histograms of single experiments of single experiments
            for exp_id in 1:length(experiments)
                h_abs = fit(Histogram, abs_l, nbins=bins)
                h_ims = fit(Histogram, ims_l, nbins=bins)
                h_res = fit(Histogram, res_l, nbins=bins)
                counts_abs[exp_id] += h_abs.weights
                counts_ims[exp_id] += h_ims.weights
                counts_res[exp_id] += h_res.weights
            end
    
            # look at change of the parameters
            for id in 2:length(experiments) # for all experiments
                if id == 2
                    change_abs = [abs.(abs_λs[id] - abs_λs[id-1])]
                else
                    push!(change_abs, abs.(abs_λs[id] - abs_λs[id-1]))
                end
            end
    
            # percentual change
            for id in 2:length(experiments) # for all experiments
                if id == 2
                    change_abs_p = [abs.(abs_λs[id] ./ abs_λs[id-1])]
                else
                    push!(change_abs_p, abs.(abs_λs[id] ./ abs_λs[id-1]))
                end
            end
        end
    
        arr_change_abs = combinedims(change_abs, 1)
        arr_change_abs_p = combinedims(change_abs_p, 1)
    
        idcs = findall(x -> x >= 2, arr_change_abs_p)
        for id in idcs
            arr_change_abs_p[id] = 0
        end
    
        if d_part == 1
            # plot the change of some of the imoprtant eigenvalues
            val, idx = findmax(sum(arr_change_abs, dims=1))
            val_p, idx_p = findmax(sum(arr_change_abs_p, dims=1))
    
            # specific plot, to show change of the most changing evs, also just the firs time this is run
            a = plot(abs_ls_sort[:, idx[2], idx[3]], xlabel="experiment", ylabel="b.E.")
            b = plot(abs_ls_sort[:, idx_p[2], idx_p[3]], xlabel="experiment", ylabel="%")
    
            p = plot(a, b, layout=2,
                size=(1000, 500),
                legend=nothing,
                margin=10mm,
                title=["absolute change" "percentual change"],
                plot_title="detailed changes of absolute values", plot_titlevspan=0.05)
            savefig(d, "../Figures/bursting_neuron_all/ev_t_detail.png")
        end
    
        if d_part == 1
            h_abs_all = fit(Histogram, arr_change_abs, nbins=bins)
            h_abs_p = fit(Histogram, arr_change_abs_p, nbins=bins)
    
            global counts_abs_p = h_abs_p.weights
            global counts_abs_all = h_abs_all.weights
            global edges_abs_p = h_abs_p.edges[1]
            global edges_abs_all = h_abs_all.edges[1]
        else
            h_abs_all = fit(Histogram, arr_change_abs, nbins=bins)
            h_abs_p = fit(Histogram, arr_change_abs_p, nbins=bins)
    
            counts_abs_p += h_abs_p.weights
            counts_abs_all += h_abs_all.weights
        end
    end

    HIST_abs_p = Histogram(edges_abs_p, counts_abs_p)
    HIST_abs_all = Histogram(edges_abs_all, counts_abs_all)

    HISTS_abs = [Histogram(edges_[i], counts_[i]) for i in 1:length(experiments)]
    HISTS_ims = [Histogram(edges_[i], counts_[i]) for i in 1:length(experiments)]
    HISTS_res = [Histogram(edges_[i], counts_[i]) for i in 1:length(experiments)]

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

# for debugging reasons
if pwd() != "/home/patrick/.julia/dev/empirical_evd/eval"
    cd("eval")
end
