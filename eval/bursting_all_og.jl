
include("../src/empirical_evd.jl")
using .empirical_evd
using Plots
using BPTT
using LinearAlgebra
using Measures
using SplitApplyCombine

bins = 100
exp_str = "bursting_neuron"
dir = "../Figures/$(exp_str)_all"
function remove!(a, item)
    deleteat!(a, findall(x -> x == item, a))
end

experiments = readdir("../Results")
remove!(experiments, "default")
deleteat!(experiments, findall(x -> !occursin(exp_str, x), experiments))

# load all the experiments
for exp in experiments
    model_path = load_result_path(exp, 5000, 1)
    m = BPTT.load_model(model_path)

    # copy all images
    mkpath(dir)
    dst_path = "$dir/$(exp)_traj.png"
    # cp(load_result_path(exp, 5000, 1; img=true), dst_path, force=true)

    # safe the latent dim of the plrnn
    dim = length(m.A)

    # loop over all possible matrices D
    all_Ds = ([0, 1] for i in 1:dim)
    global all_λs = nothing
    global all_λs = Vector{Vector{ComplexF32}}()
    for (i, D) in enumerate(collect(Iterators.product(all_Ds...)))
        Dt = diagm(collect(D))

        # compose single matrix A+WD
        AWD = diagm(m.A) + m.W * Dt

        # get evd of AWD
        EVD = eigen(AWD)
        λs = EVD.values
        Vs = EVD.vectors
        # add the same for Vs by doing the same with labda->v
        push!(all_λs, complex(λs))
    end
    try
        push!(exp_λs, all_λs)
    catch e
        if typeof(e) == UndefVarError
            global exp_λs = [all_λs]
        else
            rethrow(e)
        end
    end
end

model_path = load_result_path(experiments[1], 5000, 1)
m = BPTT.load_model(model_path)
dim = length(m.A)
#just a few checks
@assert length(exp_λs) == length(experiments)
@assert typeof(exp_λs) == Vector{Vector{Vector{ComplexF32}}}
@assert [length(exp_λs[i]) for i in 1:length(experiments)] == [2^length(m.A) for i in 1:length(experiments)]

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

# convert to array
arr_l = combinedims(exp_λs)
arr_l = combinedims(arr_l)
exp_λs = splitdims(arr_l)
arr_l = permutedims(arr_l, (3, 2, 1))

# sort all the arrays
abs_ls_sort = zeros(6, 2^dim, dim)
res_ls_sort = zeros(6, 2^dim, dim)
ims_ls_sort = zeros(6, 2^dim, dim)

for (ex, all_lambdas) in enumerate(exp_λs)
    abs_ls = abs.(all_lambdas)
    res_ls = real(all_lambdas)
    ims_ls = imag(all_lambdas)
    v = nothing
    for i in 1:size(all_lambdas)[2]
        v = sortperm(abs_ls[:, i])
        abs_ls_sort[ex, i, :] = abs_ls[v, i]
        res_ls_sort[ex, i, :] = res_ls[v, i]
        ims_ls_sort[ex, i, :] = ims_ls[v, i]
    end
end

# Plots of single experiments
# make it a vector for easy iteration
abs_λs = splitdims(abs_ls_sort, 1)
ims_λs = splitdims(ims_ls_sort, 1)
res_λs = splitdims(res_ls_sort, 1)


# histogram of all abs values
hists_abs_val = (histogram(collect(Iterators.flatten(abs_val)), yaxis=:log, bins=bins, legend=nothing, xlabel="Ev value [b.E.]", ylabel="Hits") for abs_val in abs_λs)

# im and real part histograms
hists_ims = (histogram(collect(Iterators.flatten(ims_val)), yaxis=:log, bins=bins, xlabel="Ev value [b.E.]") for ims_val in ims_λs)
hists_res = (histogram(collect(Iterators.flatten(res_val)), yaxis=:log, bins=bins, xlabel="Ev value [b.E.]") for res_val in res_λs)

make_subplots(exp_str, hists_abs_val, "histogram of absolute values")
make_subplots(exp_str, hists_ims, "histogram of imaginary part")
make_subplots(exp_str, hists_res, "historam of real parts")

_, index = findmax(abs_λs[1])

# abs values together for one case
if exp_str == "lorenz"
    c = [:red, :red, :turquoise, :green, :black, :black]
    label = ["limit cycle", "limit cycle", "limit cycle", "right after bifurcation", "chaotic regime", "chaotic regime"]
else
    c = [:orange, :red, :darkred, :purple, :black]
    label = ["cycle", "1 loop", "4 loops", "almost bursting", "bursting"]
end
for (i, abs_val) in enumerate(abs_λs)
    try
        plot!(abs_plot, abs_val[index[1], :], 1:dim, c=c[i], xlabel="Ev value [b.E.]",
            ylabel="Ev No", label=label[i])
    catch
        global abs_plot = plot(abs_val[index[1], :], 1:dim,
            c=c[i],
            label=label[i],
            xlabel="Ev value [b.E.]",
            ylabel="Ev No",
            plot_title="Eigenvalue evolution of quadrant $(index[1])", plot_titlevspan=0.05)
    end
end
savefig(abs_plot, "../Figures/lorenz_all/change_in_abs_vals.png")

# look at change of the parameters
for id in 2:length(experiments) # for all experiments
    try
        push!(change_abs, abs.(abs_λs[id] - abs_λs[id-1]))
    catch e
        global change_abs = [abs.(abs_λs[id] - abs_λs[id-1])]
    end
end

# percentual change
for id in 2:length(experiments) # for all experiments
    try
        push!(change_abs_p, abs.(abs_λs[id] ./ abs_λs[id-1]))
    catch e
        global change_abs_p = [abs.(abs_λs[id] ./ abs_λs[id-1])]
    end
end
arr_change_abs = combinedims(change_abs, 1)
arr_change_abs_p = combinedims(change_abs_p, 1)

idcs = findall(x -> x >= 2, arr_change_abs_p)
for id in idcs
    arr_change_abs_p[id] = 0
end

# hist of overall change in all eigenvalues
hist_abs_p = histogram(collect(Iterators.flatten(arr_change_abs_p .- 1)),
    title="percentual change of all eigenvalues",
    yaxis=:log,
    legend=nothing, bins=bins)

# hist of overall change in all eigenvalues absou change
hist_abs = histogram(collect(Iterators.flatten(arr_change_abs)),
    yaxis=:log,
    title="absolute change of all eigenvalues",
    legend=nothing, bins=bins)


savefig(hist_abs_p, "../Figures/lorenz_all/overall_cahnge_abs_val_p.png")
savefig(hist_abs, "../Figures/lorenz_all/overall_cahnge_abs_val.png")

# plot the change of some of the imoprtant eigenvalues
val, idx = findmax(sum(arr_change_abs, dims=1))
val_p, idx_p = findmax(sum(arr_change_abs_p, dims=1))

# specific plot, to show change of the most changing evs
a = plot(abs_ls_sort[:, idx[2], idx[3]], xlabel="experiment", ylabel="b.E.")
b = plot(abs_ls_sort[:, idx_p[2], idx_p[3]], xlabel="experiment", ylabel="%")

d = plot(a, b, layout=2,
    size=(1000, 500),
    legend=nothing,
    margin=10mm,
    title=["most absolute change" "most percentual change"],
    plot_title="detailed changes of absolute values", plot_titlevspan=0.05)
savefig(d, "../Figures/lorenz_all/ev_t_detail.png")
