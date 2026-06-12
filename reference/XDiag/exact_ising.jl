using XDiag
using LinearAlgebra
using Printf
using Random
using CairoMakie
using MakiePublication

const PAULI_LABELS = ("I", "X", "Y", "Z")

function apply_pauli_op(psi, site::Int, pauli_idx::Int)
    if pauli_idx == 1 # I
        return psi
    elseif pauli_idx == 2 # X = S+ + S-
        return apply(Op("S+", site), psi) + apply(Op("S-", site), psi)
    elseif pauli_idx == 3 # Y = -i(S+ - S-)
        return apply(-1.0im * Op("S+", site), psi) + apply(1.0im * Op("S-", site), psi)
    elseif pauli_idx == 4 # Z = 2Sz
        return apply(2.0 * Op("Sz", site), psi)
    else
        error("Invalid Pauli index")
    end
end

function compute_state_from_indices(psi0, indices)
    state = psi0
    for (site, p_idx) in enumerate(indices)
        if p_idx != 1
            state = apply_pauli_op(state, site, p_idx)
        end
    end
    return state
end

function build_ising_ops(N::Int, h::Float64; J::Float64 = 1.0)
    ops = OpSum()
    for i in 1:N
        j = mod1(i + 1, N)
        ops += -4.0 * J * Op("SzSz", [i, j])
        ops += -1.0 * h * Op("S+", [i])
        ops += -1.0 * h * Op("S-", [i])
    end
    return ops
end

function ground_state_ising(N::Int, h::Float64; J::Float64 = 1.0)
    block = Spinhalf(N)
    ops = build_ising_ops(N, h; J = J)
    energy, psi0 = eig0(ops, block)
    return energy, psi0
end

function run_exact_ising(; N::Int = 8, h::Float64 = 1.0, verbose::Bool = true)
    verbose && println("Running Exact Enumeration for N=$N, h=$h...")

    _, psi0 = ground_state_ising(N, h)

    stabilizer_sum = 0.0
    n_paulis = 4^N

    # Stabilizer Rényi entropy: M₂ = -log(2^-N * Σ_P |<ψ|P|ψ>|^4).
    indices = ones(Int, N)

    for k in 0:(n_paulis - 1)
        temp_k = k
        for i in 1:N
            indices[i] = (temp_k % 4) + 1
            temp_k ÷= 4
        end

        state = compute_state_from_indices(psi0, indices)
        expectation = dot(psi0, state)
        stabilizer_sum += abs2(abs2(expectation))
    end

    magic_density = -log(stabilizer_sum / 2^N) / N

    if verbose
        println("\nStabilizer Sum: $stabilizer_sum")
        println("2^-N Normalized Sum: $(stabilizer_sum / 2^N)")
        println("Exact Magic Density M₂/N: $magic_density")
    end

    return magic_density
end

function central_derivative(xs, ys)
    n = length(xs)
    dydx = similar(ys)
    dydx[1] = (ys[2] - ys[1]) / (xs[2] - xs[1])
    dydx[n] = (ys[n] - ys[n - 1]) / (xs[n] - xs[n - 1])
    for i in 2:(n - 1)
        dydx[i] = (ys[i + 1] - ys[i - 1]) / (xs[i + 1] - xs[i - 1])
    end
    return dydx
end

function vector_ground_state_ising(N::Int, h::Float64; J::Float64 = 1.0)
    _, psi0 = ground_state_ising(N, h; J = J)
    return ComplexF64.(vector(psi0))
end

function pauli_weight(psi::Vector{ComplexF64}, p_idx::Integer, N::Int)::Float64
    dim = length(psi)
    x_mask = 0
    z_mask = 0
    y_count = 0
    temp = p_idx

    for site in 0:(N - 1)
        code = temp & 0x3
        temp >>= 2
        if code == 1
            x_mask |= 1 << site
        elseif code == 2
            x_mask |= 1 << site
            z_mask |= 1 << site
            y_count += 1
        elseif code == 3
            z_mask |= 1 << site
        end
    end

    y_phase = (-1.0im)^y_count
    expectation = 0.0 + 0.0im
    @inbounds for basis in 0:(dim - 1)
        flipped = basis ⊻ x_mask
        sign = isodd(count_ones(basis & z_mask)) ? -1.0 : 1.0
        expectation += conj(psi[basis + 1]) * psi[flipped + 1] * sign * y_phase
    end
    return abs2(expectation)
end

function pauli_string_label(p_idx::Integer, N::Int)::String
    labels = Vector{String}(undef, N)
    temp = p_idx
    for site in 1:N
        labels[site] = PAULI_LABELS[(temp & 0x3) + 1]
        temp >>= 2
    end
    return join(labels)
end

function pauli_string_weights(psi::Vector{ComplexF64}, N::Int)
    labels = String[]
    weights = Float64[]
    contributions = Float64[]
    for p_idx in 0:(4^N - 1)
        w = pauli_weight(psi, p_idx, N)
        push!(labels, pauli_string_label(p_idx, N))
        push!(weights, w)
        push!(contributions, w^2)
    end
    return labels, weights, contributions
end

function plot_pauli_string_histogram(; N::Int = 2, h::Float64 = 1.0, J::Float64 = 1.0, outdir = "fig")
    isdir(outdir) || mkpath(outdir)

    energy, psi0 = ground_state_ising(N, h; J = J)
    psi = ComplexF64.(vector(psi0))
    labels, weights, contributions = pauli_string_weights(psi, N)

    println("Pauli string histogram for N=$N, h=$h, J=$J")
    println("Ground-state energy: $energy")
    @printf "%6s  %18s  %18s\n" "P" "|<P>|^2" "|<P>|^4"
    for (label, weight, contribution) in zip(labels, weights, contributions)
        @printf "%6s  %18.12g  %18.12g\n" label weight contribution
    end
    @printf "stabilizer_sum = %.12g\n" sum(contributions)
    @printf "M2/N = %.12g\n" (-log(sum(contributions) / 2^N) / N)

    with_theme(theme_web()) do
        fig = Figure(size = (980, 460), backgroundcolor = RGBf(0.985, 0.982, 0.975))
        ax = Axis(
            fig[1, 1],
            xlabel = "Pauli string",
            ylabel = L"|\langle P \rangle|^2",
            title = "TFIM Pauli-string weights, N = $N, h = $h, J = $J",
            xticks = (1:length(labels), labels),
            xticklabelrotation = pi / 3,
            backgroundcolor = RGBf(0.985, 0.982, 0.975),
            xgridcolor = (:gray30, 0.08),
            ygridcolor = (:gray30, 0.12),
            topspinevisible = false,
            rightspinevisible = false
        )
        colors = [cgrad(:magma, 256)[0.25 + 0.65 * sqrt(clamp(w, 0.0, 1.0))] for w in weights]
        barplot!(
            ax,
            1:length(weights),
            weights;
            color = colors,
            strokecolor = (:black, 0.35),
            strokewidth = 0.8,
            gap = 0.25
        )
        ylims!(ax, 0, 1.08 * maximum(weights))
        for (idx, weight) in enumerate(weights)
            weight > 1e-10 || continue
            text!(
                ax,
                idx,
                weight + 0.035 * maximum(weights);
                text = @sprintf("%.2g", weight),
                align = (:center, :bottom),
                fontsize = 11,
                color = (:black, 0.75)
            )
        end
        for ext in ("pdf", "png")
            path = joinpath(outdir, "pauli_string_histogram_N$(N)_h$(h_cache_key(h))_J$(h_cache_key(J)).$ext")
            save(path, fig)
            println("Saved $path")
        end
    end

    return labels, weights, contributions
end

function exact_magic_density_bitwise(psi::Vector{ComplexF64}, N::Int)::Float64
    stabilizer_sum = 0.0
    for p_idx in 0:(4^N - 1)
        w = pauli_weight(psi, p_idx, N)
        stabilizer_sum += w^2
    end
    return -log(stabilizer_sum / 2^N) / N
end

function mcmc_magic_density(
    psi::Vector{ComplexF64},
    N::Int;
    samples::Int = 100_000,
    thermalization::Int = 10_000,
    rng::AbstractRNG = Random.default_rng(),
    checkpoints::Vector{Int} = Int[]
)
    current_p = rand(rng, 0:(4^N - 1))
    current_weight = pauli_weight(psi, current_p, N)
    while current_weight < 1e-14
        current_p = rand(rng, 0:(4^N - 1))
        current_weight = pauli_weight(psi, current_p, N)
    end

    checkpoint_set = Set(checkpoints)
    checkpoint_values = Dict{Int, Float64}()
    weight_sum = 0.0
    n_acc = 0

    for step in 1:(thermalization + samples)
        proposed_p = rand(rng, 0:(4^N - 1))
        proposed_weight = pauli_weight(psi, proposed_p, N)
        ratio = proposed_weight / current_weight
        if rand(rng) < min(1.0, ratio)
            current_p = proposed_p
            current_weight = proposed_weight
            step > thermalization && (n_acc += 1)
        end

        if step > thermalization
            sample_idx = step - thermalization
            weight_sum += current_weight
            if sample_idx in checkpoint_set
                checkpoint_values[sample_idx] = -log(weight_sum / sample_idx) / N
            end
        end
    end

    magic_density = -log(weight_sum / samples) / N
    return magic_density, n_acc / samples, checkpoint_values
end

function mcmc_magic_density(psi::Vector{ComplexF64}, N::Int, samples::Int, thermalization::Int; rng::AbstractRNG = Random.default_rng())
    magic_density, _, _ = mcmc_magic_density(psi, N; samples = samples, thermalization = thermalization, rng = rng)
    return magic_density
end

function tfim_magic_benchmark(; samples::Int = 50_000, thermalization::Int = 5_000, outdir = "fig")
    isdir(outdir) || mkpath(outdir)

    hs = collect(0.5:0.1:1.5)
    exact8 = Float64[]
    mcmc8 = Float64[]
    mcmc12 = Float64[]

    Random.seed!(20260518)
    for h in hs
        psi8 = vector_ground_state_ising(8, h)
        push!(exact8, exact_magic_density_bitwise(psi8, 8))
        push!(mcmc8, mcmc_magic_density(psi8, 8, samples, thermalization))

        psi12 = vector_ground_state_ising(12, h)
        push!(mcmc12, mcmc_magic_density(psi12, 12, samples, thermalization))
    end

    with_theme(theme_web()) do
        fig = Figure()
        ax = Axis(
            fig[1, 1],
            xlabel = L"h/J",
            ylabel = L"M_2/N",
            title = "MCMC stabilizer Rényi magic density"
        )

        plots = Any[]
        labels = String[]
        push!(plots, scatterlines!(ax, hs, exact8))
        push!(labels, "N = 8 exact")
        push!(plots, scatterlines!(ax, hs, mcmc8))
        push!(labels, "N = 8 MCMC")
        push!(plots, scatterlines!(ax, hs, mcmc12))
        push!(labels, "N = 12 MCMC")
        vlines!(ax, [1.0], color = :gray, linestyle = :dash)
        add_axislegend_if_labeled!(ax, plots, labels, position = :lt)

        for ext in ("pdf", "png")
            path = joinpath(outdir, "tfim_magic_benchmark.$ext")
            save(path, fig)
            println("Saved $path")
        end
    end

    return nothing
end

function h_cache_key(h::Real; digits::Int = 12)
    return round(Float64(h); digits = digits)
end

function exact_magic_cache_entry_key(N::Int, h::Real)
    return "N=$(N),h=$(h_cache_key(h))"
end

function parse_exact_magic_cache_line(line::AbstractString)
    m = match(r"^\{\"N\":(\d+),\"h\":([^,]+),\"magic_density\":([^}]+)\}$", strip(line))
    m === nothing && return nothing

    N = parse(Int, m.captures[1])
    h = h_cache_key(parse(Float64, m.captures[2]))
    value = parse(Float64, m.captures[3])
    return N, h, value
end

function load_exact_magic_cache(path::AbstractString, Ns)
    cache = Dict{Tuple{Int, Float64}, Float64}()
    isfile(path) || return cache

    for line in eachline(path)
        isempty(strip(line)) && continue
        parsed = parse_exact_magic_cache_line(line)
        parsed === nothing && continue
        N, h, value = parsed
        N in Ns || continue
        cache[(N, h)] = value
    end

    println("Loaded cache from $path")
    return cache
end

function append_exact_magic_cache!(path::AbstractString, new_values::Dict{Tuple{Int, Float64}, Float64})
    isempty(new_values) && return nothing
    open(path, "a") do io
        for ((N, h), value) in sort(collect(new_values); by = item -> (item[1][1], item[1][2]))
            println(io, "{\"N\":$N,\"h\":$h,\"magic_density\":$value}")
        end
    end
    println("Updated cache $path")
    return path
end

function save_exact_magic_csv(path::AbstractString, hs, magic::Dict{Int, Vector{Float64}}, Ns)
    open(path, "w") do io
        print(io, "h")
        for N in Ns
            print(io, ",N$(N)_magic_density")
        end
        println(io)

        for (idx, h) in enumerate(hs)
            print(io, h)
            for N in Ns
                print(io, ",", magic[N][idx])
            end
            println(io)
        end
    end
    println("Saved $path")
    return path
end

function exact_magic_sweep(; Ns = [4, 6, 8], hs = collect(0.5:0.05:1.5), outdir = nothing, csv_name = "exact_ising_magic.csv", cache_name = "exact_ising_magic_cache.jsonl", use_cache::Bool = true)
    hs = collect(Float64, hs)
    csv_path = outdir === nothing ? nothing : joinpath(outdir, csv_name)
    cache_path = outdir === nothing ? nothing : joinpath(outdir, cache_name)
    cache = use_cache && cache_path !== nothing ? load_exact_magic_cache(cache_path, Ns) : Dict{Tuple{Int, Float64}, Float64}()
    new_cache_values = Dict{Tuple{Int, Float64}, Float64}()

    magic = Dict{Int, Vector{Float64}}()
    for N in Ns
        vals = Float64[]
        println("\n=== Exact transverse-field sweep for N=$N ===")
        for h in hs
            key = (N, h_cache_key(h))
            if use_cache && haskey(cache, key)
                val = cache[key]
                @printf "N=%2d  h=%5.2f  magic_density=%10.7f  cache\n" N h val
            else
                val = run_exact_ising(N = N, h = Float64(h), verbose = false)
                new_cache_values[key] = val
                @printf "N=%2d  h=%5.2f  magic_density=%10.7f  computed\n" N h val
            end
            push!(vals, val)
        end
        magic[N] = vals
    end
    if csv_path !== nothing
        isdir(outdir) || mkpath(outdir)
        cache_path !== nothing && append_exact_magic_cache!(cache_path, new_cache_values)
        save_exact_magic_csv(csv_path, hs, magic, Ns)
    end
    return hs, magic
end

function add_axislegend_if_labeled!(ax::Axis, plots, labels; kwargs...)
    isempty(labels) || axislegend(ax, plots, labels; kwargs...)
    return nothing
end

function make_thesis_exact_ising_figures(; Ns = [4, 6, 8], hs = collect(0.5:0.05:1.5), outdir = "fig", use_cache::Bool = true)
    isdir(outdir) || mkpath(outdir)
    hs, magic = exact_magic_sweep(Ns = Ns, hs = hs, outdir = outdir, use_cache = use_cache)

    with_theme(theme_web()) do
        fig = Figure()
        ax = Axis(
            fig[1, 1],
            xlabel = L"h/J",
            ylabel = L"M_2/N",
            title = "Exact stabilizer Rényi magic density"
        )

        plots = Any[]
        labels = String[]
        for N in Ns
            push!(plots, scatterlines!(ax, hs, magic[N]))
            push!(labels, "N = $N")
        end
        vlines!(ax, [1.0], color = :gray, linestyle = :dash)
        add_axislegend_if_labeled!(ax, plots, labels, position = :lt)

        for ext in ("pdf", "png")
            path = joinpath(outdir, "exact_ising_magic_thesis.$ext")
            save(path, fig)
            println("Saved $path")
        end
    end

    return hs, magic
end

if abspath(PROGRAM_FILE) == @__FILE__
    make_thesis_exact_ising_figures()
end
