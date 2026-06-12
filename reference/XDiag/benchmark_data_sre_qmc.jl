include("exact_ising.jl")

function load_l8_exact_from_sweep_csv(path::AbstractString)
    isfile(path) || return nothing

    lines = readlines(path)
    isempty(lines) && return nothing
    header = split(lines[1], ",")
    h_col = findfirst(==("h"), header)
    val_col = findfirst(==("N8_magic_density"), header)
    (h_col === nothing || val_col === nothing) && return nothing

    hs = Float64[]
    vals = Float64[]
    for line in lines[2:end]
        isempty(strip(line)) && continue
        fields = split(line, ",")
        length(fields) >= max(h_col, val_col) || continue
        push!(hs, parse(Float64, fields[h_col]))
        push!(vals, parse(Float64, fields[val_col]))
    end

    println("Loaded N=8 exact data from $path")
    return hs, vals
end

function load_l8_exact_from_jsonl_cache(path::AbstractString)
    cache = load_exact_magic_cache(path, [8])
    isempty(cache) && return nothing

    pairs = sort([(h, value) for ((N, h), value) in cache if N == 8]; by = first)
    isempty(pairs) && return nothing
    hs = [p[1] for p in pairs]
    vals = [p[2] for p in pairs]
    println("Loaded N=8 exact data from $path")
    return hs, vals
end

function load_l8_exact_cache(; outdir = "fig", cache_name = "exact_ising_magic_cache.jsonl", csv_name = "exact_ising_magic.csv")
    jsonl_data = load_l8_exact_from_jsonl_cache(joinpath(outdir, cache_name))
    jsonl_data !== nothing && return jsonl_data
    return load_l8_exact_from_sweep_csv(joinpath(outdir, csv_name))
end

function load_l8_benchmark_csv(path::AbstractString)
    isfile(path) || return nothing

    lines = readlines(path)
    isempty(lines) && return nothing
    split(lines[1], ",") == ["J_over_h", "exact_L8_magic_density", "qmc_L8_magic_density", "qmc_L8_error"] || return nothing

    exact_J_over_hs = Float64[]
    exact_vals = Float64[]
    qmc_vals = Float64[]
    qmc_errs = Float64[]
    for line in lines[2:end]
        isempty(strip(line)) && continue
        fields = split(line, ",")
        length(fields) == 4 || return nothing
        push!(exact_J_over_hs, parse(Float64, fields[1]))
        push!(exact_vals, parse(Float64, fields[2]))
        push!(qmc_vals, parse(Float64, fields[3]))
        push!(qmc_errs, parse(Float64, fields[4]))
    end

    println("Loaded $path")
    return exact_J_over_hs, exact_vals, qmc_vals, qmc_errs
end

function thermal_ising_magic_density(N::Int, beta::Float64; h::Float64 = 1.0, J::Float64 = 1.0)
    block = Spinhalf(N)
    ops = build_ising_ops(N, h; J = J)
    H = Hermitian(Matrix(matrix(ops, block)))
    spectrum = eigen(H)
    shifted_energies = spectrum.values .- minimum(spectrum.values)
    boltzmann = exp.(-beta .* shifted_energies)
    rho = spectrum.vectors * Diagonal(boltzmann / sum(boltzmann)) * spectrum.vectors'

    stabilizer_sum = 0.0
    dim = size(rho, 1)
    for p_idx in 0:(4^N - 1)
        x_mask = 0
        z_mask = 0
        temp = p_idx
        for site in 0:(N - 1)
            code = temp & 0x3
            temp >>= 2
            if code == 1
                x_mask |= 1 << site
            elseif code == 2
                x_mask |= 1 << site
                z_mask |= 1 << site
            elseif code == 3
                z_mask |= 1 << site
            end
        end

        expectation = 0.0 + 0.0im
        @inbounds for basis in 0:(dim - 1)
            flipped = basis ⊻ x_mask
            sign = isodd(count_ones(basis & z_mask)) ? -1.0 : 1.0
            expectation += rho[basis + 1, flipped + 1] * sign
        end
        stabilizer_sum += abs2(expectation)^2
    end

    return -log(stabilizer_sum / 2^N) / N
end

function load_l8_thermal_benchmark_csv(path::AbstractString)
    isfile(path) || return nothing

    lines = readlines(path)
    isempty(lines) && return nothing
    split(lines[1], ",") == ["beta_over_L", "exact_beta", "exact_L8_magic_density", "qmc_L8_magic_density", "qmc_L8_error"] || return nothing

    beta_over_Ls = Float64[]
    exact_betas = Float64[]
    exact_vals = Float64[]
    qmc_vals = Float64[]
    qmc_errs = Float64[]
    for line in lines[2:end]
        isempty(strip(line)) && continue
        fields = split(line, ",")
        length(fields) == 5 || return nothing
        push!(beta_over_Ls, parse(Float64, fields[1]))
        push!(exact_betas, parse(Float64, fields[2]))
        push!(exact_vals, parse(Float64, fields[3]))
        push!(qmc_vals, parse(Float64, fields[4]))
        push!(qmc_errs, parse(Float64, fields[5]))
    end

    println("Loaded $path")
    return beta_over_Ls, exact_betas, exact_vals, qmc_vals, qmc_errs
end

function read_number_list(path::AbstractString)
    return parse.(Float64, readlines(path))
end

function find_or_clone_data_repo(; repo_url = "git@github.com:ymd-physics/Data-SRE-QMC.git", data_dir = joinpath(@__DIR__, "Data-SRE-QMC"))
    if isdir(data_dir)
        return data_dir
    end

    parent = dirname(data_dir)
    isdir(parent) || mkpath(parent)
    run(`git clone --depth 1 $repo_url $data_dir`)
    return data_dir
end

function compare_exact_ising_with_data_sre_qmc(; data_dir = joinpath(@__DIR__, "Data-SRE-QMC"), outdir = "fig", use_cache::Bool = true, use_exact_cache_only::Bool = true)
    isdir(outdir) || mkpath(outdir)
    csv_path = joinpath(outdir, "exact_ising_vs_data_sre_qmc_L8.csv")
    cached = use_cache ? load_l8_benchmark_csv(csv_path) : nothing

    if cached === nothing
        data_dir = find_or_clone_data_repo(data_dir = data_dir)

        qmc_dir = joinpath(data_dir, "SE_1D")
        qmc_J_over_hs_full = read_number_list(joinpath(qmc_dir, "L8_J_list.dat"))
        qmc_vals_full = read_number_list(joinpath(qmc_dir, "L8_val_list.dat"))
        qmc_errs_full = read_number_list(joinpath(qmc_dir, "L8_err_list.dat"))

        exact_J_over_hs = collect(0.8:0.02:1.2)
        exact_vals = [exact_magic_density_bitwise(vector_ground_state_ising(8, 1.0 / J_over_h), 8) for J_over_h in exact_J_over_hs]
        qmc_J_over_hs = Float64[]
        qmc_vals = Float64[]
        qmc_errs = Float64[]
        for J_over_h in exact_J_over_hs
            idx = argmin(abs.(qmc_J_over_hs_full .- J_over_h))
            push!(qmc_J_over_hs, qmc_J_over_hs_full[idx])
            push!(qmc_vals, qmc_vals_full[idx])
            push!(qmc_errs, qmc_errs_full[idx])
        end
        open(csv_path, "w") do io
            println(io, "J_over_h,exact_L8_magic_density,qmc_L8_magic_density,qmc_L8_error")
            for (J_over_h, exact_val, qmc_val, qmc_err) in zip(exact_J_over_hs, exact_vals, qmc_vals, qmc_errs)
                println(io, "$(J_over_h),$(exact_val),$(qmc_val),$(qmc_err)")
            end
        end
        println("Saved $csv_path")
    else
        exact_J_over_hs, exact_vals, qmc_vals, qmc_errs = cached
        qmc_J_over_hs = exact_J_over_hs
    end

    with_theme(theme_web()) do
        fig = Figure()
        ax = Axis(
            fig[1, 1],
            xlabel = L"J/h",
            ylabel = L"M_2/N",
            title = "Exact Ising vs Data-SRE-QMC, N = 8"
        )

        plots = Any[]
        labels = String[]
        push!(plots, lines!(ax, qmc_J_over_hs, qmc_vals))
        push!(labels, "Data-SRE-QMC L=8")
        errorbars!(ax, qmc_J_over_hs, qmc_vals, qmc_errs)
        push!(plots, scatterlines!(ax, exact_J_over_hs, exact_vals))
        push!(labels, "XDiag exact L=8")
        vlines!(ax, [1.0], color = :gray, linestyle = :dash)
        add_axislegend_if_labeled!(ax, plots, labels, position = :lt)

        for ext in ("pdf", "png")
            path = joinpath(outdir, "exact_ising_vs_data_sre_qmc_L8.$ext")
            save(path, fig)
            println("Saved $path")
        end
    end

    return exact_J_over_hs, exact_vals, qmc_J_over_hs, qmc_vals, qmc_errs
end

function compare_thermal_ising_with_data_sre_qmc(; data_dir = joinpath(@__DIR__, "Data-SRE-QMC"), outdir = "fig", L::Int = 8, beta::Real = L, use_cache::Bool = true)
    L == 8 || error("Only L=8 is wired to the available benchmark data in this script")
    isdir(outdir) || mkpath(outdir)

    beta = Float64(beta)
    beta_over_L = beta / L
    csv_path = joinpath(outdir, "thermal_ising_vs_data_sre_qmc_L$(L)_beta$(h_cache_key(beta)).csv")
    cached = use_cache ? load_l8_thermal_benchmark_csv(csv_path) : nothing

    if cached === nothing
        data_dir = find_or_clone_data_repo(data_dir = data_dir)
        qmc_dir = joinpath(data_dir, "SE_Thermal_2D")
        qmc_beta_over_Ls_full = read_number_list(joinpath(qmc_dir, "L$(L)_beta_list.dat"))
        qmc_vals_full = read_number_list(joinpath(qmc_dir, "L$(L)_val_list.dat"))
        qmc_errs_full = read_number_list(joinpath(qmc_dir, "L$(L)_err_list.dat"))

        qmc_idx = argmin(abs.(qmc_beta_over_Ls_full .- beta_over_L))
        qmc_beta_over_L = qmc_beta_over_Ls_full[qmc_idx]
        qmc_val = qmc_vals_full[qmc_idx]
        qmc_err = qmc_errs_full[qmc_idx]
        exact_val = thermal_ising_magic_density(L, beta)

        open(csv_path, "w") do io
            println(io, "beta_over_L,exact_beta,exact_L8_magic_density,qmc_L8_magic_density,qmc_L8_error")
            println(io, "$(qmc_beta_over_L),$(beta),$(exact_val),$(qmc_val),$(qmc_err)")
        end
        println("Saved $csv_path")
        beta_over_Ls = [qmc_beta_over_L]
        exact_betas = [beta]
        exact_vals = [exact_val]
        qmc_vals = [qmc_val]
        qmc_errs = [qmc_err]
    else
        beta_over_Ls, exact_betas, exact_vals, qmc_vals, qmc_errs = cached
    end

    with_theme(theme_web()) do
        fig = Figure()
        ax = Axis(
            fig[1, 1],
            xlabel = L"\beta/L",
            ylabel = L"M_2/N",
            title = "Thermal Ising vs Data-SRE-QMC, L = $L, beta = $beta"
        )

        plots = Any[]
        labels = String[]
        push!(plots, scatter!(ax, beta_over_Ls, qmc_vals))
        push!(labels, "Data-SRE-QMC L=$L")
        errorbars!(ax, beta_over_Ls, qmc_vals, qmc_errs)
        push!(plots, scatter!(ax, beta_over_Ls, exact_vals))
        push!(labels, "XDiag thermal exact beta=$(only(exact_betas))")
        add_axislegend_if_labeled!(ax, plots, labels, position = :lt)

        for ext in ("pdf", "png")
            path = joinpath(outdir, "thermal_ising_vs_data_sre_qmc_L$(L)_beta$(h_cache_key(beta)).$ext")
            save(path, fig)
            println("Saved $path")
        end
    end

    return beta_over_Ls, exact_betas, exact_vals, qmc_vals, qmc_errs
end

if abspath(PROGRAM_FILE) == @__FILE__
    compare_exact_ising_with_data_sre_qmc()
end
