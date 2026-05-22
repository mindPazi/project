"""
Phase 5: Numerical vs theoretical analysis.
Compares Monte Carlo results against analytical predictions.
"""

include("src/network.jl")
include("src/swapping.jl")
include("src/metrics.jl")

using .Network, .Swapping, .Metrics
using Plots, Statistics, Printf
using Distributions: Geometric

const M_RUNS = 1000

# --- Theoretical formulas ---

harmonic(n) = sum(1.0/k for k in 1:n)

"""
Exact expected value of the maximum of n i.i.d. Geometric(p) random variables.
E[max] = Σ_{t=0}^{∞} [1 - (1 - (1-p)^t)^n]
"""
function expected_max_geometric(n, p; max_terms=5000)
    q = 1 - p
    s = 0.0
    for t in 0:max_terms
        term = 1.0 - (1.0 - q^t)^n
        s += term
        term < 1e-12 && break
    end
    s
end

"""
Approximate fidelity from the project plan: F ≈ (1 - 4/3 p_w)^(T·N).
Rough worst-case: assumes all qubits decohere for the full T.
"""
fidelity_approx(N, T, p_w) = max(0.25, (1 - 4/3 * p_w) ^ (T * N))

"""
Werner-state fidelity for a chain of n_links = N+1 links.
Each link i has both qubits depolarized for time wait_i = T - t_i.
Single-qubit depol parameter: p(Δt) = 1 - (1-p_w)^Δt.
Werner parameter per link: w_i = (1 - 4/3 p(wait_i))^2  (two independent qubits).
Swapping multiplies Werner parameters: w_total = ∏ w_i.
F = (1 + 3·w_total) / 4.
"""
function fidelity_werner(gen_times, p_w)
    T = maximum(gen_times)
    w_total = 1.0
    for t_i in gen_times
        wait = T - t_i
        p_depol = 1.0 - (1.0 - p_w)^wait
        w_link = (1.0 - 4.0/3.0 * p_depol)^2
        w_total *= w_link
    end
    (1.0 + 3.0 * w_total) / 4.0
end

"""
Run M single runs and return vectors of (fidelity_mc, fidelity_werner, fidelity_approx, time).
"""
function monte_carlo_with_theory(N::Int, M::Int; p_success=1.0, p_w=0.0)
    n_links = N + 1
    f_mc  = Vector{Float64}(undef, M)
    f_wer = Vector{Float64}(undef, M)
    f_app = Vector{Float64}(undef, M)
    times = Vector{Float64}(undef, M)
    for i in 1:M
        f, t = Metrics.single_run(N; p_success=p_success, p_w=p_w)
        gen_times = p_success < 1.0 ?
            [rand(Geometric(p_success)) + 1 for _ in 1:n_links] :
            ones(Int, n_links)
        T_th = maximum(gen_times)
        f_mc[i]  = f
        f_wer[i] = fidelity_werner(gen_times, p_w)
        f_app[i] = fidelity_approx(N, T_th, p_w)
        times[i] = t
    end
    (f_mc, f_wer, f_app, times)
end

# ============================================================
#  1. Distribution time: MC vs exact vs harmonic approximation
# ============================================================
function analysis_distribution_time()
    println("=" ^60)
    println("  1. DISTRIBUTION TIME: numerical vs theory")
    println("=" ^60)

    p_values = 0.1:0.1:0.9
    N_values = [1, 3, 5]

    println(@sprintf("\n%-4s  %-6s  %10s  %10s  %10s  %8s  %8s",
        "N", "p_s", "MC mean", "Exact", "Harmonic", "Err_ex%", "Err_ha%"))
    println("-" ^68)

    results = Dict()
    for N in N_values
        n_links = N + 1
        for ps in p_values
            _, _, t_mc, _ = Metrics.monte_carlo(N, M_RUNS; p_success=ps, p_w=0.0)
            t_exact = expected_max_geometric(n_links, ps)
            t_harmonic = harmonic(n_links) / ps
            err_ex = abs(t_mc - t_exact) / t_exact * 100
            err_ha = abs(t_mc - t_harmonic) / t_harmonic * 100
            println(@sprintf("N=%d   p=%.1f  %10.3f  %10.3f  %10.3f  %7.2f%%  %7.2f%%",
                N, ps, t_mc, t_exact, t_harmonic, err_ex, err_ha))
            results[(N, ps)] = (t_mc, t_exact, t_harmonic)
        end
    end

    # Plot
    plt = plot(xlabel="p_success", ylabel="Distribution time",
               title="Distribution time: MC vs Theory", legend=:topright)
    colors = [:blue, :red, :green]
    for (idx, N) in enumerate(N_values)
        ps_vec = collect(p_values)
        mc_vals = [results[(N, p)][1] for p in ps_vec]
        exact_vals = [results[(N, p)][2] for p in ps_vec]
        harm_vals = [results[(N, p)][3] for p in ps_vec]
        plot!(plt, ps_vec, mc_vals, marker=:circle, ms=4, label="MC N=$N", color=colors[idx])
        plot!(plt, ps_vec, exact_vals, ls=:solid, lw=2, label="Exact N=$N", color=colors[idx])
        plot!(plt, ps_vec, harm_vals, ls=:dash, lw=1, label="H(N+1)/p N=$N", color=colors[idx], alpha=0.6)
    end
    savefig(plt, "analysis_time_comparison.png")
    println("\n  -> analysis_time_comparison.png\n")
end

# ============================================================
#  2. Fidelity: MC vs Werner (per-run) vs plan approximation
# ============================================================
function analysis_fidelity()
    println("=" ^60)
    println("  2. FIDELITY: MC vs Werner vs plan approx")
    println("=" ^60)

    p_values = 0.1:0.1:1.0
    pw_values = [0.01, 0.05, 0.10]
    N_fixed = 3

    println(@sprintf("\n%-6s  %-6s  %8s  %8s  %8s  %7s  %7s",
        "p_w", "p_s", "MC", "Werner", "Approx", "Ew%", "Ea%"))
    println("-" ^62)

    results = Dict()
    for pw in pw_values
        for ps in p_values
            f_mc, f_wer, f_app, _ = monte_carlo_with_theory(N_fixed, M_RUNS; p_success=ps, p_w=pw)
            fm = mean(f_mc)
            fw = mean(f_wer)
            fa = mean(f_app)
            ew = fm > 0.26 ? abs(fm - fw) / fm * 100 : 0.0
            ea = fm > 0.26 ? abs(fm - fa) / fm * 100 : 0.0
            println(@sprintf("pw=%.2f  p=%.1f  %8.4f  %8.4f  %8.4f  %6.1f%%  %6.1f%%",
                pw, ps, fm, fw, fa, ew, ea))
            results[(pw, ps)] = (fm, fw, fa)
        end
    end

    ps_vec = collect(p_values)
    colors = [:blue, :orange, :green]

    # Plot: MC vs Werner
    plt1 = plot(xlabel="p_success", ylabel="Fidelity",
                title="MC vs Werner (N=$N_fixed)", legend=:bottomright)
    for (idx, pw) in enumerate(pw_values)
        mc_vals = [results[(pw, p)][1] for p in ps_vec]
        wer_vals = [results[(pw, p)][2] for p in ps_vec]
        plot!(plt1, ps_vec, mc_vals, marker=:circle, ms=4, label="MC p_w=$pw", color=colors[idx])
        plot!(plt1, ps_vec, wer_vals, ls=:dash, lw=2, label="Werner p_w=$pw", color=colors[idx])
    end

    # Plot: MC vs plan approx
    plt2 = plot(xlabel="p_success", ylabel="Fidelity",
                title="MC vs Plan approx (N=$N_fixed)", legend=:bottomright)
    for (idx, pw) in enumerate(pw_values)
        mc_vals = [results[(pw, p)][1] for p in ps_vec]
        app_vals = [results[(pw, p)][3] for p in ps_vec]
        plot!(plt2, ps_vec, mc_vals, marker=:circle, ms=4, label="MC p_w=$pw", color=colors[idx])
        plot!(plt2, ps_vec, app_vals, ls=:dot, lw=2, label="Approx p_w=$pw", color=colors[idx])
    end

    plt = plot(plt1, plt2, layout=(1, 2), size=(1200, 450))
    savefig(plt, "analysis_fidelity_comparison.png")
    println("\n  -> analysis_fidelity_comparison.png\n")
end

# ============================================================
#  3. Fidelity scaling with N
# ============================================================
function analysis_scaling_N()
    println("=" ^60)
    println("  3. FIDELITY SCALING WITH N")
    println("=" ^60)

    N_range = 1:7
    ps_fixed = 0.5
    pw_fixed = 0.05

    println(@sprintf("\n%-4s  %8s  %8s  %8s  %7s  %7s",
        "N", "MC F", "Werner", "Approx", "Ew%", "Ea%"))
    println("-" ^52)

    mc_vals = Float64[]
    wer_vals = Float64[]
    app_vals = Float64[]
    for N in N_range
        f_mc, f_wer, f_app, _ = monte_carlo_with_theory(N, M_RUNS; p_success=ps_fixed, p_w=pw_fixed)
        fm = mean(f_mc)
        fw = mean(f_wer)
        fa = mean(f_app)
        ew = fm > 0.26 ? abs(fm - fw) / fm * 100 : 0.0
        ea = fm > 0.26 ? abs(fm - fa) / fm * 100 : 0.0
        println(@sprintf("N=%d   %8.4f  %8.4f  %8.4f  %6.1f%%  %6.1f%%", N, fm, fw, fa, ew, ea))
        push!(mc_vals, fm)
        push!(wer_vals, fw)
        push!(app_vals, fa)
    end

    ns = collect(N_range)
    plt = plot(ns, mc_vals, marker=:circle, ms=5, label="Monte Carlo",
               xlabel="N (repeaters)", ylabel="Fidelity",
               title="Fidelity scaling (p_s=$ps_fixed, p_w=$pw_fixed)", legend=:topright)
    plot!(plt, ns, wer_vals, ls=:dash, lw=2, marker=:square, ms=3, label="Werner")
    plot!(plt, ns, app_vals, ls=:dot, lw=2, marker=:diamond, ms=3, label="Plan approx")
    savefig(plt, "analysis_scaling_N.png")
    println("\n  -> analysis_scaling_N.png\n")
end

# ============================================================
#  4. N=1 closed form validation
# ============================================================
function analysis_N1_closed_form()
    println("=" ^60)
    println("  4. N=1 CLOSED FORM VALIDATION")
    println("=" ^60)
    println("\nFor N=1, T = max(Geo(p), Geo(p)).")
    println("E[T] = (3 - 2p) / (p(2-p))  [exact for max of 2 geometrics]\n")

    ps_values = 0.1:0.1:0.9

    println(@sprintf("%-6s  %10s  %10s  %8s", "p_s", "MC T", "Closed T", "Err%"))
    println("-" ^40)

    for ps in ps_values
        _, _, tm, _ = Metrics.monte_carlo(1, M_RUNS; p_success=ps, p_w=0.0)
        # Exact: E[max(X1,X2)] for X_i ~ Geo(p) starting at 1
        t_closed = expected_max_geometric(2, ps)
        err = abs(tm - t_closed) / t_closed * 100
        println(@sprintf("p=%.1f   %10.3f  %10.3f  %7.2f%%", ps, tm, t_closed, err))
    end
    println()
end

# ============================================================
#  5. Emergent effects: fidelity distribution shape
# ============================================================
function analysis_emergent()
    println("=" ^60)
    println("  5. EMERGENT EFFECTS")
    println("=" ^60)

    N = 3
    ps = 0.3
    pw = 0.05
    M = 2000

    fidelities = Float64[]
    times = Float64[]
    for _ in 1:M
        f, t = Metrics.single_run(N; p_success=ps, p_w=pw)
        push!(fidelities, f)
        push!(times, t)
    end

    println(@sprintf("\nN=%d, p_s=%.1f, p_w=%.2f, M=%d runs:", N, ps, pw, M))
    println(@sprintf("  Fidelity:  mean=%.4f  std=%.4f  min=%.4f  max=%.4f",
        mean(fidelities), std(fidelities), minimum(fidelities), maximum(fidelities)))
    println(@sprintf("  Time:      mean=%.2f  std=%.2f  min=%d  max=%d",
        mean(times), std(times), minimum(times), maximum(times)))
    println(@sprintf("  Correlation(F, T) = %.4f  (expected: negative)", cor(fidelities, times)))

    # Fidelity distribution histogram
    p1 = histogram(fidelities, bins=30, xlabel="Fidelity", ylabel="Count",
                    title="Fidelity distribution (N=$N, p_s=$ps, p_w=$pw)",
                    legend=false, fillalpha=0.7)

    # Scatter: fidelity vs time
    p2 = scatter(times, fidelities, xlabel="Distribution time", ylabel="Fidelity",
                 title="Fidelity vs Time (N=$N)", legend=false, ms=2, alpha=0.3)

    # Time distribution histogram
    p3 = histogram(times, bins=30, xlabel="Distribution time", ylabel="Count",
                   title="Time distribution (N=$N, p_s=$ps)",
                   legend=false, fillalpha=0.7)

    plt = plot(p1, p2, p3, layout=(1, 3), size=(1500, 400))
    savefig(plt, "analysis_emergent.png")
    println("\n  -> analysis_emergent.png\n")
end

# ============================================================
#  Main
# ============================================================
function main()
    println("\n" * "=" ^60)
    println("  PHASE 5: NUMERICAL vs THEORETICAL ANALYSIS")
    println("  M = $M_RUNS Monte Carlo runs per data point")
    println("=" ^60 * "\n")

    analysis_distribution_time()
    analysis_fidelity()
    analysis_scaling_N()
    analysis_N1_closed_form()
    analysis_emergent()

    println("=" ^60)
    println("  ANALYSIS COMPLETE")
    println("=" ^60)
end

main()
