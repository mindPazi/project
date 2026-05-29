# Monte Carlo of the noisy chain.
# Sweeps p_success and p_w, then writes the figures under figures/simulation/.

include("src/network.jl")
include("src/swapping.jl")
include("src/metrics.jl")
include("src/plots.jl")

using .Network, .Swapping, .Metrics, .PlotGeneration

# Sweep parameters
const P_SUCCESS_RANGE = 0.1:0.1:1.0
const PW_VALUES = [0.01, 0.05, 0.10]
const N_VALUES = [1, 3, 5]
const M_RUNS = 500

function main()
    # Distribution time vs p_success, one curve per chain length N
    println("Distribution time vs p_success")
    results_time = Dict{Tuple{Int,Float64}, NTuple{4,Float64}}()
    for N in N_VALUES
        for ps in P_SUCCESS_RANGE
            fm, fs, tm, ts = Metrics.monte_carlo(N, M_RUNS; p_success=ps, p_w=0.0)
            results_time[(N, ps)] = (fm, fs, tm, ts)
            println("  N=$N  p_success=$ps  T=$(round(tm, digits=2)) ± $(round(ts, digits=2))")
        end
    end
    PlotGeneration.plot_time_vs_psuccess(results_time, P_SUCCESS_RANGE; N_values=N_VALUES)
    println("  saved figures/simulation/plot_time_vs_psuccess.png\n")

    # Fidelity vs p_success at fixed N, one curve per depolarization rate p_w
    println("Fidelity vs p_success (N=3)")
    N_fixed = 3
    results_fidelity = Dict{Tuple{Float64,Float64}, NTuple{4,Float64}}()
    for pw in PW_VALUES
        for ps in P_SUCCESS_RANGE
            fm, fs, tm, ts = Metrics.monte_carlo(N_fixed, M_RUNS; p_success=ps, p_w=pw)
            results_fidelity[(pw, ps)] = (fm, fs, tm, ts)
            println("  p_w=$pw  p_success=$ps  F=$(round(fm, digits=4)) ± $(round(fs, digits=4))")
        end
    end
    PlotGeneration.plot_fidelity_vs_psuccess(results_fidelity, P_SUCCESS_RANGE; N=N_fixed, pw_values=PW_VALUES)
    println("  saved figures/simulation/plot_fidelity_vs_psuccess.png\n")

    # Fidelity vs chain length N at fixed p_success and p_w
    println("Fidelity vs N (p_success=0.5, p_w=0.05)")
    N_RANGE = 1:7
    results_vs_N = Dict{Int, NTuple{4,Float64}}()
    for n in N_RANGE
        fm, fs, tm, ts = Metrics.monte_carlo(n, M_RUNS; p_success=0.5, p_w=0.05)
        results_vs_N[n] = (fm, fs, tm, ts)
        println("  N=$n  F=$(round(fm, digits=4)) ± $(round(fs, digits=4))")
    end
    PlotGeneration.plot_fidelity_vs_N(results_vs_N, N_RANGE; p_success=0.5, p_w=0.05)
    println("  saved figures/simulation/plot_fidelity_vs_N.png\n")

    # Fidelity over the (p_success, p_w) grid for the heatmap
    println("Fidelity heatmap (N=3)")
    PW_RANGE = 0.01:0.01:0.15
    results_heat = Dict{Tuple{Float64,Float64}, Float64}()
    for ps in P_SUCCESS_RANGE
        for pw in PW_RANGE
            fm, _, _, _ = Metrics.monte_carlo(3, M_RUNS; p_success=ps, p_w=pw)
            results_heat[(ps, pw)] = fm
        end
        println("  p_success=$ps done")
    end
    PlotGeneration.plot_heatmap_fidelity(results_heat, P_SUCCESS_RANGE, PW_RANGE; N=3)
    println("  saved figures/simulation/plot_heatmap_fidelity.png")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
