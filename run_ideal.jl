# Noiseless sanity check: swapping N repeaters over perfect links must give F = 1.

include("src/network.jl")
include("src/swapping.jl")
include("src/metrics.jl")

using .Network, .Swapping, .Metrics

function main()
    println("Noiseless swap, expecting F=1:")
    for N in [1, 2, 3, 5]
        result = Metrics.single_run(N; p_success=1.0, p_w=0.0, ideal=true)
        status = result.fidelity ≈ 1.0 ? "OK" : "FAIL"
        println("  N=$N  F=$(result.fidelity)  T=$(result.dist_time)  $status")
    end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
