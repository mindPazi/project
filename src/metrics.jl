module Metrics

using ..Network, ..Swapping
using QuantumSavory
using Statistics

"""
Compute the fidelity of the final state with respect to |Φ+⟩.
F = (1 + ⟨XX⟩ + ⟨ZZ⟩ - ⟨YY⟩) / 4
"""
function compute_fidelity(alice_slot, bob_slot)
    xx = real(observable((alice_slot, bob_slot), X⊗X))
    zz = real(observable((alice_slot, bob_slot), Z⊗Z))
    yy = real(observable((alice_slot, bob_slot), Y⊗Y))
    (1 + xx + zz - yy) / 4
end

"""
Run a single simulation (ideal or noisy).
Returns (fidelity, distribution_time).
"""
function single_run(N::Int; p_success=1.0, p_w=0.0)
    net = Network.create_network(N; depolarization_rate=p_w)

    if p_success >= 1.0
        Network.generate_entanglement_ideal!(net, N)
        dist_time = 0
    else
        dist_time = Network.generate_entanglement_probabilistic!(net, N, p_success)
    end

    Swapping.perform_swapping!(net, N)

    fidelity = compute_fidelity(net[1][1], net[N + 2][1])
    (fidelity, dist_time)
end

"""
Run M Monte Carlo iterations and collect statistics.
Returns (fidelity_mean, fidelity_std, time_mean, time_std).
"""
function monte_carlo(N::Int, M::Int; p_success=1.0, p_w=0.0)
    fidelities = Vector{Float64}(undef, M)
    times = Vector{Float64}(undef, M)

    for i in 1:M
        f, t = single_run(N; p_success=p_success, p_w=p_w)
        fidelities[i] = f
        times[i] = t
    end

    (mean(fidelities), std(fidelities), mean(times), std(times))
end

end # module
