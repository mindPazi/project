using QuantumSavory
using QuantumSavory.CircuitZoo: EntanglementSwap
using Graphs

println("=== Exercise 6: Depolarization ===\n")

# Depolarization(τ): probability of error after Δt is p = 1 - exp(-Δt/τ)
function run_with_noise(τ, wait_time)
    regs = if τ == Inf
        [
            Register([Qubit()], [QuantumOpticsRepr()]),
            Register([Qubit(), Qubit()], [QuantumOpticsRepr(), QuantumOpticsRepr()]),
            Register([Qubit()], [QuantumOpticsRepr()])
        ]
    else
        [
            Register([Qubit()], [QuantumOpticsRepr()], [Depolarization(τ)]),
            Register([Qubit(), Qubit()], [QuantumOpticsRepr(), QuantumOpticsRepr()], [Depolarization(τ), Depolarization(τ)]),
            Register([Qubit()], [QuantumOpticsRepr()], [Depolarization(τ)])
        ]
    end

    net = RegisterNet(grid([3]), regs)
    bell = (Z1⊗Z1 + Z2⊗Z2) / sqrt(2)

    initialize!((net[1][1], net[2][1]), bell; time=0.0)
    initialize!((net[2][2], net[3][1]), bell; time=0.0)

    uptotime!((net[1][1], net[2][1], net[2][2], net[3][1]), wait_time)

    swap = EntanglementSwap()
    swap(net[2][1], net[1][1], net[2][2], net[3][1])

    xx = real(observable((net[1][1], net[3][1]), X⊗X))
    zz = real(observable((net[1][1], net[3][1]), Z⊗Z))
    yy = real(observable((net[1][1], net[3][1]), Y⊗Y))
    return (1 + xx + zz - yy) / 4
end

F_ideal = run_with_noise(Inf, 50.0)
F_low   = run_with_noise(1000.0, 50.0)
F_med   = run_with_noise(100.0, 50.0)
F_high  = run_with_noise(10.0, 50.0)

println("τ=Inf   wait=50 → F=$(round(F_ideal, digits=4))  $(F_ideal ≈ 1.0 ? "✓" : "✗")")
println("τ=1000  wait=50 → F=$(round(F_low, digits=4))")
println("τ=100   wait=50 → F=$(round(F_med, digits=4))")
println("τ=10    wait=50 → F=$(round(F_high, digits=4))")

ok = F_ideal ≥ F_low ≥ F_med ≥ F_high
println("\nF_ideal ≥ F_low ≥ F_med ≥ F_high → $ok  $(ok ? "✓" : "✗")")

println("\nVarying wait time (τ=100):")
for t in [0.0, 10.0, 50.0, 200.0]
    F = run_with_noise(100.0, t)
    println("  wait=$t → F=$(round(F, digits=4))")
end
