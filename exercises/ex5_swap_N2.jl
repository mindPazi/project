using QuantumSavory
using QuantumSavory.CircuitZoo: EntanglementSwap
using Graphs

println("=== Exercise 5: Chain N=2 ===\n")

net = RegisterNet(grid([4]), [
    Register([Qubit()], [QuantumOpticsRepr()]),
    Register([Qubit(), Qubit()], [QuantumOpticsRepr(), QuantumOpticsRepr()]),
    Register([Qubit(), Qubit()], [QuantumOpticsRepr(), QuantumOpticsRepr()]),
    Register([Qubit()], [QuantumOpticsRepr()])
])

bell = (Z1⊗Z1 + Z2⊗Z2) / sqrt(2)
initialize!((net[1][1], net[2][1]), bell)
initialize!((net[2][2], net[3][1]), bell)
initialize!((net[3][2], net[4][1]), bell)

swap = EntanglementSwap()

# After swap at R1, Alice's qubit becomes entangled with R2's slot 1
xm1, zm1 = swap(net[2][1], net[1][1], net[2][2], net[3][1])
println("Swap R1: xm=$xm1, zm=$zm1")

# So the second swap uses net[1][1] as remoteL for net[3][1]
xm2, zm2 = swap(net[3][1], net[1][1], net[3][2], net[4][1])
println("Swap R2: xm=$xm2, zm=$zm2")

xx = real(observable((net[1][1], net[4][1]), X⊗X))
zz = real(observable((net[1][1], net[4][1]), Z⊗Z))
yy = real(observable((net[1][1], net[4][1]), Y⊗Y))
F = (1 + xx + zz - yy) / 4

println("\nAfter 2 swaps (Alice ↔ Bob):")
println("  F = $F  (expected: 1.0)  $(F ≈ 1.0 ? "✓" : "✗")")
