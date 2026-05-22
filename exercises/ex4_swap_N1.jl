using QuantumSavory
using QuantumSavory.CircuitZoo: EntanglementSwap
using Graphs

println("=== Exercise 4: Entanglement swap (N=1) ===\n")

net = RegisterNet(grid([3]), [
    Register([Qubit()], [QuantumOpticsRepr()]),
    Register([Qubit(), Qubit()], [QuantumOpticsRepr(), QuantumOpticsRepr()]),
    Register([Qubit()], [QuantumOpticsRepr()])
])

bell = (Z1⊗Z1 + Z2⊗Z2) / sqrt(2)
initialize!((net[1][1], net[2][1]), bell)
initialize!((net[2][2], net[3][1]), bell)

println("Before swap:")
println("  Alice-R1 ⟨ZZ⟩ = $(real(observable((net[1][1], net[2][1]), Z⊗Z)))")
println("  R1-Bob   ⟨ZZ⟩ = $(real(observable((net[2][2], net[3][1]), Z⊗Z)))")

# EntanglementSwap()(localL, remoteL, localR, remoteR)
swap = EntanglementSwap()
xm, zm = swap(net[2][1], net[1][1], net[2][2], net[3][1])
println("\nBSM outcome: xmeas=$xm, zmeas=$zm")

xx = real(observable((net[1][1], net[3][1]), X⊗X))
zz = real(observable((net[1][1], net[3][1]), Z⊗Z))
yy = real(observable((net[1][1], net[3][1]), Y⊗Y))
F = (1 + xx + zz - yy) / 4

println("\nAfter swap (Alice ↔ Bob):")
println("  ⟨XX⟩=$xx  ⟨ZZ⟩=$zz  ⟨YY⟩=$yy")
println("  F = $F  (expected: 1.0)  $(F ≈ 1.0 ? "✓" : "✗")")
