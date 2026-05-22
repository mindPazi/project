using QuantumSavory
using Graphs

println("=== Exercise 3: RegisterNet ===\n")

regs = [
    Register([Qubit()], [QuantumOpticsRepr()]),
    Register([Qubit(), Qubit()], [QuantumOpticsRepr(), QuantumOpticsRepr()]),
    Register([Qubit()], [QuantumOpticsRepr()])
]
net = RegisterNet(grid([3]), regs)

println("Nodes: $(nv(net.graph))  (expected: 3)  $(nv(net.graph) == 3 ? "✓" : "✗")")
println("Edges: $(ne(net.graph))  (expected: 2)  $(ne(net.graph) == 2 ? "✓" : "✗")")

println("\nSlots per node:")
for i in 1:3
    n = nsubsystems(net[i])
    label = i == 1 ? "Alice" : i == 3 ? "Bob" : "R1"
    expected = i == 2 ? 2 : 1
    println("  net[$i] ($label): $n slots  (expected: $expected)  $(n == expected ? "✓" : "✗")")
end

bell = (Z1⊗Z1 + Z2⊗Z2) / sqrt(2)
initialize!((net[1][1], net[2][1]), bell)
zz = real(observable((net[1][1], net[2][1]), Z⊗Z))
println("\nBell pair via RegisterNet: ⟨ZZ⟩ = $zz  (expected: 1.0)  $(zz ≈ 1.0 ? "✓" : "✗")")
