using QuantumSavory

println("=== Exercise 1: Register basics ===\n")

reg = Register([Qubit(), Qubit()], [QuantumOpticsRepr(), QuantumOpticsRepr()])

initialize!(reg[1], Z1)
val = real(observable(reg[1], Z))
println("|0⟩ → ⟨Z⟩ = $val  (expected: 1.0)  $(val ≈ 1.0 ? "✓" : "✗")")
traceout!(reg[1])

initialize!(reg[1], Z2)
val = real(observable(reg[1], Z))
println("|1⟩ → ⟨Z⟩ = $val  (expected: -1.0)  $(val ≈ -1.0 ? "✓" : "✗")")
traceout!(reg[1])

initialize!(reg[1], X1)
val_x = real(observable(reg[1], X))
val_z = real(observable(reg[1], Z))
println("|+⟩ → ⟨X⟩ = $val_x  (expected: 1.0)  $(val_x ≈ 1.0 ? "✓" : "✗")")
println("      ⟨Z⟩ = $val_z  (expected: 0.0)  $(abs(val_z) < 1e-10 ? "✓" : "✗")")
traceout!(reg[1])

initialize!(reg[1], X2)
val = real(observable(reg[1], X))
println("|-⟩ → ⟨X⟩ = $val  (expected: -1.0)  $(val ≈ -1.0 ? "✓" : "✗")")
