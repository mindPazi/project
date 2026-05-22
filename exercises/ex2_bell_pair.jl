using QuantumSavory

println("=== Exercise 2: Bell pair ===\n")

reg_a = Register([Qubit()], [QuantumOpticsRepr()])
reg_b = Register([Qubit()], [QuantumOpticsRepr()])

bell = (Z1⊗Z1 + Z2⊗Z2) / sqrt(2)
initialize!((reg_a[1], reg_b[1]), bell)

xx = real(observable((reg_a[1], reg_b[1]), X⊗X))
zz = real(observable((reg_a[1], reg_b[1]), Z⊗Z))
yy = real(observable((reg_a[1], reg_b[1]), Y⊗Y))

println("⟨XX⟩ = $xx  (expected: 1.0)   $(xx ≈ 1.0 ? "✓" : "✗")")
println("⟨ZZ⟩ = $zz  (expected: 1.0)   $(zz ≈ 1.0 ? "✓" : "✗")")
println("⟨YY⟩ = $yy  (expected: -1.0)  $(yy ≈ -1.0 ? "✓" : "✗")")

F = real((1 + xx + zz - yy) / 4)
println("\nF = (1 + ⟨XX⟩ + ⟨ZZ⟩ - ⟨YY⟩)/4 = $F  (expected: 1.0)  $(F ≈ 1.0 ? "✓" : "✗")")
