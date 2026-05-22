module Swapping

using QuantumSavory
using QuantumSavory.CircuitZoo: EntanglementSwap

"""
Perform entanglement swapping across the entire chain of N repeaters.
For each repeater k (node k+1): BSM on its 2 local slots,
with Alice (net[1][1]) as remoteL and the next node's slot 1 as remoteR.
After all swaps, Alice's qubit is entangled with Bob's.
"""
function perform_swapping!(net, N::Int)
    swap = EntanglementSwap()
    for k in 1:N
        node = k + 1
        swap(net[node][1], net[1][1], net[node][2], net[node + 1][1])
    end
end

end # module
