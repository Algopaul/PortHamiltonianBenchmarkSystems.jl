module PortHamiltonianBenchmarkSystems

using LinearAlgebra
using SparseArrays

abstract type BenchmarkConfig end


struct PHSystem{TE,TJ,TR,TQ,TG,TP,TS,TN}
    E::TE
    J::TJ
    R::TR
    Q::TQ
    G::TG
    P::TP
    S::TS
    N::TN
    function PHSystem(
        E::TE,
        J::TJ,
        R::TR,
        Q::TQ,
        G::TG,
        P::TP,
        S::TS,
        N::TN,
    ) where {TE,TJ,TR,TQ,TG,TP,TS,TN}
        n = size(E, 1)
        m = size(G, 2)
        @assert size(E) == (n, n)
        @assert size(J) == (n, n)
        @assert size(R) == (n, n)
        @assert size(Q) == (n, n)
        @assert size(G) == (n, m)
        @assert size(P) == (n, m)
        @assert size(S) == (m, m)
        @assert size(N) == (m, m)
        return new{TE,TJ,TR,TQ,TG,TP,TS,TN}(E, J, R, Q, G, P, S, N)
    end
end

include("./SingleMSDChain.jl")
include("./PoroModel.jl")
include("./RCLLadders.jl")
include("./DampedWaveNet.jl")

export PHSystem, construct_system
end
