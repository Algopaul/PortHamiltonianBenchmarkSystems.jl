@testset "DampedWaveNet" begin
    # Test transfer functions, sizes implicitly
    ids = ["pipe","fork","diamond"]
    Hs = [
        [ 1.8289-0.6127im -0.9100-1.5194im
          0.9100+1.5194im  1.6749-0.9510im],

        [ 1.0329+0.0537im  0.1914+0.1994im  0.0050-0.0004im
          0.1914+0.1994im  0.7863+0.2329im  0.0029-0.0052im
         -0.0050+0.0004im -0.0029+0.0052im  0.9922-0.0956im],

        [ 0.5666+3.0991im -0.8491+0.0570im
          0.8491-0.0570im  0.0463-0.0826im]
    ]

    Hp1s = [Set(Set.(eachrow(H))) for H in Hs]
    Hp2s = [Set(Set.(eachcol(H))) for H in Hs]

    for (id,Hp1,Hp2) in zip(ids,Hp1s,Hp2s)
        s = PHSystem(PortHamiltonianBenchmarkSystems.DampedWaveNet(id))
        H = (s.G+s.P)'*s.Q*(5.0*im*s.E-(s.J-s.R)*s.Q)\Array((s.G-s.P)')

        @test Set(Set.(eachrow(H))) == Hp1
        @test Set(Set.(eachcol(H))) == Hp2
    end
end