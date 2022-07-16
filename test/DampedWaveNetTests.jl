@testset "DampedWaveNet" begin
    ids = ["pipe","fork","diamond"]
    ref_shapes = [
      (n_x=22,n_u=2,n_y=2)
      (n_x=325,n_u=3,n_y=3)
      (n_x=7012,n_u=2,n_y=2)
    ]
    ref_Hs = [
        [ 1.8289-0.6127im -0.9100-1.5194im
          0.9100+1.5194im  1.6749-0.9510im],

        [ 1.0329+0.0537im  0.1914+0.1994im  0.0050-0.0004im
          0.1914+0.1994im  0.7863+0.2329im  0.0029-0.0052im
         -0.0050+0.0004im -0.0029+0.0052im  0.9922-0.0956im],

        [ 0.5666+3.0991im -0.8491+0.0570im
          0.8491-0.0570im  0.0463-0.0826im]
    ]

    for (id,ref_H,(n_x,n_u,n_y)) in zip(ids,ref_Hs,ref_shapes)
        p     = DampedWaveNet(id)
        E,A,B = construct_system(p)
        s     = PHSystem(p)
        H     = B' * ((5.0*im*E - A)\Array(B))

        #Correct system shapes
        @test size(E ) == (n_x,n_x)
        @test size(A ) == (n_x,n_x)
        @test size(B ) == (n_x,n_u)
        @test size(B') == (n_y,n_x)

        #Correspondence between natural and pH form
        @test s.E == E
        @test (s.J-s.R)*s.Q == A
        @test s.G-s.P == B
        @test (s.G+s.P)'*s.Q == B'
        @test (s.S+s.N) == spzeros(n_y,n_u)

        #Correct system transfer function
        @test round.(H,digits=4) == ref_H
    end
end