@testset "RCL-Ladders" begin
    # No reference available, thus only testing sizes.
    for ns in [10, 20]
        for m in [1, 2]
            E1, J1, R1, Q1, G1 = setup_DAE1_RCL_LadderNetwork_sparse(ns = ns, m = m)
            E2, J2, R2, G2 = setup_DAE2_RCL_LadderNetwork_sparse(ns = ns, m = m)
            n = 3ns + 2m
            @test size(E1) == (n, n)
            @test size(E2) == (n, n)
            @test size(J1) == (n, n)
            @test size(J2) == (n, n)
            @test size(R1) == (n, n)
            @test size(R2) == (n, n)
            @test size(Q1) == (n, n)
            @test size(G1) == (n, m)
            @test size(G2) == (n, m)
        end
    end
end
