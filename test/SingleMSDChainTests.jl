@testset "GugercinMSDChain" begin
    # Test for correct dimension
    for n_cells in [1, 10]
        for m in [1, 2]
            if m > n_cells
                @test_throws AssertionError(
                    "number of inputs and outputs must be less than or equal to the number of cells",
                ) SingleMSDConfig(n_cells, m)
            else
                config = SingleMSDConfig(n_cells, m)
                J, R, Q, B = construct_system(config)
                @test size(J) == (n_cells * 2, n_cells * 2)
                @test size(R) == (n_cells * 2, n_cells * 2)
                @test size(Q) == (n_cells * 2, n_cells * 2)
                @test size(B) == (n_cells * 2, m)
            end
        end
    end

    # Test for correct values of transfer function
    begin
        config = SingleMSDConfig(50, 2, 1.0, 4.0, 4.0)
        J, R, Q, B = construct_system(config)
        H(s) = B' * Q * ((s * I - (J - R) * Q) \ B)
        # Values taken from baseline-implementation (not publicly available)
        H1 = [
            0.125712+0.0732172im 0.122624+0.0506278im
            0.122624+0.0506278im 0.120132+0.0531872im
        ]
        s1 = 0.1im
        H2 = [
            0.21075-0.0901108im 0.0225277-0.197312im
            0.0225277-0.197312im 0.0493281+0.00563193im
        ]
        s2 = 1.0im
        @test norm(H(s1) - H1) < 1e-6
        @test norm(H(s2) - H2) < 1e-6
    end

    # Test if using vectors for constants leads to the same systems
    begin
        n_cells_1 = 10
        io_dim_1 = 2
        k_1 = 1.5
        c_1 = 2.5
        m_1 = 3.5
        config_1 = SingleMSDConfig(n_cells_1, io_dim_1, c_1, m_1, k_1)
        system_1 = construct_system(config_1)
        n_cells_2 = 10
        io_dim_2 = 2
        k_2 = 1.5 * ones(2n_cells_2)
        c_2 = 2.5 * ones(2n_cells_2)
        m_2 = 3.5 * ones(2n_cells_2)
        config_2 = SingleMSDConfig(n_cells_2, io_dim_2, c_2, m_2, k_2)
        system_2 = construct_system(config_2)
        @test system_1 == system_2
    end

    # Test direct construction of PHSystem object
    begin
        config = SingleMSDConfig()
        J, R, Q, B = construct_system(config)
        system = PHSystem(config)
        @test system.J == J
        @test system.R == R
        @test system.Q == Q
        @test system.G == B
        @test system.P == zero(B)
        @test norm(system.E - one(R)) == 0
        @test norm(system.S) == 0
        @test norm(system.N) == 0
    end

    # Test setup of control plant
    begin
        for n_cells in [5, 10]
            A, B, C, D, nz, nw = generate_MSD_plant(n_cells)
            @test size(A) == (n_cells * 2, n_cells * 2)
            @test size(B) == (n_cells * 2, 4)
            @test size(C) == (4, n_cells * 2)
            @test size(D) == (4, 4)
            @test nz == 2
            @test nw == 2
        end
    end
end
