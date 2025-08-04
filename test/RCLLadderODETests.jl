@testset "RCLLadderODE" begin
    
    # Test for correct dimension
    for n_cells in [1, 10]
        for m in [1, 2, 3]
            if m > 2
                @test_throws AssertionError(
                    "input and output dimension must be either 1 or 2",
                ) RCLLadderConfig(n_cells, m)
            elseif  m > n_cells
                @test_throws AssertionError(
                    "number of inputs and outputs must be less than or equal to the number of cells",
                ) RCLLadderConfig(n_cells, m)
            else
                config = RCLLadderConfig(n_cells, m)
                J, R, Q, B = construct_system(config)
                @test size(J) == (n_cells * 2, n_cells * 2)
                @test size(R) == (n_cells * 2, n_cells * 2)
                @test size(Q) == (n_cells * 2, n_cells * 2)
                @test size(B) == (n_cells * 2, m)
            end
        end
    end

    # Test if using vectors for constants leads to the same systems
    begin
        n_cells = 10
        io_dim = 2
        R = 1
        C = 2
        L = 3
        config_1 = RCLLadderConfig(n_cells, io_dim, R, C, L)
        system_1 = construct_system(config_1)
        
        R_2 = R * ones(n_cells + 1)
        C_2 = C * ones(n_cells)
        L_2 = L * ones(n_cells)
        config_2 = RCLLadderConfig(n_cells, io_dim, R_2, C_2, L_2)
        system_2 = construct_system(config_2)
        @test system_1 == system_2
    end
end