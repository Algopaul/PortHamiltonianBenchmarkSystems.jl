
"""
This struct configures port-Hamiltonian ODE RCL ladder network described in [PS10](@cite) and [GPBS12](@cite).

# Arguments
- `n_cells::Int`: The number of cells in the ladder network
- `io_dim::Int`: The input and output dimension of the system
- `R::Vector{T}`: The resistances (Vector of length `n_cells + 1`)
- `C::Vector{T}`: The capacitances (Vector of length `n_cells`)
- `L::Vector{T}`: The inductances (Vector of length `n_cells`)
# Outputs
- `config`: The configuration struct for the system. The system can subsequently be created with `construct_system(config)`
"""
struct RCLLadderConfig{T} <: BenchmarkConfig
    n_cells::Int
    io_dim::Int
    R::Vector{T}
    C::Vector{T}
    L::Vector{T}
    function RCLLadderConfig(
        n_cells::Int=50,
        io_dim::Int=1,
        R::Vector{T}=[0.2 * ones(n_cells)...; 0.4],
        C::Vector{T}=ones(n_cells),    
        L::Vector{T}=ones(n_cells),
    ) where {T}
    @assert n_cells > 0 "number of cells must be positive"
    @assert io_dim == 1 || io_dim == 2 "input and output dimension must be either 1 or 2"
    @assert io_dim <= n_cells "number of inputs and outputs must be less than or equal to the number of cells"
    @assert length(R) == n_cells + 1 "length of R must be n_cells + 1"
    @assert length(C) == n_cells "length of C must be n_cells"
    @assert length(L) == n_cells "length of L must be n_cells"
    @assert all(R .>= 0) "R cannot have any negative entries"
    @assert all(C .>= 0) "C cannot have any negative entries"
    @assert all(L .>= 0) "L cannot have any negative entries"

    return new{T}(n_cells, io_dim, R, C, L)
    end
end

function RCLLadderConfig(n_cells::Int, io_dim::Int, R::Union{Number, AbstractVector}, C::Union{Number, AbstractVector}, L::Union{Number, AbstractVector})
    T = promote_type(eltype(R), eltype(C), eltype(L))
    
    if isa(R, Number)
        R = fill(convert(T, R), n_cells + 1)
    else
        R = convert(Vector{T}, R)
    end

    if isa(C, Number)
        C = fill(convert(T, C), n_cells)
    else
        C = convert(Vector{T}, C)
    end

    if isa(L, Number)
        L = fill(convert(T, L), n_cells)
    else
        L = convert(Vector{T}, L)
    end

    return RCLLadderConfig(n_cells, io_dim, R, C, L)
end


"""
    RCLLadderConfig(id::String)

External constructor providing various default instances of RCLLadderConfig.
# Arguments
- `id`: The identifier of the desired configuration. Use 
    - `"PS10"` for the setup in [PS10](@cite) or 
    - `"GPBS12"` for the setup in [GPBS12](@cite).
# Outputs
- `config`: Instance of `RCLLadderConfig`.
"""
function RCLLadderConfig(id::String)
    if id == "PS10"
        n_cells = 50
        io_dim = 1
        R = [0.2 * ones(n_cells)...; 0.4]
        C = ones(n_cells)
        L = ones(n_cells)
        return RCLLadderConfig(n_cells, io_dim, R, C, L)
    elseif id == "GPBS12"
        n_cells = 50
        io_dim = 2
        R = [3 * ones(n_cells)...; 1]
        C = 0.1 * ones(n_cells)
        L = 0.1 * ones(n_cells)
        return RCLLadderConfig(n_cells, io_dim, R, C, L)
    end
end

function construct_system(config::RCLLadderConfig)
    n = 2 * config.n_cells
    m = config.io_dim

    # J = diagm(1 => [-1*ones(n-2); 1], -1 => [ones(n-2); -1])
    J = diagm(1 => -1*ones(n-1), -1 => ones(n-1))

    diagR = zeros(n)
    index_diagR = 2:2:n
    diagR[index_diagR] = config.R[1:end-1]
    diagR[end] = diagR[end] + config.R[end]
    R = diagm(diagR)
    
    diagQ = zeros(n)
    index_diagQ_C = 1:2:n-1
    index_diagQ_L = 2:2:n
    diagQ[index_diagQ_C] .= 1 ./ config.C
    diagQ[index_diagQ_L] .= 1 ./ config.L
    Q = diagm(diagQ)

    # G
    G = zeros(n,m); 
    if m == 1
        G[1, 1] = 1;
    elseif m == 2
        G[1, 1] = 1;
        G[end, 2] = 1;
    else
        @error "Only m=1 and m=2 are supported, not m=$m"
    end
    
    return J, R, Q, G    
end

export RCLLadderConfig
