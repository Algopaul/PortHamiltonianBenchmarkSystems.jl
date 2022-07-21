using BlockArrays
using IterTools
using LinearAlgebra
using SparseArrays

export DampedWaveNetConfig

"""
Composite type descibing port-Hamiltonian, pressure wave conducting pipe systems, as
described in Egger et al. 'Structure-Preserving Model Reduction for Damped Wave Propagation in Transport Networks'.
# Arguments
- `incidence_matrix`: Sparse incidence matrix describing the pipe network
- `edge_parameters`: Named tuple containing vectors `a`, `b`, `d`, `l`, `n`, respectively containing the parameters ``a_e,\\ b_e,\\ d_e``,
                     the length and the number of FEM elements for each pipe (ordered as in `incidence_matrix`)
- `boundary_conditions`: Vector of chars `'p'`, `'m'`, determining the boundary condition type at each boundary
                         vertex (ordered as in `incidence_matrix`)
"""
struct DampedWaveNetConfig <: BenchmarkConfig
    incidence_matrix::SparseMatrixCSC{Int8,Int64}
    edge_parameters::NamedTuple{
        (:a, :b, :d, :l, :n),
        Tuple{
            Vector{Float64},
            Vector{Float64},
            Vector{Float64},
            Vector{Float64},
            Vector{Int64},
        },
    }
    boundary_conditions::Vector{Char}

    function DampedWaveNetConfig(imat, epar, bcon)
        imat, epar, bcon = convert(Tuple{fieldtypes(DampedWaveNetConfig)...}, (imat, epar, bcon))

        @assert all(in.(nonzeros(imat), Ref([-1, 1]))) "Invalid incidence matrix: found value(s) other than {-1,0,1}"
        @assert all(nnz.(eachcol(imat)) .== 2) "Invalid incidence matrix: found column(s) with other than 2 entries"
        @assert all(sum(abs.(imat), dims = 2) .> 0) "Invalid incidence matrix: found disconnected vertices"
        @assert sum(sum(abs.(imat), dims = 2) .== 1) > 0 "Invalid incidence matrix: found no boundary vertices"
        @assert all(length.(values(epar)) .== size(imat)[2]) "Invalid edge parameters: need same number as edges"
        @assert all([all(v .> 0) for v in values(epar)]) "Invalid edge parameters: all parameters must be positive"
        @assert length(bcon) == sum(sum(abs.(imat), dims = 2) .== 1) "Invalid boundary conditions: need same number as boundary nodes"
        @assert all(in.(bcon, Ref(['p', 'm']))) "Invalid boundary conditions: found identifier other than {p,m}"

        return new(imat, epar, bcon)
    end
end

"""
External constructor providing various default instances of DampedWaveNetConfig.
# Arguments
- `id`: String to identify a default configuration, with possible values: `"pipe"`, `"fork"`, `"diamond"`
"""
function DampedWaveNetConfig(id::String)
    if id == "pipe"
        imat = reshape([-1; 1], :, 1)
        epar = (a = [1], b = [1], d = [1], l = [1], n = [10])
        bcon = ['p', 'm']
    elseif id == "fork"
        imat = [
            -1 0 0
            0 1 0
            0 0 1
            1 -1 -1
        ]
        epar =
            (a = [1, 1, 1], b = [1, 1, 1], d = [1, 1, 1], l = [2, 1, 10], n = [40, 30, 90])
        bcon = ['p', 'p', 'm']
    elseif id == "diamond"
        imat = [
            -1 0 0 0 0 0 0
            0 0 0 0 0 0 1
            1 -1 -1 0 0 0 0
            0 1 0 -1 -1 0 0
            0 0 1 1 0 -1 0
            0 0 0 0 1 1 -1
        ]
        epar = (
            a = [4, 4, 1, 1, 1, 4, 4],
            b = [1, 1, 4, 4, 4, 1, 1] ./ 4,
            d = [1, 1, 8, 8, 8, 1, 1] ./ 80,
            l = [1, 1, 1, 1, 1, 1, 1],
            n = [1, 1, 1, 1, 1, 1, 1] .* 500,
        )
        bcon = ['p', 'm']
    else
        throw("Config id \'" * id * "\' not recognized")
    end
    return DampedWaveNetConfig(imat, epar, bcon)
end

"""
Method for constructing the 'natural' DAE system.
# Arguments
- `config`: `DampedWaveNetConfig` instance
# Output
- `system`: Named tuple containing sparse matrices `E`, `A`, `B`
"""
function construct_system(config::DampedWaveNetConfig)
    #Convenience
    imat = sparse(config.incidence_matrix')
    epar = config.edge_parameters
    bcon = config.boundary_conditions

    #Index calculations
    n_p = sum(epar.n)             #Number of pressure variables
    n_m = sum(epar.n .+ 1)        #Number of mass flow variables
    n_b = length(bcon)            #Number of boundary conditions
    n_bp = sum(bcon .== 'p')       #Number of boundary conditions for p
    n_bm = sum(bcon .== 'm')       #Number of boundary conditions for m
    n_im = size(imat)[2] - n_b     #Number of internal conditions for m
    n_x = n_p + n_m + n_im + n_bm #Number of state variables

    i_ep = collect(eachblock(BlockArray(1:n_p, epar.n))) #Edge indices for p
    i_em = collect(eachblock(BlockArray(1:n_m, epar.n .+ 1))) #Edge indices for m

    #Global FEM system
    E = spzeros(n_x, n_x)
    A = spzeros(n_x, n_x)
    B = spzeros(n_x, n_b)

    Es = PseudoBlockArray(E, [n_p, n_m, n_im + n_bm], [n_p, n_m, n_im + n_bm])
    As = PseudoBlockArray(A, [n_p, n_m, n_im, n_bm], [n_p, n_m, n_im, n_bm])
    Bs = PseudoBlockArray(B, [n_p, n_m, n_im, n_bm], [n_bp, n_bm])

    E11 = view(Es, Block(1, 1))
    E22 = view(Es, Block(2, 2))
    A12 = view(As, Block(1, 2))
    A21 = view(As, Block(2, 1))
    A22 = view(As, Block(2, 2))
    A23 = view(As, Block(2, 3))
    A24 = view(As, Block(2, 4))
    A32 = view(As, Block(3, 2))
    A42 = view(As, Block(4, 2))
    B21 = view(Bs, Block(2, 1))
    B42 = view(Bs, Block(4, 2))

    #Conditions
    i_im, i_bp, i_bm = (1, 1, 1) #Condition counters

    for v in eachcol(imat)
        #Indices & directions of edges connected to v, m indices corresponding to v
        es, ns = (rowvals(v), nonzeros(v))
        i_vm = ifelse.(ns .< 0, first.(i_em[es]), last.(i_em[es]))

        if length(es) > 1 #Internal conditions for m
            A32[i_im, i_vm] = ns
            i_im += 1
        elseif bcon[i_bp + i_bm - 1] == 'p' #Boundary conditions for p
            B21[i_vm, i_bp] = -ns
            i_bp += 1
        else #Boundary conditions for m
            A42[i_bm, i_vm] = ns
            B42[i_bm, i_bm] = 1
            i_bm += 1
        end
    end

    #Physics
    for (e, h) in enumerate(epar.l ./ epar.n)
        #Local matrices
        Mp_loc = [1] * h
        Mm_loc = [2 1; 1 2] * h / 6
        Gm_loc = [-1 1]

        #Contributions for each element
        for (i_p, i_m) in zip(eachrow(i_ep[e]), collect.(partition(i_em[e], 2, 1)))
            E11[i_p, i_p] += Mp_loc .* epar.a[e]
            E22[i_m, i_m] += Mm_loc .* epar.b[e]
            A12[i_p, i_m] -= Gm_loc
            A22[i_m, i_m] -= Mm_loc .* epar.d[e]
        end
    end

    A21 .= -A12'
    A23 .= -A32'
    A24 .= -A42'

    #Packaging
    return (E = E, A = A, B = B)
end

function PHSystem(config::DampedWaveNetConfig)
    E, A, B = construct_system(config)

    J = (A - A') ./ 2
    R = -(A + A') ./ 2
    Q = sparse(1.0I, size(A)...)
    G = B
    P = spzeros(size(B)...)
    S = spzeros(size(B)[2], size(B)[2])
    N = spzeros(size(S)...)

    return PHSystem(E, J, R, Q, G, P, S, N)
end
