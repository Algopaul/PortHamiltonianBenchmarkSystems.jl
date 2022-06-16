using SparseArrays
using BlockArrays
using IterTools

struct DampedWaveNet
	incidence_matrix    ::SparseMatrixCSC{Int8,Int64}
	edge_parameters     ::NamedTuple{(:a,:b,:d,:l,:n),
	                           Tuple{Vector{Float64},
	                                 Vector{Float64},
	                                 Vector{Float64},
	                                 Vector{Float64},
	                                 Vector{Int64  }}}
	boundary_conditions ::Vector{Char}

	function DampedWaveNet(inc_mat,ep,bc)
		inc_mat, ep, bc = convert(Tuple{fieldtypes(DampedWaveNet)...},(inc_mat,ep,bc))

		@assert all(in.(nonzeros(inc_mat),Ref([-1,1]))) "" *
		"Invalid incidence matrix: found value(s) other than {-1,0,1}"

		@assert all(nnz.(eachcol(inc_mat)) .== 2) "" *
		"Invalid incidence matrix: found column(s) with other than 2 entries"

		@assert all(sum(abs.(inc_mat),dims=2) .> 0) "" *
		"Invalid incidence matrix: found disconnected vertices"

		@assert sum(sum(abs.(inc_mat),dims=2) .== 1) > 0 "" *
		"Invalid incidence matrix: found no boundary vertices"

		@assert all(length.(values(ep)) .== size(inc_mat)[2]) "" *
		"Invalid edge parameters: need same number as edges"

		@assert all(ep.d .>= 0) "" *
		"Invalid edge parameters: found damping coefficient(s) < 0"

		@assert all(ep.l .> 0) "" *
		"Invalid edge parameters: found length(s) <= 0"

		@assert all(ep.n .> 0) "" *
		"Invalid edge parameters: found cell number(s) < 1"

		@assert length(bc) == sum(sum(abs.(inc_mat),dims=2) .== 1) "" *
		"Invalid boundary conditions: need same number as boundary nodes"

		@assert all(in.(bc,Ref(['p','m']))) "" *
		"Invalid boundary conditions: found identifier other than {p,m}"

		return new(inc_mat,ep,bc)
	end
end

function DampedWaveNet(id::String)
	if      id == "pipe"
		inc_mat = reshape([ 1;-1],:,1)
		ep = (a=[1],b=[1],d=[1],l=[1],n=[10])
		bc = ['p','p']

	elseif id == "fork"
		inc_mat = [-1  0  0;
		       0  1  0;
		       0  0  1;
		       1 -1 -1]
		ep = (a=[1  ,1  ,1 ],
		      b=[1  ,1  ,1 ],
		      d=[1  ,1  ,1 ],
		      l=[2  ,1  ,10],
		      n=[40 ,30 ,90])
		bc = ['p','p','p']

	elseif id == "diamond"
		inc_mat = [-1  0  0  0  0  0  0;
		       0  0  0  0  0  0  1;
		       1 -1 -1  0  0  0  0;
		       0  1  0 -1 -1  0  0;
		       0  0  1  1  0 -1  0;
		       0  0  0  0  1  1 -1]
		ep = (a=[4  ,4  ,1  ,1  ,1  ,4  ,4  ]     ,
		      b=[1/4,1/4,1  ,1  ,1  ,1/4,1/4]     ,
		      d=[1  ,1  ,8  ,8  ,8  ,1  ,1  ]./80 ,
		      l=[1  ,1  ,1  ,1  ,1  ,1  ,1  ]     ,
		      n=[1  ,1  ,1  ,1  ,1  ,1  ,1  ].*500)
		bc = ['p','p']

	else
		throw("Config id \'"*id*"\' not recognized")
	end

	return DampedWaveNet(inc_mat,ep,bc)
end

function build(problem::DampedWaveNet)
	#Convenience
	inc_mat = sparse(problem.incidence_matrix')
	ep = problem.edge_parameters
	bc = problem.boundary_conditions

	#Index calculations
	n_p  = sum(ep.n)                 #Number of pressure variables
	n_m  = sum(ep.n .+ 1)            #Number of mass flow variables
	n_b  = length(bc)                #Number of boundary conditions
	n_ip = sum(nnz.(eachcol(inc_mat)).-1) #Number of internal conditions for p
	n_im = size(im)[2] - n_b         #Number of internal conditions for m
	n_x  = n_p+n_m+n_ip+n_im+n_b     #Number of state variables
	n_u  = n_b                       #Number of input variables
	n_y  = n_b                       #Number of output variables

	i_ep = collect(eachblock(BlockArray(1:n_p,ep.n   ))) #Edge indices for p
	i_em = collect(eachblock(BlockArray(1:n_m,ep.n.+1))) #Edge indices for m

	#Global FEM system
	E   = spzeros(n_x,n_x)
	A   = spzeros(n_x,n_x)
	B   = spzeros(n_x,n_u)

	Ew  = PseudoBlockArray(E,[n_p,n_m,n_ip+n_im+n_b],[n_p,n_m,n_ip+n_im+n_b])
	Aw  = PseudoBlockArray(A,[n_p,n_m,n_ip,n_im,n_b],[n_p,n_m,n_ip,n_im,n_b])
	Bw  = PseudoBlockArray(B,[n_p+n_m+n_ip+n_im,n_b],[n_b])

	Mp  = view(Ew,Block(1,1)) #Mass matrix for p
	Mm  = view(Ew,Block(2,2)) #Mass matrix for m
	Gp  = view(Aw,Block(2,1)) #Negative gradient matrix for p
	Gm  = view(Aw,Block(1,2)) #Negative gradient matrix for m
	Dm  = view(Aw,Block(2,2)) #Damping matrix for m
	Cip = view(Aw,Block(3,1)) #Internal condition matrix for p
	Cim = view(Aw,Block(4,2)) #Internal condition matrix for m
	Cbp = view(Aw,Block(5,1)) #Boundary condition matrix for p
	Cbm = view(Aw,Block(5,2)) #Boundary condition matrix for m
	Bu  = view(Bw,Block(2,1)) #Input matrix

	#Edge physics
	for (e,h) in enumerate(ep.l./ep.n)
		#Local FEM system on edge e
		Mp_loc = [ 1]*h*ep.a[e]
		Mm_loc = [ 2  1; 1  2]*h/6*ep.b[e]
		Gp_loc = [-1; 1]
		Gm_loc = [ 1 -1]
		Dm_loc = [-2 -1;-1 -2]*h/6*ep.d[e]

		#Contributions for each element
		for (i_p,i_m) in zip(eachrow(i_ep[e]),collect.(partition(i_em[e],2,1)))
			Mp[i_p,i_p] += Mp_loc
			Mm[i_m,i_m] += Mm_loc
			Gp[i_m,i_p] += Gp_loc
			Gm[i_p,i_m] += Gm_loc
			Dm[i_m,i_m] += Dm_loc
		end
	end

	#Algebraic conditions
	i_ip, i_im, i_b = (1,1,1) #Counter for each type of condition
                        
	for v in eachcol(inc_mat)
		#Indices, directions of edges connected to v, corresponding variable indices
		es, ds  = (rowvals(v),nonzeros(v))
		js(i_e) = ifelse.(ds .< 0,first.(i_e[es]),last.(i_e[es]))

		if length(es) > 1
			#Internal conditions for p
			for j_pair in collect.(partition(js(i_ep),2,1))
				Cip[i_ip,j_pair] = [1,-1]
				i_ip += 1
			end

			#Internal conditions for m
			Cim[i_im,js(i_em)] .= ds
			i_im += 1
		else
			#Boundary conditions
			Cb, j       = (bc[i_b] == 'p') ? (Cbp,js(i_ep)) : (Cbm,js(i_em))
			Cb[i_b,j  ] = ds
			Bu[i_b,i_b] = 1
			i_b += 1
		end
	end

	#Apply Lagrange multiplier: L .= C'
	view(Aw,Block(1,3)) .= Cip'
	view(Aw,Block(2,4)) .= Cim'
	view(Aw,Block(1,5)) .= Cbp'
	view(Aw,Block(2,5)) .= Cbm'

	#Packaging
	return (shape=(n_x,n_u,n_y),E=E,A=A,B=B)
end
