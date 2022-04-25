using LinearAlgebra, SparseArrays
"""
Index 2 RCL Ladder Network

Description:
  This demo provides a semi-explicit index-2 port-Hamiltonian DAE system 
  derived from a simple RCL ladder network with buffer capacitor

Input Arguments:
  - ns:           Number of loops in network
  - r, c, l:      Resistances, capacitances and inductances vectors with length ns+1 / ns / ns
  - m:            Number of inputs (1: SISO or 2: MIMO)

Output Arguments:
  - sys_DAE:      Index-2 PH-DAE model of the RCL ladder network

References:   R. W. Freund. Structure-Preserving Model Order Reduction of
              RCL Circuit Equations, 2008.

  Author: Tim Moser
  E-Mail: tim.moser@tum.de
  Date:   2021/12/14

"""
function  setup_DAE2_RCL_LadderNetwork_sparse(ns, r::AbstractVector, c::AbstractVector, l::AbstractVector, m)
  # Transform scalar inputs to vectors
  # Dimensions
  nk = 2*ns+m; # Number of nodes for network analysis 
  nR = ns+1
  nL = ns
  nC = ns
  # Sub-matrices
  R = spdiagm(0 => r)
  L = spdiagm(0 => l)
  C = spdiagm(0 => c)
  # Initialization
  Ar = spzeros(nk,nR)
  Ac = spzeros(nk,nC)
  Al = spzeros(nk,nL)
  # Ar
  ptr = 1
  for i=1:nR-1
    Ar[ptr,i]=1
    Ar[ptr+1,i]=-1
    ptr = ptr+2
  end
  Ar[nk-m+1,nR] = 1
  if m == 2
    Ar[end,end] = -1
  end
  # Ac
  ptr = 1
  for i=1:nC
    Ac[ptr,i] = -1
    ptr = ptr+2
  end
  # Al
  ptr = 2
  for i=1:nL
    Al[ptr,i]=1
    Al[ptr+1,i]=-1
    ptr = ptr+2
  end
  # Av 
  Av = spzeros(nk,1)
  Av[1] = -1.0
  if m == 2
    Av = hcat(Av, zeros(nk,1))
    Av[end,end] = 1
  end
  # Bring to semi-explicit index-2 pH-Form
  A11 = -Ar*(R\Ar')
  E11 = Ac*C*Ac'; 
  E = blockdiag(sparse(E11),sparse(L),spzeros(m, m))
  J = [spzeros(size(A11, 1), size(A11, 2)) -Al -Av;
       Al' spzeros(size(Al',1),size(Al,2)+size(Av,2));
       Av' spzeros(size(Av',1),size(Al,2)+size(Av,2))]
  R = [-A11 spzeros(size(A11,1),size(Al,2)+size(Av,2));
       spzeros(size(Al',1)+size(Av',1),size(A11,2)+size(Al,2)+size(Av,2))]
  Q = sparse(I(size(A11,1)+size(Al',1)+size(Av',1)))
  B = spzeros(size(A11,1)+size(Al',1)+size(Av',1),m)
  B[end-m+1:end, end-m+1:end] = -sparse(I(m))
  return E, J, R, Q, B
end
