# Planar Elasticity (Work in progress)

## Description
This benchmark is the planar elastodynamic problem presented in
[brugnoli2020](#References) and discretized using the Arnold-Falk-Whinter weakly sysmmetric element [arnold2014]{#References}. It corresponds to a mixed Hamiltonian formulation of linearized elasticity. Let ``\Omega \subseteq \mathbb{R}^2`` be a Lipschitz domain and ``\mathbb{T} = [0,T]`` for ``T \in (0,\infty)``.
Consider the system of coupled partial differential equations
```math
\begin{aligned}
 \rho \frac{\partial}{\partial t} \bm{v}(t,\bm{\xi}) &= \mathrm{Div} \bm{\Sigma} + \bm{f}(t, \bm{\xi}) \quad \text{in } (0,T] \times \Omega, \\ 
 \mathcal{C}\frac{\partial}{\partial t} \bm{\Sigma} &= \mathrm{sym}\nabla \bm{v} \quad \text{in } (0,T] \times \Omega.
\end{aligned}
```
The operators appearing in this formulation are the symmetric gradient ``\mathrm{sym}\nabla = \frac{1}{2}(\nabla + \nabla^\top)`` and the row-wise tensor divergence ``(\mathrm{Div} \bm{\Sigma})_i = \sum_j \partial_j \bm{Sigma}_{ij}``.
The velocity field ``\bm{v}: \mathbb{T} \times \Omega \to \R^2`` and the symmetric
stress tensor field ``\bm{\Sigma}: \mathbb{T} \times \Omega \to \R^{2\times 2}_{\text{sym}}`` are the unknown of the problem. The physical parameters are
- ``\rho`` the density,
- ``\mathcal{C}= \frac{1}{2\mu}\left[(\cdot) - \frac{\lambda}{2\mu + 2\lambda} \mathrm{Tr}(\cdot) \mathbf{I}`` the isotropic compliance tensor, where ``\mu, \; \lambda`` are the Lam\`e coefficients ``\mathbf{I}`` is the identity in ``\R^2``. 
The forcing term is given by a gaussian spatial distribution and a control term depending only on time, i.e.
```math
\bm{f} = g(\bm{\xi})\mathbf{u}(t), \qquad \mathbf{u}(t) \in \R^2.
```
The system can be compactly rewritten as
```math
\mathcal{E}\partial_t \bm{x} = \mathcal{J}\bm{x} + \mathcal{B}\bm{f} \quad \text{in } (0,T] \times \Omega.
```
where ``\mathcal{E}`` is a bounded, symmetric and uniformly positive operator, ``\mathcal{J}`` is a formally skew-adjoint operator and ``\mathcal{B}`` is a bounded operator for the distributed forcing.

An homogenous Dirichlet condition holds at the boundary
```math
 \bm{v}(t,\xi) = \bm{0}, \quad \text{on} (0,T] \times \partial \Omega
```
as well initial conditions ``\bm{v}(0,\cdot) = \bm{v}^0 : \partial \Omega \to \R^2``, ``\bm{\Sigma}(0,\cdot) = \bm{\Sigma}^0 : \partial \Omega \to \R^{d \times d}_{\mathrm{sym}}``. 
The construction of conforming finite elements for symmetric tensors is involved. Therefore, a weakly symmetric finite element formulation will be employed, based on the following decomposition of the symmetric gradient 
```math
 \mathrm{sym}\nabla \bm{v} = \nabla \bm{v} - \mathrm{skw}(\nabla \bm{v})
```
where ``skw \bm{A} = (\bm{A} - \bm{A}^\top)/2 \in \R^{2\times 2}_{\mathrm{skw}}`` is the skew-symmetric part of a matrix. 
To introduce the weak form, define the Hilbert spaces
```math
\begin{aligned}
{V} := L^2(\Omega; \R^2), \\
{M} := H^{\mathrm{div}}(\Omega; \R^{2 \times 2}) = \{\bm{A} \in L^2(\Omega; \R^{2\times 2})|  \mathrm{div}\bm{A} \in L^2(\Omega; \R^{2})\}, \\
{K} := L^2(\Omega; \R^{2\times 2}_{\mathrm{skw}}), 
\end{aligned}
```
To determine the weak form of the PDE, the first equation is 
multiplied by a test function ``\bm{\psi}_v \in {V}``. The second equation is multiplied by ``\bm{\Psi}_{\Sigma} \in M`` and integrated by parts. The homogeneous Dirichlet boundary condition is included naturally in this formulation. A third equation imposes the symmetry weakly by taking the inner product between a skew-symmetric test function and the stress tensor. The weak formulation therefore reads: find ``\bm{v} \in {V}, \; \bm{Sigma} \in {M}, \; \bm{R} \in {K}`` such that
```math
\begin{aligned}
\int_\Omega \bm{\psi}_v \cdot \partial_t \bm{v} dx &= \int_\Omega \bm{\psi}_v \cdot \mathrm{Div} \bm{\Sigma} dx + \int_\Omega \bm{\psi}_v \cdot \bm{f} dx \qquad \forall \bm{\psi}_v \in \mathcal{V}, \\
\int_\Omega \bm{\Psi}_\Sigma : \partial_t \bm{\Sigma} dx + \int_\Omega \bm{\Psi}_\Sigma : \partial_t \bm{R} dx &= -\int_\Omega \mathrm{Div} \bm{\psi}_\Sigma \cdot \bm{v} dx \qquad \forall \bm{\Psi}_{\Sigma} \in \mathcal{M}, \\
\int_{\Omega} \bm{\Psi}_R : \partial_t \bm{\Sigma} &= 0 \qquad \forall \bm{\Psi}_R \in \mathcal{K}
\end{aligned}
```
where ``\bm{A} : \bm{B} = \Sum_{i,j} \bm{A}_{ij} \bm{B}_{ij}`` denotes the tensor contraction. Consider a regular triangulation ``\mathcal{T}_h`` with elements ``T``. At the boundary ``\partial\Omega``. The following conforming finite element spaces are used to discretize the weak formulation 
```math
\begin{aligned}
V_h &= \{\bm{v}_h \in L^2(\Omega; \R^2)|  \; \forall T \in \mathcal{T}_h, \bm{v}_h|_T \in \mathrm{DG}_{k-1}^2\}, \\
M_h &= \{\bm{\Sigma}_h \in H^{\mathrm{div}}(\Omega; \R^{2 \times 2})|  \; \forall T \in \mathcal{T}_h \; \text{ and for } i=1,2, (\bm{\Sigma}_{i1, h}, \; \bm{\Sigma}_{j2, h})|_T \in \mathrm{BDM}_{k}\},
K_h &= \{\bm{R}_h \in L^2(\Omega; \R^{2\times 2}_{\mathrm{skw}})|  \; \forall T \in \mathcal{T}_h, \bm{R}_h|_T \in \mathrm{DG}_{k-1}\}, \\
\end{aligned}
```
where ``k \ge 1`` is the polynomial degree of the finite elements and the acronyms stand for
- DG: the discontinous galerkin finite element space;
- BDM: the Brezzi-Douglas-Marini space;

Inserting the finite element approximation into the weak formulation, the following port-Hamiltonian system is obtained
```math
\begin{aligned}
 \mathbf{E} \dot{\mathbf{x}}(t) &= \mathbf{J}\mathbf{x}(t) + \mathbf{B} \mathbf{u}(t), \\
  \mathrm{y}(t) &= \mathbf{B}^\mathsf{T} \mathbf{x}(t)
\end{aligned}
```
with 
```math
\begin{aligned}
 \mathbf{E} &:= \begin{bmatrix}
 \mathbf{M}_\rho & 0 & 0 \\ 0 & \mathbf{M}_{\mathcal{C}} & \mathbf{A}_{\lambda}^\mathsf{T} \\ 0 & \mathbf{A}_{\lambda} & 0
 \end{bmatrix}, \\
  \mathbf{J} &:= \begin{bmatrix}
 0 & \mathbf{D}_{\mathrm{div}} & 0 \\ - \mathbf{D}_{\mathrm{div}}^\mathsf{T} & 0 & 0 \\  0 & 0 & 0
 \end{bmatrix}, \\ 
 \mathbf{B} := \begin{bmatrix} \mathbf{B}_f \\ 0 \\ 0 \end{bmatrix}.
\end{aligned}
```

## Parameters
For this benchmark, a square domain ``\Omega = [0,1]^2`` is considered. The forcing term enters the system via a gaussian distribution ``g(\xi)`` given by
```math
g(\bm{\xi}) = \exp(-\frac{100}{2}||\bm{\xi}- \frac{1}{2}||^2).
```
Moreover, different discretization levels are available,  
resulting in systems with state-space dimensions ``n=1260`` (5 elements per side and ``k=2``), ``n = 1880`` (10 elements per side and ``k=1``), and ``n = 4920`` (10 elements per side and ``k=2``). These discretizations have been obtained using the `python` 
interface of `FEniCS`. The following fixed parameters have been chosen:
- ``\lambda = 20``,
-  ``\mu = 4``
- ``\rho = 1``,

## Interface

```@docs
Elasticity2DAFWConfig
```
```@docs
Elasticity2DAFWConfig()
```

## References

```latex
@article{Arnold2014,
   author  = {D. Arnold and J. Lee},
   issue   = {6},
   journal = {SIAM Journal on Numerical Analysis},
   pages   = {2743-2769},
   title   = {Mixed Methods for Elastodynamics with Weak Symmetry},
   volume  = {52},
   url     = {https://epubs.siam.org/doi/10.1137/13095032X},
   year    = {2014}
   }

@phdthesis{Brugnoli2020,
  author = {A.~Brugnoli},
  title  = {A port-{H}amiltonian formulation of flexible structures. Modelling and structure-preserving finite element discretization},
  school = {Universit\'e de Toulouse, ISAE-SUPAERO, France},
  year   = {2020}
  }
```
