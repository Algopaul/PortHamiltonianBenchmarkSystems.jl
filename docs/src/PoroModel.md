# Poroelastic Network Model

## Description
This benchmark is the poroelastic network model presented in 
[Altmann2021](#References). The model is derived from Biot's consolidation model 
for poroelastic elasticity. Let ``\Omega \subseteq \mathbb{R}^d`` be a Lipschitz 
domain with ``d~\in~\{2,\,3\}`` and ``\mathbb{T} = [0,T]`` for ``T \in (0,\infty)``.
Consider the system of coupled partial differential equations
```math
\begin{aligned}
 \rho \frac{\partial^2}{\partial t^2} u(t,\xi) - \nabla \sigma(u(t,\xi)) + 
 \nabla (\alpha p(t,\xi)) &= \widehat{f}(t,\xi) \quad \text{in } (0,T] \times 
\Omega, \\ 
 \frac{\partial}{\partial t} \left( \alpha \nabla\cdot u(t,\xi) + 
 \frac{1}{M} p(t,\xi)\right) - \nabla \cdot \left( \frac{\kappa}{\nu}\nabla 
p(t,\xi) \right) &= \widehat{g}(t,\xi) \quad \text{in } (0,T] \times \Omega.
\end{aligned}
```
Here, the displacement field ``u: \mathbb{T} \times \Omega \to \R`` and the 
pressure field ``p: \mathbb{T} \times \Omega \to \R`` are searched solution 
functions. Moreover, the stress-strain constitute relation
```math
\sigma(u(t,\xi)) = 2\mu\varepsilon(u(t,\xi)) + 
\lambda(\nabla \cdot u(t,\xi)) \mathcal{I} \quad \text{with} \quad 
\varepsilon(u(t,\xi)) = \frac{1}{2}\left( \nabla u(t,\xi) + (\nabla u(t,\xi))^\mathsf{T} 
\right)
```
is satisfied, where
- ``\mu`` and ``\lambda`` are the Lamé coefficients,
- ``\mathcal{I}`` denotes the identity tensor,
- ``\alpha`` is the Biot-Willes fluid solid coupling coefficient,
- ``M`` is the Biot modulus,
- ``\kappa`` is the permeability,
- ``\rho`` is the density,
- ``\mu`` is the fluid viscosity,
- ``\widehat{f}: (0,T] \times \Omega \to \R^d`` are the volume-distributed forces,
- ``\widehat{g}: (0,T] \times \Omega \to \R`` is the external injection.


The PDE system is equipped with homogeneous Dirichlet boundary conditions
```math
 u(t,\xi) = 0, \quad p(t,\xi) = 0 \quad \text{on 
} (0,T] \times \partial \Omega
```
as well initial conditions ``p(0,\cdot) = p^0 : \partial \Omega \to \R``, 
``u(0,\cdot) = u^0 : \partial \Omega \to \R^d``, and ``\frac{\partial}{\partial 
t}u(0,\cdot) = \dot{u}^0 : \partial \Omega \to \R^d``. Define the Hilbert spaces
```math
 \mathcal{V} := \left[H_0^1(\Omega)\right]^d,\quad \mathcal{H}_{\mathcal{V}} := 
\left[L^2(\Omega)\right]^d, \quad \mathcal{Q} := H_0^1(\Omega),\quad 
\mathcal{H}_{\mathcal{Q}} := L^2(\Omega)
```
and the operators
```math
\begin{aligned}
\mathcal{Y}: \mathcal{H}_{\mathcal{V}} \to \mathcal{H}_{\mathcal{V}}^*,& \quad \left\langle \mathcal{Y}u,v \right\rangle := \int_\Omega \rho u 
v\,\mathrm{d}\xi, \\
\mathcal{M}: \mathcal{H}_{\mathcal{Q}} \to \mathcal{H}_{\mathcal{Q}}^*,& 
\quad \left\langle \mathcal{M}p,q \right\rangle := \int_\Omega \frac{1}{M} 
pq\,\mathrm{d}\xi, \\
\mathcal{A}: \mathcal{V} \to \mathcal{V}^*,& 
\quad \left\langle \mathcal{A}u,v \right\rangle := \int_\Omega 
\sigma(u): \varepsilon(v)\,\mathrm{d}\xi, \\
\mathcal{K}: \mathcal{Q} \to \mathcal{Q}^*,& 
\quad \left\langle \mathcal{K}p,q \right\rangle := \int_\Omega 
\frac{\kappa}{\nu} \nabla p \cdot \nabla q\,\mathrm{d}\xi, \\
\mathcal{D}: \mathcal{V} \to \mathcal{H}_\mathcal{Q}^*,& 
\quad \left\langle \mathcal{D}u,q \right\rangle := \int_\Omega \alpha(\nabla 
\cdot u)q \,\mathrm{d}\xi.
\end{aligned}
```

Note that ``\mathcal{Y}``, ``\mathcal{M}``, ``\mathcal{A}``, and ``\mathcal{K}`` are 
positive definite.
To determine the weak form of the PDE, the first equation is 
multiplied by a test function ``v \in \mathcal{V}`` while the second equation is 
multiplied with ``q \in \mathcal{Q}``. Further we introduce the linear forms
```math
 f(t) := \int_\Omega \widehat{f}(t) \cdot \,\mathrm{d} \xi \in 
\mathcal{H}_{\mathcal{V}}^*, \quad 
 g(t) := \int_\Omega \widehat{g}(t) \cdot \,\mathrm{d} \xi \in 
\mathcal{H}_{\mathcal{Q}}^*.
```
Then for initial conditions ``p^0 \in \mathcal{H}_\mathcal{Q}``, ``u^0 \in 
\mathcal{V}``, and ``\dot{u}^0 \in \mathcal{H}_{\mathcal{V}}`` and right-hand 
sides ``f \in L^2(\mathbb{T},\mathcal{H}_\mathcal{V})`` and ``f \in 
L^2(\mathbb{T},\mathcal{H}_\mathcal{Q})`` one aims to find ``u \in 
L^2(\mathbb{T},\mathcal{V})`` and ``p \in 
L^2(\mathbb{T},\mathcal{Q})`` with ``\dot{u} \in 
L^2(\mathbb{T},\mathcal{H}_\mathcal{V})``, ``\ddot{u} \in 
L^2(\mathbb{T},\mathcal{V}^*)``, and ``\dot{p} \in 
L^2(\mathbb{T},\mathcal{Q}^*)`` such that
```math
\begin{aligned}
 \mathcal{Y} \ddot{u}(t) + \mathcal{A} \dot{u}(t) - \mathcal{D}^* u(t) &= f(t) 
\quad \text{in }  \mathcal{V}^*, \\
\mathcal{D} \dot{u}(t) + \mathcal{M} \dot{p}(t) + \mathcal{K} p(t) &= g(t) 
\quad \text{in } \mathcal{Q}^*
\end{aligned}
```
for almost all ``t \in (0,T)``, where ``\mathcal{D}^*`` denotes the dual operator 
of ``\mathcal{D}``. By introducing the auxiliary variable ``w := \dot{u}``, this 
operator equation can be written in first-order form as
```math
 \begin{bmatrix}
 \mathcal{Y} & 0 & 0 \\ 0 & \mathcal{A} & 0 \\ 0 & 0 & \mathcal{M}
 \end{bmatrix} \begin{pmatrix} \dot{w}(t) \\ \dot{u}(t) \\ \dot{p}(t) 
\end{pmatrix} = 
 \begin{bmatrix}
 0 & -\mathcal{A} & \mathcal{D}^* \\ \mathcal{A}^* & 0 & 0 \\ -\mathcal{D} & 0 
& -\mathcal{K}
 \end{bmatrix} \begin{pmatrix} w(t) \\ u(t) \\ p(t) \end{pmatrix} + 
\begin{pmatrix} f(t) \\ 0 \\ g(t) \end{pmatrix}.
```
Writing the inhomogeneity as
```math
 \begin{pmatrix} f(t) \\ 0 \\ g(t) \end{pmatrix} = \begin{bmatrix} 
\operatorname{id} & 0 \\ 0 & 0 \\ 0 & \operatorname{id} \end{bmatrix} 
\begin{pmatrix} f(t) \\ g(t) \end{pmatrix}
```
and defining the output
```math
 \mathbf{y}(t) := \begin{pmatrix} w(t) \\ p(t) \end{pmatrix} = \begin{bmatrix} 
\operatorname{id} & 0 & 0 \\ 0 & 0 & \operatorname{id} \end{bmatrix} 
\begin{pmatrix} w(t) \\ u(t) \\ p(t) \end{pmatrix},
```
we obtain the port-Hamiltonian system
```math
\begin{aligned}
 \mathcal{E} \dot{\mathbf{x}}(t) &= (\mathcal{J} - \mathcal{R}) \mathbf{x}(t) + 
\mathcal{B}   \mathbf{v}(t), \\
 \mathbf{y}(t) &= \mathcal{B}^* \mathbf{x}(t)
\end{aligned}
```
with 
```math
 \mathcal{E} := \begin{bmatrix}
 \mathcal{Y} & 0 & 0 \\ 0 & \mathcal{A} & 0 \\ 0 & 0 & \mathcal{M}
 \end{bmatrix}, \quad \mathcal{J} := \begin{bmatrix}
 0 & -\mathcal{A} & \mathcal{D}^* \\ \mathcal{A}^* & 0 & 0 \\ -\mathcal{D} & 0 
& 0
 \end{bmatrix}, \quad \mathcal{R} := \begin{bmatrix}
 0 & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 0 
& \mathcal{K}
 \end{bmatrix}, \quad \mathcal{B} := \begin{bmatrix} 
\operatorname{id} & 0 \\ 0 & 0 \\ 0 & \operatorname{id} \end{bmatrix},
```
where ``\mathbf{x}(t) := \left[\begin{smallmatrix} w(t) \\ u(t) \\ p(t) 
\end{smallmatrix}\right]`` and ``\mathbf{v}(t) := \left[\begin{smallmatrix} f(t) 
\\ g(t) \end{smallmatrix}\right]``. 

Discretizing this system with standard ``\mathcal{P}_1`` Lagrange finite elements 
results in the finite-dimensional port-Hamiltonian system 
```math
\begin{aligned}
 E \dot{x}(t) &= (J - R) x(t) + Bv(t), \\
         y(t) &= B^\mathsf{T} x(t)
\end{aligned}
```
with 
```math
\begin{aligned}
 E &:= \begin{bmatrix}
 \rho M_u & 0 & 0 \\ 0 & K_u(\mu,\lambda) & 0 \\ 0 & 0 & \frac{1}{M} M_p
 \end{bmatrix}, \quad J := \begin{bmatrix}
 0 & -K_u(\mu,\lambda) & \alpha D^\mathsf{T} \\ K_u(\mu,\lambda)^\mathsf{T} & 0 & 0 \\ 
-\alpha D & 0 & 0
 \end{bmatrix}, \\ R &:= \begin{bmatrix}
 0 & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 0 
& \frac{\kappa}{\nu} K_p
 \end{bmatrix}, \quad \mathcal{B} := \begin{bmatrix} B_f & 0 \\ 0 & 0 \\ 0 & 
B_g \end{bmatrix}.
\end{aligned}
```

## Parameters
For this benchmark, the domain ``\Omega = [0,1]^2`` with ``d=2`` has been chosen 
and the volume-distributed forces ``\widehat{f}`` and injection ``\widehat{g}`` are 
spatially independent resulting in two inputs, i.\,e., ``B \in \R^{n \times m}`` 
with ``m = 2``. Moreover, different discretization levels are available,  
resulting in systems with state-space dimensions ``n=320``, ``n = 980``, and ``n = 
1805``. These discretizations have been obtained using the `python` 
interface of `FEniCS`. The following fixed parameters have been chosen:
- ``\lambda = 12``,
- ``\mu = 6``.
The following parameters are variable with the default values
- ``\rho = 10^{-3}``,
- ``\alpha = 0.79``,
- ``\frac{1}{M} = 7.80\cdot 10^3``,
- ``\frac{\kappa}{\nu} = 633.33``.

## Interface

The system matrices ``E, J, R,`` and ``B`` can be generated by the following function call.
```julia
using PortHamiltonianBenchmarkSystems
E, J, R, B = poro_elasticity_model()
```

The free parameters are given as named arguments. Note that ``n \in \{ 320, 980, 1805 \}``.
```julia
using PortHamiltonianBenchmarkSystems
E, J, R, B = poro_elasticity_model(n = 320, eta = 1e-3)
H(s) = B'*((s*E-(J-R))\B)
```
Here `H` is the transfer function.

```@docs
poro_elasticity_model
```

## References

```latex
@misc{Altmann2021,
      title={Port-{H}amiltonian formulations of poroelastic network models}, 
      author={R. Altmann and V. Mehrmann and B. Unger},
      year={2021},
      eprint={2012.01949},
      archivePrefix={arXiv},
      primaryClass={math.DS}
}
```
