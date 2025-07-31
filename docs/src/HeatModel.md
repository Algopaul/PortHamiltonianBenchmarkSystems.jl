# Heat equation with Neumann boundary control

## Description

This example considers the temperature $T$ of a 2D domain $\Omega \subset \mathbb{R}^2$. 
We follow [SHM19](@cite) for the derivation of the port-Hamiltonian formulation and [SHM19a](@cite) for the subsequent structure-preserving discretization.  
Denoting $C_v$ the heat capacity (at constant volume), $\rho$ the mass density and $\lambda$ the heat conductivity, a positive definite tensor, leads to the following well-known *heat equation*

```math
\rho(x) C_v(x) \frac{\partial}{\partial t} T(t,x) - {\rm div} \left( \lambda(x) \cdot {\rm grad} \left( T(t,x) \right) \right) = 0, \quad t \ge 0, \, x \in \Omega,
```

together with *Neumann boundary control* 

```math
\ - \left( \lambda(x) \cdot {\rm grad} \left( T(t,x) \right) \right) \cdot \mathbf{n} = u_\partial(t,x), \quad t \ge 0, \, x \in \partial \Omega,
```

where $\mathbf{n}$ is the outward normal to $\Omega$.

The **Hamiltonian** is taken as the usual $L^2$ functional, despite its lack of thermodynamical meaning

```math
\mathcal{H}(t) := \frac{1}{2} \int_\Omega \rho(x) C_v(x) \left( T(t,x) \right)^2 {\rm d}x, \qquad t \ge 0.
```

Taking the *internal energy density* $\alpha_u := u = C_v T$ as **energy variable** (with Dulong-Petit model), the Hamiltonian rewrites

```math
\mathcal{H}(t) = \mathcal{H}(\alpha_u(t,\cdot)) = \frac{1}{2} \int_\Omega \rho(x) \frac{\alpha_u^2(t,x)}{C_v(x)} {\rm d}x.
```

The **co-energy variable** is the variational derivatives of the Hamiltonian, with respect to the weighted $L^2$-product with weight $\rho$

```math
e_u := \delta^\rho_{\alpha_u} \mathcal{H} = \frac{\alpha_u}{C_v} = T,
```

*i.e.* the *temperature*. This equality is the first **constitutive relation**.

Denoting $J_Q$ the heat flux, the first principle of thermodynamics reads: $\frac{\partial}{\partial t} \rho u + {\rm div} \left( J_Q \right) = 0$. Fourier's law gives the second constitutive relation: $J_Q = - \lambda \cdot {\rm grad} \left( T \right)$. Then one can write

```math
\begin{bmatrix}
\rho \frac{\partial}{\partial t} \alpha_u \\
f_Q
\end{bmatrix} =
\begin{bmatrix}
0 & -{\rm div} \\
? & 0
\end{bmatrix}
\begin{bmatrix}
e_u \\
J_Q
\end{bmatrix}, 
\qquad \left\lbrace
\begin{array}{rcl}
e_u &=& \frac{\alpha_u}{C_v}, \\
J_Q &=& -\lambda \cdot {\rm grad} \left( T \right),
\end{array}\right.
\qquad \left\lbrace
\begin{array}{rcl}
u_\partial &=& J_Q \cdot \mathbf{n}, \\
y_\partial &=& e_u\mid_{\partial \Omega},
\end{array}\right.
```

As port-Hamiltonian systems deal with *formally skew-symmetric operator*, the example of the wave equation leads us to complete this system with $-{\rm grad}$ in order to obtain the heat equation as a **port-Hamiltonian system**

```math
\begin{bmatrix}
\rho \frac{\partial}{\partial t} \alpha_u \\
f_Q
\end{bmatrix} =
\begin{bmatrix}
0 & -{\rm div} \\
-{\rm grad} & 0
\end{bmatrix}
\begin{bmatrix}
e_u \\
J_Q
\end{bmatrix}, 
\qquad \left\lbrace
\begin{array}{rcl}
e_u &=& \frac{\alpha_u}{C_v}, \\
J_Q &=& \lambda \cdot f_Q,
\end{array}\right.
\qquad \left\lbrace
\begin{array}{rcl}
u_\partial &=& J_Q \cdot \mathbf{n}, \\
y_\partial &=& e_u\mid_{\partial \Omega},
\end{array}\right.
```

The **power balance** satisfied by the Hamiltonian is

```math
\frac{\rm d}{ {\rm d}t} \mathcal{H} = - \int_\Omega J_Q \cdot \lambda^{-1} \cdot J_Q + \langle u_\partial, y_\partial \rangle_{H^{-\frac12}(\partial \Omega),H^\frac12(\partial \Omega)}.
```

To get rid of the first algebraic constraint induced by the constitutive relation $e_u = \frac{\alpha_u}{C_v}$, one rewrites $\rho C_v \frac{\partial}{\partial t} T$. Furthermore, we also include Fourier's law as $\lambda^{-1} \cdot J_Q = f_Q$ inside the Dirac structure. The port-Hamiltonian system then reads

```math
\begin{bmatrix}
\rho C_v \frac{\partial}{\partial t} T \\
\lambda^{-1} \cdot J_Q
\end{bmatrix} =
\begin{bmatrix}
0 & -{\rm div} \\
-{\rm grad} & 0
\end{bmatrix}
\begin{bmatrix}
T \\
J\_Q
\end{bmatrix},
\qquad \left\lbrace
\begin{array}{rcl}
u_\partial &=& J_Q \cdot \mathbf{n}, \\
y_\partial &=& e_u\mid_{\partial \Omega}.
\end{array}\right.
```

### The Partitioned Finite Element Method



#### Weak formulation

Let $\varphi_T$, $\phi_Q$ and $\psi$ be scalar-valued, vector-valued and boundary scalar-valued test functions respectively. The weak formulation reads

```math
\left\lbrace
\begin{array}{rcl}
\displaystyle \int_\Omega \varphi_T \rho C_v \frac{\partial}{\partial t} T 
&=& \displaystyle - \int_\Omega \varphi_T {\rm div} \left( J_Q \right), \\
\displaystyle \int_\Omega \phi_Q \cdot \lambda^{-1} \cdot J_Q 
&=& \displaystyle - \int_\Omega \phi_Q \cdot {\rm grad} \left( T \right), \\
\displaystyle \int_{\partial \Omega} \psi y_\partial &=& \displaystyle \int_{\partial \Omega} \psi T.
\end{array}\right.
```

#### Integration by parts

The integration by parts of the first line makes $u_\partial = \mathbf{J}\_Q \cdot \mathbf{n}$ appear

```math
\left\lbrace
\begin{array}{rcl}
\displaystyle \int_\Omega \varphi_T \rho C_v \frac{\partial}{\partial t} T 
&=& \displaystyle \int_\Omega {\rm grad} \left( \varphi_T \right) J_Q 
-\int_{\partial \Omega} u_\partial \varphi_T, \\
\displaystyle \int_\Omega \phi_Q \cdot \lambda^{-1} \cdot J_Q 
&=& \displaystyle - \int_\Omega \phi_Q \cdot {\rm grad} \left( T \right), \\
\displaystyle \int_{\partial \Omega} \psi y_\partial &=& \displaystyle \int_{\partial \Omega} \psi T.
\end{array}\right.
```

#### Projection

Let $(\varphi_T^i)\_{1 \le i \le N_T}$, $(\phi_Q^j)\_{1 \le j \le N_Q}$ and $(\psi^k)\_{1 \le k \le N_\partial}$ be finite element families for $u$-type, $Q$-type and boundary-type variables. Variables are approximated in their respective finite element family

```math
T^d(t,x) := \sum_{i=1}^{N_T} T^i(t) \varphi_T^i(x),
\qquad J_Q^d(t,x) := \sum_{j=1}^{N_Q} J_Q^j(t) \phi_Q^j(x),
```

```math
u_\partial^d(t,x) := \sum_{k=1}^{N_\partial} u_\partial^k(t) \psi^k(x),
\qquad y_\partial^d(t,x) := \sum_{k=1}^{N_\partial} y_\partial^k(t) \psi^k(x).
```

Denoting $\underline{\star}$ the (time-varying) vector of coordinates of $\star^d$ in its respective finite element family, the discrete system reads

```math
\underset{M}{\underbrace{\begin{bmatrix}
M_T & 0 & 0 \\
0 & M_Q & 0 \\
0 & 0 & M_\partial
\end{bmatrix} } }
\begin{bmatrix}
\frac{\rm d}{ {\rm d}t} \underline{T}(t) \\
\underline{J_Q}(t) \\
-\underline{y_\partial}(t)
\end{bmatrix} =
\underset{J}{\underbrace{\begin{bmatrix}
0 & D & B \\
-D^\top & 0 & 0 \\
-B^\top & 0 & 0
\end{bmatrix} } }
\begin{bmatrix}
\underline{T}(t) \\
\underline{J_Q}(t) \\
\underline{u_\partial}(t)
\end{bmatrix}
```

where

```math
(M_T)_{ij} := \int_\Omega \varphi_T^i \rho C_v \varphi_T^j,
\qquad 
(M_Q)_{ij} := \int_\Omega \phi_Q^i \cdot \lambda^{-1} \cdot \phi_Q^j,
\qquad 
(M_\partial)_{ij} := \int_{\partial \Omega} \psi^i \psi^j,
```

and

```math
(D)_{ij} := \int_\Omega {\rm grad} \left( \varphi_T^i \right) \cdot \phi_Q^j,
\qquad
(B)_{jk} := \int_{\partial \Omega} \varphi_T^j \psi^k.
```


#### Discrete Hamiltonian

By definition, the discrete Hamiltonian is equal to the continuous Hamiltonian evaluated in the approximated variables. Recalling that

```math
\mathcal{H} = \frac{1}{2} \int_\Omega \rho C_v (T)^2.
```

Then, the discrete Hamiltonian is defined as

```math
\mathcal{H}^d := \frac{1}{2} \int_\Omega \rho C_v (T^d)^2.
```

After straightforward computations, it comes

```math
\mathcal{H}^d(t) = \frac{1}{2} \underline{T}(t)^\top M_T \underline{T}(t),
```

and the **discrete power balance** follows

```math
\frac{\rm d}{ {\rm d}t} \mathcal{H}^d(t) = - \underline{J_Q}(t)^\top M_Q \underline{J_Q}(t) + \underline{u_\partial}(t)^\top M_\partial \underline{y_\partial}(t).
```

## Interface

We have generated models in domains of three different shapes `disc`, `rectangle`, and `L` each for three different resolutions `small`, `medium`, and `large`. The models can be configured via
```jldoctest; output = false
using PortHamiltonianBenchmarkSystems
config = HeatModelConfigFactory(shape="disc", size="small")

# output
HeatModelConfig(21701)
```
and the data is loaded using `construct_system(config)`.

```@docs
HeatModelConfigFactory
HeatModelConfig
construct_system(::HeatModelConfig)
```

## References
```@bibliography
Pages = ["HeatModel.md"]
Canonical = false
```
