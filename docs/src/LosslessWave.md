# The lossless wave equation with Neumann boundary control

## The model

Let us consider the vertical deflection from equilibrium $w$ of a 2D membrane $\Omega \subset \mathbb{R}^2$. Denoting $\rho$ the mass density and $T$ the Young modulus of the membrane, a positive definite tensor, leads to the following well-known *wave equation*

```math
\rho(x) \frac{\partial^2}{\partial t^2} w(t,x) - {\rm div} \left( T(x) \cdot {\rm grad} \left( w(t,x) \right) \right) = 0, \quad t \ge 0, \, x \in \Omega,
```

together with *Neumann boundary control* 

```math
\left( T(x) \cdot {\rm grad} \left( w(t,x) \right) \right) \cdot \mathbf{n} = u_\partial(t,x), \quad t \ge 0, \, x \in \partial \Omega,
```

where $\mathbf{n}$ is the outward normal to $\Omega$.

The **Hamiltonian** is the total mechanical energy, given as the sum of potential and kinetic energies

```math
\mathcal{H}(t) := \frac{1}{2} \int_{\Omega} \left( {\rm grad} \left( w(t,x) \right) \right)^\top \cdot T(x) \cdot {\rm grad} \left( w(t,x) \right) {\rm d}x + \frac{1}{2} \int_{\Omega} \rho(x) \left( \frac{\partial}{\partial t} w(t,x) \right)^2 {\rm d}x, \qquad t \ge 0.
```

Taking the *strain* and the *linear momentum*

```math
\alpha_q := {\rm grad} \left( w \right), \qquad \alpha_p := \frac{\partial}{\partial t} w,
```

as **energy variables**, the Hamiltonian rewrites

```math
\mathcal{H}(t) = \mathcal{H}(\alpha_q (t,\cdot), \alpha_p(t,\cdot)) = \frac{1}{2} \int_{\Omega} \left( \alpha_q(t,x) \right)^\top \cdot T(x) \cdot \alpha_q(t,x) {\rm d}x + \frac{1}{2} \int_{\Omega} \frac{\alpha_p^2(t,x)}{\rho(x)} {\rm d}x.
```

The **co-energy variables** are by definition the variational derivatives of the Hamiltonian

```math
e_q := \delta_{\alpha_q} \mathcal{H} = T \cdot \alpha_q, 
\qquad e_p := \delta_{\alpha_p} \mathcal{H} = \frac{\alpha_p}{\rho},
```

*i.e.* the *stress* and the *velocity* respectively. These equality are the **constitutive relation** which close the system.

Thanks to these variables, the wave equation writes as a **port-Hamiltonian system**

```math
\begin{pmatrix}
\frac{\partial}{\partial t} \alpha_q \\
\frac{\partial}{\partial t} \alpha_p
\end{pmatrix} =
\begin{bmatrix}
0 & {\rm grad} \\
{\rm div} & 0
\end{bmatrix}
\begin{pmatrix}
e_q \\
e_p
\end{pmatrix}, 
\qquad \left\lbrace
\begin{array}{rcl}
e_q &=& T \cdot \alpha_q, \\
e_p &=& \frac{\alpha_p}{\rho},
\end{array}\right.
\qquad \left\lbrace
\begin{array}{rcl}
u_\partial &=& e_q \cdot \mathbf{n}, \\
y_\partial &=& e_p\mid_{\partial \Omega},
\end{array}\right.
```

The **power balance** satisfied by the Hamiltonian is

```math
\frac{\rm d}{ {\rm d}t} \mathcal{H} = \langle u_\partial, y_\partial \rangle_{H^{-\frac12}(\partial \Omega),H^\frac12(\partial \Omega)}
```

To get rid of the algebraic constraints induced by the constitutive relations, one rewrites the port-Hamiltonian system as

```math
\begin{bmatrix}
T^{-1} & 0 \\
0 & \rho
\end{bmatrix}
\begin{pmatrix}
\frac{\partial}{\partial t} e_q \\
\frac{\partial}{\partial t} e_p
\end{pmatrix} =
\begin{bmatrix}
0 & {\rm grad} \\
{\rm div} & 0
\end{bmatrix}
\begin{pmatrix}
e_q \\
e_p
\end{pmatrix}, 
\qquad \left\lbrace
\begin{array}{rcl}
u_\partial &=& e_q \cdot \mathbf{n}, \\
y_\partial &=& e_p\mid_{\partial \Omega},
\end{array}\right.
```

known as the **co-energy formulation**. This allows to get a simple Ordinary Differential Equation at the discrete level (instead of a Differential Algebraic Equation in general).


## Structure-preserving discretization


### Weak formulation

Let $\phi_q$, $\varphi_p$ and $\psi$ be vector-valued, scalar-valued and boundary scalar-valued test functions respectively. The weak formulation reads

```math
\left\lbrace
\begin{array}{rcl}
\displaystyle \int_{\Omega} \phi_q \cdot T^{-1} \cdot \frac{\partial}{\partial t} e_q 
&=& \displaystyle \int_{\Omega} \phi_q \cdot {\rm grad} \left( e_p \right), \\
\displaystyle \int_{\Omega} \varphi_p \rho \frac{\partial}{\partial t} e_p 
&=& \displaystyle \int_{\Omega} \varphi_p {\rm div} \left( e_q \right), \\
\displaystyle \int_{\partial \Omega} \psi y_\partial &=& \displaystyle \int_{\partial \Omega} \psi e_p.
\end{array}\right.
```



### Integration by parts

The integration by parts of the second line makes $u_\partial = e_q \cdot \mathbf{n}$ appear

```math
\left\lbrace
\begin{array}{rcl}
\displaystyle \int_{\Omega} \phi_q \cdot T^{-1} \cdot \frac{\partial}{\partial t} e_q 
&=& \displaystyle \int_{\Omega} \phi_q \cdot {\rm grad} \left( e_p \right), \\
\displaystyle \int_{\Omega} \varphi_p \rho \frac{\partial}{\partial t} e_p 
&=& \displaystyle - \int_{\Omega} {\rm grad} \left( \varphi_p \right) \cdot e_q + \int_{\partial \Omega} \varphi_p u_\partial, \\
\displaystyle \int_{\partial \Omega} \psi y_\partial &=& \displaystyle \int_{\partial \Omega} \psi e_p.
\end{array}\right.
```



### Projection

Let $(\phi_q^i)\_{1 \le i \le N_q}$, $(\varphi_p^j)\_{1 \le j \le N_p}$ and $(\psi^k)\_{1 \le k \le N_\partial}$ be finite element families for $q$-type, $p$-type and boundary-type variables. Variables are approximated in their respective finite element family

```math
e_q^d(t,x) := \sum_{i=1}^{N_q} e_q^i(t) \phi_q^i(x),
\qquad e_p^d(t,x) := \sum_{j=1}^{N_p} e_p^j(t) \varphi_p^j(x),
```

```math
u_\partial^d(t,x) := \sum_{k=1}^{N_\partial} u_\partial^k(t) \psi^k(x),
\qquad y_\partial^d(t,x) := \sum_{k=1}^{N_\partial} y_\partial^k(t) \psi^k(x).
```

Denoting $\underline{\star}$ the (time-varying) vector of coordinates of the discretisation $\star^d$ of $\star$ in its respective finite element family, the discrete system reads

```math
\underset{M}{\underbrace{\begin{bmatrix}
M_q & 0 & 0 \\
0 & M_p & 0 \\
0 & 0 & M_\partial
\end{bmatrix} } }
\begin{pmatrix}
\frac{\rm d}{ {\rm d}t} \underline{e_q}(t) \\
\frac{\rm d}{ {\rm d}t} \underline{e_p}(t) \\
\ - \underline{y_\partial}(t)
\end{pmatrix} =
\underset{J}{\underbrace{\begin{bmatrix}
0 & D & 0 \\
\ -D^\top & 0 & B \\
0 & -B^\top & 0
\end{bmatrix} } }
\begin{pmatrix}
\underline{e_q}(t) \\
\underline{e_p}(t) \\
\underline{u_\partial}(t)
\end{pmatrix}
```

where

```math
(M_q)\_{ij} := \int_{\Omega} \phi_q^i \cdot T^{-1} \cdot \phi_q^j,
\qquad 
(M_p)\_{ij} := \int_{\Omega} \varphi_p^i \rho \varphi_p^j,
\qquad 
(M_\partial)\_{ij} := \int_{\partial \Omega} \psi^i \psi^j,
```

and

```math
(D)\_{ij} := \int_{\Omega} \phi_q^i \cdot {\rm grad} \left( \varphi_p^j \right),
\qquad
(B)\_{jk} := \int_{\partial \Omega} \varphi_p^j \psi^k.
```



### Discrete Hamiltonian

By definition, the discrete Hamiltonian is equal to the continuous Hamiltonian evaluated in the approximated variables. As we are working with the co-energy formulation, a first step is to restate the Hamiltonian in terms of co-energy variables

```math
\mathcal{H} = \frac{1}{2} \int_{\Omega} e_q \cdot T^{-1} \cdot e_q + \frac{1}{2} \int_\Omega \rho (e_p)^2.
```

Then, the discrete Hamiltonian is defined as

```math
\mathcal{H}^d := \frac{1}{2} \int_{\Omega} e_q^d \cdot T^{-1} \cdot e_q^d + \frac{1}{2} \int_{\Omega} \rho (e_p^d)^2.
```

After straightforward computations, it comes

```math
\mathcal{H}^d(t) = \frac{1}{2} \underline{e_q}(t)^\top M_q \underline{e_q}(t) + \frac{1}{2} \underline{e_p}(t)^\top M_p \underline{e_p}(t),
```

and the **discrete power balance** follows

```math
\frac{\rm d}{ {\rm d}t} \mathcal{H}^d(t) = \underline{u_\partial}(t)^\top M_\partial \underline{y_\partial}(t).
```

## Interface

We have generated models in domains of three different shapes `disc`, `rectangle`, and `L` each for three different resolutions `small`, `medium`, and `large`. The models can be configured via
```jldoctest; output = false
using PortHamiltonianBenchmarkSystems
config = LosslessWaveModelConfigFactory(shape="disc", size="small")

# output
LosslessWaveModelConfig(21701)
```
and the data is loaded using `construct_system(config)`.

```@docs
LosslessWaveModelConfigFactory
LosslessWaveModelConfig
construct_system(::LosslessWaveModelConfig)
```

## References


```
@incollection{Serhani2019a,
author = {Serhani, Anass and Matignon, Denis and Haine, Ghislain},
booktitle = {Geometric Science of Information},
editor = { {Nielsen, Frank} and {Barbaresco, Fr{\'e}d{\'e}ric} },
pages = {549--558},
publisher = {Springer},
address ={Cham},
series = {Lecture Notes in Computer Science},
title = { {A Partitioned Finite Element Method for the Structure-Preserving Discretization of Damped Infinite-Dimensional Port-Hamiltonian Systems with Boundary Control} },
volume = {11712},
year = {2019}
}
```

```
@article{Serhani2019b,
author={Serhani, Anass and Matignon, Denis and Haine, Ghislain},
journal = {IFAC-PapersOnLine},
title={ {Partitioned Finite Element Method for port-Hamiltonian systems with Boundary Damping: Anisotropic Heterogeneous 2-D wave equations} },
volume = {52},
number = {2},
pages = {96--101},
year = {2019},
note = {3rd IFAC Workshop on Control of Systems Governed by Partial Differential Equations (CPDE)}
}
```
