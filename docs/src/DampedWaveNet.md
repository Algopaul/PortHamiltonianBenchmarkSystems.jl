# Damped Wave Net

## Description
This benchmark is a model for wave propagation in gas pipeline networks, as presented in [EKL+18](@cite). 
The network is modeled as directed, connected graph ``\mathcal{G}(\mathcal{V},\mathcal{E})``, with vertices ``v\in\mathcal{V}``, edges ``e\in\mathcal{E}`` and at least one boundary vertex ``v\in\mathcal{V}_b\subseteq\mathcal{V}``, connected to a single edge.

On each edge, the following 1D damped wave equation holds, with presssure ``p_e\ (\text{kg}\text{m}^{-1}\text{s}^{-2})``, mass flow ``m_e\ (\text{kg}\text{s}^{-1})`` and pipe constants ``a_e\ (\text{s}^2)``, ``b_e\ (\text{m}^{-2})``, ``d_e\ (\text{m}^{-2}\text{s}^{-1})``.
```math
\begin{align*}
	a_e\partial_tp_e &= -\partial_xm_e &&\forall e\in\mathcal{E},\\
	b_e\partial_tm_e &= -\partial_xp_e-d_em_e &&\forall e\in\mathcal{E},
\end{align*}
```
At each inner vertex ``v\in\mathcal{V}_i\equiv\mathcal{V}\setminus\mathcal{V}_b``, the following pressure continuity and mass conservation conditions hold, where ``p_i|_v`` is the pressure at ``v`` and ``n_e|_v`` is the direction of edge ``e`` at ``v`` (+1: incoming, -1: outgoing)
```math
\begin{align*}
	&p_e|_v \equiv p_{i}|_v &&\forall e\in\mathcal{E}(v) &&\forall v\in\mathcal{V}_i,\\
	&\sum_{e\in\mathcal{E}(v)} n_{e} m_e|_v = 0 &&\forall v\in\mathcal{V}_i.
\end{align*}
```
At each boundary vertex, either a pressure or mass flow must be fixed ``(p_{u}|_v,\ m_{u}|_v)``, leaving the other quantity to be solved for ``(p_{y}|_v,\ m_{y}|_v)``:
```math
\begin{align*}
	&p_{\mathcal{E}(v)}|_v = p_{u}|_v &n_{\mathcal{E}(v)}m_{\mathcal{E}(v)}|_v = m_{y}|_v &&\forall v \in\mathcal{V}_{b,p},\\
	&p_{\mathcal{E}(v)}|_v = p_{y}|_v &-n_{\mathcal{E}(v)}m_{\mathcal{E}(v)}|_v = m_{u}|_v &&\forall v \in\mathcal{V}_{b,m}.
\end{align*}
```
## Discretization
The Galerkin variational form of the damped wave equations can be formulated as shown below, where ``{p_e\in\text{span}(\mathcal{P}_e),\ \mathcal{P}_e=\{\pi_1\dots\pi_n\}}`` and ``m_e\in\text{span}(\mathcal{M}_e),\ \mathcal{M}_e=\{\mu_1\dots\mu_n\}``. In our implementation ``\mathcal{P}_e`` and ``\mathcal{M}_e`` are respectively discontinuous element-wise constant and continuous element-wise linear function spaces, but the shown approach holds in general:
```math
\begin{align*}
	a_e\int_e\partial_tp_e\pi_e\ \mathrm{d}x &= -\int_e\partial_xm_e\pi_e\ \mathrm{d}x &&\forall \pi_e\in\mathcal{P}_e,\\
	b_e\int_e\partial_tm_e\mu_e\ \mathrm{d}x &= -[p_e\mu_e]_{\partial e} + \int_ep_e\partial_x\mu_e\ \mathrm{d}x -d_e\int_em_e\mu_e\ \mathrm{d}x &&\forall \mu_e\in\mathcal{M}_e.
\end{align*}
```
The boundary terms ``-[p_e\mu_e]_{\partial e}`` are the result of integration by parts and can be grouped by vertex as follows:
```math
\begin{align*}
	-n_{e} \mu_ep_{i}|_v && \forall \mu_e \in\mathcal{M}_e\quad \forall e\in\mathcal{E}(v)\quad \forall v\in\mathcal{V}_i,\\
	-n_{\mathcal{E}(v)} \mu_{\mathcal{E}(v)}p_{u}|_v &&\forall \mu_e \in\mathcal{M}_e\quad \forall v \in\mathcal{V}_{b,p},\\
	-n_{\mathcal{E}(v)} \mu_{\mathcal{E}(v)}p_y|_v &&\forall \mu_e \in\mathcal{M}_e\quad \forall v \in\mathcal{V}_{b,m}.
\end{align*}
```
It now becomes apparent that in matrix form, the linear operators in several pairs of terms are each other's (negative) transpose. Hence, we can write our system of equations as the following linear DAE:
```math
\begin{align*}
    \underbrace{
    \begin{bmatrix}
        A_pM_p & & & \\
        & B_mM_m & & \\
        & & 0 & \\
        & & & 0
    \end{bmatrix}}_{E}
    \begin{bmatrix}
        \dot{p}\\
        \dot{m}\\
        \dot{p}_i\\
        \dot{p}_y
    \end{bmatrix} &=
    \underbrace{
    \begin{bmatrix}
        & -G_m & & \\
        G_m^\mathsf{T}& -D_mM_m &-C_{m}^\mathsf{T} & U_{m}^\mathsf{T}\\
        & C_{m}& & \\
        & -U_{m}& &
    \end{bmatrix}}_{A}
    \begin{bmatrix}
        p\\
        m\\
        p_i\\
        p_y
    \end{bmatrix} +
    \underbrace{
    \begin{bmatrix}
        0 &\\
        Y_{m}^\mathsf{T} &\\
        & 0 \\
        & I
    \end{bmatrix}}_{B}
    \begin{bmatrix}
        p_u\\
        m_u
    \end{bmatrix},\\
    \begin{bmatrix}
        m_y\\
        p_y
    \end{bmatrix} &=
    \underbrace{
    \begin{bmatrix}
        0  &Y_{m} & &\\
        & & 0  & I
    \end{bmatrix}}_{B^\mathsf{T}} 
    \begin{bmatrix}
        p\\
        m\\
        p_i\\
        p_y
    \end{bmatrix}.
\end{align*}
```
- ``M_p,\ M_m``: mass matrices for ``p,\ m``
- ``A_p,\ B_m,\ D_m``: diagonal matrices containing the edge parameters ``a_e,\ b_e,\ d_e``
- ``G_m``: Galerkin variational operator for ``\partial_xm``
- ``C_m``: mass conservation conditions for ``m``
- ``U_m,\ Y_m``: matrices selecting ``m_u,\ m_y`` from ``m``

Since ``p`` contains all the pressure variables, ``p_i`` and ``p_y`` are redundant in the solution vector. However, they are not explicitly tied to ``p`` in the system. It can be proven that the system has a unique solution and that this constrains ``p_i`` and ``p_y`` to be equal to their counterparts in ``p``, ensuring that the original variational problem is solved (cf. [EKL+18](@cite)).

Finally, the system can be written in linear port-Hamiltonian form as follows:
```math
\begin{align*}
    E\dot{x} &= (J-R)Qx + (G-P)u,\\
    y &= (G+P)^\mathsf{T}Qx + (S+N)u,
\end{align*}
```
where ``J = \frac{1}{2}(A-A^\mathsf{T})``, ``R = -\frac{1}{2}(A+A^\mathsf{T})``, ``Q = I``, ``G = B``, ``P = 0``, ``S = N = 0``.

## Interface
```@docs
DampedWaveNetConfig
```
```@docs
DampedWaveNetConfig(id::String)
```
```@docs
construct_system(problem::DampedWaveNetConfig)
```

## References
```@bibliography
Pages = ["DampedWaveNet.md"]
Canonical = false
```
