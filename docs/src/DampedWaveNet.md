# Damped Wave Net

## Description

This benchmark is a model for wave propagation in gas pipeline networks, as presented in ([Egger, 2018](#References)). The network is modeled as directed, connected graph ``\mathcal{G}(\mathcal{V},\mathcal{E})``, with vertices ``v\in\mathcal{V}``, edges ``e\in\mathcal{E}`` and at least one boundary vertex ``v\in\mathcal{V}_b\subseteq\mathcal{V}``, connected to a single edge.

Graph Figure

On each edge, the following 1D damped wave equation holds, with presssure ``p_e\ (\text{kg}\text{m}^{-1}\text{s}^{-2})``, mass flow ``m_e\ (\text{kg}\text{s}^{-1})`` and pipe constants ``a_e\ (\text{s}^2)``, ``b_e\ (\text{m}^{-2})``, ``d_e\ (\text{m}^{-2}\text{s}^{-1})``:
```math
\begin{align*}
	a_e\partial_tp_e &= -\partial_xm_e &&\forall e\in\mathcal{E}\\
	b_e\partial_tm_e &= -\partial_xp_e-d_em_e &&\forall e\in\mathcal{E}
\end{align*}
```
At each inner vertex ``v\in\mathcal{V}_i\equiv\mathcal{V}\setminus\mathcal{V}_b``, the following pressure continuity and mass conservation conditions hold, where ``p_i|_v`` is the pressure at ``v`` and ``n_e|_v`` is the direction of edge ``e`` at ``v`` (+1: incoming, -1: outgoing):
```math
\begin{align*}
	&p_e|_v \equiv p_{i}|_v &&\forall e\in\mathcal{E}(v) &&\forall v\in\mathcal{V}_i\\
	&\sum_{e\in\mathcal{E}(v)} n_{e} m_e|_v = 0 &&\forall v\in\mathcal{V}_i
\end{align*}
```
At each boundary vertex, either a pressure or mass flow must be fixed ``(p_{u}|_v,\ m_{u}|_v)``, leaving the other quantity to be solved for ``(p_{y}|_v,\ m_{y}|_v)``:
```math
\begin{align*}
	p_{\mathcal{E}(v)}|_v = p_{u}|_v\quad n_{\mathcal{E}(v)}m_{\mathcal{E}(v)}|_v = m_{y}|_v &&\forall v \in\mathcal{V}_{b,p}\\
	p_{\mathcal{E}(v)}|_v = p_{y}|_v\quad n_{\mathcal{E}(v)}m_{\mathcal{E}(v)}|_v = m_{u}|_v &&\forall v \in\mathcal{V}_{b,m}
\end{align*}
```
## Discretization

The Galerkin variational form of the damped wave equations can be formulated as follows, where ``p_e\in\text{span}\ \mathcal{P}_e,\ \mathcal{P}_e=\{\pi_1\dots\pi_n\}``, ``m_e\in\text{span}\ \mathcal{M}_e,\ \mathcal{M}_e=\{\mu_1\dots\mu_n\}``:
```math
\begin{align*}
	a_e\int_e\partial_tp_e\pi_e\ dx &= -\int_e\partial_xm_e\pi_e\ dx &&\forall \pi_e\in\mathcal{P}_e\\
	b_e\int_e\partial_tm_e\mu_e\ dx &= -[p_e\mu_e]_{\partial e}+\int_ep_e\partial_x\mu_e\ dx -d_e\int_em_e\mu_e\ dx &&\forall \mu_e\in\mathcal{M}_e
\end{align*}
```
The boundary terms ``-[p_e\mu_e]_{\partial e}`` are the result of integration by parts and can be grouped by vertex as follows:
```math
\begin{align*}
	-n_{e} \mu_ep_{i}|_v && \forall \mu_e \in\mathcal{M}_e\quad \forall e\in\mathcal{E}(v)\quad \forall v\in\mathcal{V}_i\\
	-n_{\mathcal{E}(v)} \mu_{\mathcal{E}(v)}p_{u}|_v &&\forall \mu_e \in\mathcal{M}_e\quad \forall v \in\mathcal{V}_{b,p}\\
	-n_{\mathcal{E}(v)} \mu_{\mathcal{E}(v)}p_y|_v &&\forall \mu_e \in\mathcal{M}_e\quad \forall v \in\mathcal{V}_{b,m}
\end{align*}
```
It now becomes apparent that in matrix form, the linear operators in several pairs of terms are each other's (negative) transpose. Hence, we can write our system of equations as the following linear DAE:
```math
\begin{align*}
    \underbrace{
    \begin{bmatrix}
        \mathbf{A}_p\mathbf{M}_p & & & \\
        & \mathbf{B}_m\mathbf{M}_m & & \\
        & & \mathbf{0}& \\
        & & & \mathbf{0}
    \end{bmatrix}}_{\mathbf{E}}
    \begin{bmatrix}
        \mathbf{\dot{p}}\\
        \mathbf{\dot{m}}\\
        \mathbf{\dot{p}}_i\\
        \mathbf{\dot{p}}_y
    \end{bmatrix} &=
    \underbrace{
    \begin{bmatrix}
        & -\mathbf{G}_m & & \\
        \mathbf{G}_m^T& -\mathbf{D}_m\mathbf{M}_m &-\mathbf{C}_{m}^T & \mathbf{U}_{m}^T\\
        & \mathbf{C}_{m}& & \\
        & -\mathbf{U}_{m}& &
    \end{bmatrix}}_{\mathbf{A}}
    \begin{bmatrix}
        \mathbf{p}\\
        \mathbf{m}\\
        \mathbf{p}_i\\
        \mathbf{p}_y
    \end{bmatrix} +
    \underbrace{
    \begin{bmatrix}
        \mathbf{0}&\\
        \mathbf{Y}_{m}^T &\\
        & \mathbf{0}\\
        & \mathbf{I}
    \end{bmatrix}}_{\mathbf{B}}
    \begin{bmatrix}
        \mathbf{p}_u\\
        \mathbf{m}_u
    \end{bmatrix}\\
    \begin{bmatrix}
        \mathbf{m}_y\\
        \mathbf{p}_y
    \end{bmatrix} &=
    \underbrace{
    \begin{bmatrix}
        \mathbf{0} &\mathbf{Y}_{m} & &\\
        & & \mathbf{0} & \mathbf{I}
    \end{bmatrix}}_{\mathbf{B}^T} 
    \begin{bmatrix}
        \mathbf{p}\\
        \mathbf{m}\\
        \mathbf{p}_i\\
        \mathbf{p}_y
    \end{bmatrix}
\end{align*}
```
- ``\mathbf{M}_p,\ \mathbf{M}_m``: mass matrices for ``p,\ m``
- ``\mathbf{A}_p,\ \mathbf{B}_m,\ \mathbf{D}_m``: diagonal matrices containing the edge parameters ``a_e,\ b_e,\ d_e``
- ``\mathbf{G}_m``: matrix coming from the gradient of ``m``
- ``\mathbf{C}_m``: matrix coming from the mass conservation conditions on ``m``
- ``\mathbf{U}_m,\ \mathbf{Y}_m``: matrices  selecting ``\mathbf{m}_u,\ \mathbf{m}_y`` from ``\mathbf{m}``

Since ``\mathbf{p}`` contains all the pressure variables, ``\mathbf{p}_i`` and ``\mathbf{p}_y`` are redundant in the solution vector. However, they are not explicitly tied to ``\mathbf{p}`` in the system. It can be proven that the system has a unique solution and that this constrains ``\mathbf{p}_i`` and ``\mathbf{p}_y`` to be equal to their counterparts in ``\mathbf{p}``, ensuring that the original variational problem is solved ([Egger, 2018](#References)).

The system can be written in standard linear port-Hamiltonian form as follows:
```math
\begin{align*}
    \mathbf{E}\mathbf{\dot{x}} &= (\mathbf{J}-\mathbf{R})\mathbf{Q}\mathbf{x} + (\mathbf{G}-\mathbf{P})\mathbf{u}\\
    \mathbf{y} &= (\mathbf{G}+\mathbf{P})^H\mathbf{Q}\mathbf{x} + (\mathbf{S}+\mathbf{N})\mathbf{u}
\end{align*}
```
- ``\mathbf{J} = \frac{1}{2}(\mathbf{A}-\mathbf{A}^H)``
- ``\mathbf{R} = -\frac{1}{2}(\mathbf{A}+\mathbf{A}^H)``
- ``\mathbf{Q} = \mathbf{I}``
- ``\mathbf{G} = \mathbf{B}``
- ``\mathbf{P} = \mathbf{S} = \mathbf{N} = \mathbf{0}``

## Interface
```@docs
DampedWaveNet
```

## References
```LaTeX
@article{
	Egger2018,
	author = {Egger, H. and Kugler, T. and Liljegren-Sailer, B. and Marheineke, N. and Mehrmann, V.},
	title = {On Structure-Preserving Model Reduction for Damped Wave Propagation in Transport Networks},
	journal = {SIAM Journal on Scientific Computing},
	volume = {40},
	number = {1},
	pages = {A331-A365},
	year = {2018},
	doi = {10.1137/17M1125303},
	URL = {https://doi.org/10.1137/17M1125303},
	eprint = {https://doi.org/10.1137/17M1125303}
}
```
