# Damped Wave Net

## Description

```math
\gdef\kett#1{\mathnormal{#1}}
```

This benchmark is a model for wave propagation in gas pipeline networks, as presented in ([EKLSMM2018](#References)). The network is modeled as directed, connected graph ``\mathcal{G}(\mathcal{V},\mathcal{E})``, with vertices ``v\in\mathcal{V}``, edges ``e\in\mathcal{E}`` and at least one boundary vertex ``v\in\mathcal{V}_b\subseteq\mathcal{V}``, connected to a single edge.

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
	&p_{\mathcal{E}(v)}|_v = p_{u}|_v &n_{\mathcal{E}(v)}m_{\mathcal{E}(v)}|_v = m_{y}|_v &&\forall v \in\mathcal{V}_{b,p}\\
	&p_{\mathcal{E}(v)}|_v = p_{y}|_v &-n_{\mathcal{E}(v)}m_{\mathcal{E}(v)}|_v = m_{u}|_v &&\forall v \in\mathcal{V}_{b,m}
\end{align*}
```
## Discretization

The Galerkin variational form of the damped wave equations can be formulated as shown below, where ``{p_e\in\text{span}(\mathcal{P}_e),\ \mathcal{P}_e=\{\pi_1\dots\pi_n\}}`` and ``m_e\in\text{span}(\mathcal{M}_e),\ \mathcal{M}_e=\{\mu_1\dots\mu_n\}``. In our implementation ``\mathcal{P}_e`` and ``\mathcal{M}_e`` are respectively discontinuous element-wise constant and continuous element-wise linear function spaces, but the shown approach holds in general:
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
        \kett{A}_p\kett{M}_p & & & \\
        & \kett{B}_m\kett{M}_m & & \\
        & & 0 & \\
        & & & 0
    \end{bmatrix}}_{\kett{E}}
    \begin{bmatrix}
        \kett{\dot{p}}\\
        \kett{\dot{m}}\\
        \kett{\dot{p}}_i\\
        \kett{\dot{p}}_y
    \end{bmatrix} &=
    \underbrace{
    \begin{bmatrix}
        & -\kett{G}_m & & \\
        \kett{G}_m^T& -\kett{D}_m\kett{M}_m &-\kett{C}_{m}^T & \kett{U}_{m}^T\\
        & \kett{C}_{m}& & \\
        & -\kett{U}_{m}& &
    \end{bmatrix}}_{\kett{A}}
    \begin{bmatrix}
        \kett{p}\\
        \kett{m}\\
        \kett{p}_i\\
        \kett{p}_y
    \end{bmatrix} +
    \underbrace{
    \begin{bmatrix}
         &\\
        \kett{Y}_{m}^T &\\
        &  \\
        & \kett{I}
    \end{bmatrix}}_{\kett{B}}
    \begin{bmatrix}
        \kett{p}_u\\
        \kett{m}_u
    \end{bmatrix}\\
    \begin{bmatrix}
        \kett{m}_y\\
        \kett{p}_y
    \end{bmatrix} &=
    \underbrace{
    \begin{bmatrix}
        0  &\kett{Y}_{m} & &\\
        & & 0  & \kett{I}
    \end{bmatrix}}_{\kett{B}^T} 
    \begin{bmatrix}
        \kett{p}\\
        \kett{m}\\
        \kett{p}_i\\
        \kett{p}_y
    \end{bmatrix}
\end{align*}
```
- ``\kett{M}_p,\ \kett{M}_m``: mass matrices for ``p,\ m``
- ``\kett{A}_p,\ \kett{B}_m,\ \kett{D}_m``: diagonal matrices containing the edge parameters ``a_e,\ b_e,\ d_e``
- ``\kett{G}_m``: Galerkin variational operator for ``\partial_xm``
- ``\kett{C}_m``: mass conservation conditions for ``m``
- ``\kett{U}_m,\ \kett{Y}_m``: matrices selecting ``\kett{m}_u,\ \kett{m}_y`` from ``\kett{m}``

Since ``\kett{p}`` contains all the pressure variables, ``\kett{p}_i`` and ``\kett{p}_y`` are redundant in the solution vector. However, they are not explicitly tied to ``\kett{p}`` in the system. It can be proven that the system has a unique solution and that this constrains ``\kett{p}_i`` and ``\kett{p}_y`` to be equal to their counterparts in ``\kett{p}``, ensuring that the original variational problem is solved ([Egger, 2018](#References)).

Finally, the system can be written in linear port-Hamiltonian form as follows:
```math
\begin{align*}
    \begin{matrix*}[l]
        \kett{E}\kett{\dot{x}} = (\kett{J}-\kett{R})\kett{Q}\kett{x} + (\kett{G}-\kett{P})\kett{u}\\[0.33em]
        \ \ \ \ \kett{y} = (\kett{G}+\kett{P})^H\kett{Q}\kett{x} + (\kett{S}+\kett{N})\kett{u}
    \end{matrix*} \quad\quad\quad\quad
    \begin{cases}
        \kett{J} = \frac{1}{2}(\kett{A}-\kett{A}^\mathsf{T})\\
        \kett{R} = -\frac{1}{2}(\kett{A}+\kett{A}^\mathsf{T})\\
        \kett{Q} = \kett{I}\\
        \kett{G} = \kett{B}\\
        \kett{P} = 0\\
        \kett{S} = 0\\
        \kett{N} = 0
    \end{cases}
\end{align*}
```

## Interface
```@docs
DampedWaveNet
```
```@docs
DampedWaveNet(id::String)
```
```@docs
construct_system(problem::DampedWaveNet)
```

## References
```LaTeX
@article{EKLSMM2018,
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
