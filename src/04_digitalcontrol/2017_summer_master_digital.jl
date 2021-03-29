### A Pluto.jl notebook ###
# v0.12.12

using Markdown
using InteractiveUtils

# ╔═╡ 98857f3e-5b0e-11eb-2db2-fb2a50161f8c
md"# Example: Digital Control in the State-Space

A system in state-space representation

$\dot{x}(t) = A x(t) + B u(t) \quad \text{with} \quad x(t_{0}) = x_{0}$

is solved via integration in time as

$x(t) = \exp \left( A(t - t_{0}) \right) x(t_{0}) + \int\limits_{t_{0}}^{t} \exp \left( A(t - \tau) \right) B ~ u(\tau) d\tau \text{.}$

This solution in time is used to derive the discrete state-space representation. If the $x(t)$ is solved iteratively only in the time interval 

$t \in \left\{ t_{0} + i~\Delta T  , t_{0} + (i+1)~\Delta T \right\}$ 

for an increasing $i \in \left\{ 0, 1, 2, \cdots, N-1 \right\}$ than the solution is found from $t_{0} $ until $ t_{final} = N ~ \Delta T$. This incremental solution is expressed by the difference equation 

$x(t_{n+1}) = \exp \left( A \Delta T \right) x(t_{n}) + \int\limits_{t_{n}}^{t_{n+1}} \exp \left( A(t_{n+1} - \tau) \right) ~ B ~ d\tau ~ u(t_{n}) 
    = A_{d} x(t_{n}) + B_{d} u(t_{n})$

with $\Delta T = t_{n+1} - t_{n}$. The linear time-invariant system 

$\begin{pmatrix}
    \dot{x}_{1}(t) \\
    \dot{x}_{2}(t) 
    \end{pmatrix}
    =
    \begin{pmatrix}
    -0.5 & 0 \\
    1 & 1 
    \end{pmatrix}
    ~
    \begin{pmatrix}
    x_{1}(t) \\
    x_{2}(t) 
    \end{pmatrix}
    +
    \begin{pmatrix}
    0.5 \\
    0.5 
    \end{pmatrix}
    ~
    u(t)$

is transformed to a [time-discrete dynamical system](https://en.wikipedia.org/wiki/Discretization) and a digital controller is designed to stabilize it.

## Stability of the continuous system

The stability is proved with position of the Eigenvalues which are calculated with

$\det \left(A - \lambda I\right) =
\det
\begin{pmatrix}
-0.5 - \lambda & 0 \\
1 & 1 - \lambda
\end{pmatrix} 
= \left(\lambda + \frac{1}{2} \right) \left(\lambda - 1 \right) = 0$

and finally one holds $\lambda_{1} = -0.5$ and $\lambda_{2} = 1.0$ . The second Eigenvalue unveils an **unstable** dynamic behaviour."

# ╔═╡ Cell order:
# ╠═98857f3e-5b0e-11eb-2db2-fb2a50161f8c
