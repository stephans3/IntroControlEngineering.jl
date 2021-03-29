### A Pluto.jl notebook ###
# v0.12.12

using Markdown
using InteractiveUtils

# ╔═╡ b6cfe490-88a4-11eb-29d9-9f8f16950a7c
using ControlSystems

# ╔═╡ de1d432e-88c4-11eb-068b-895eb2ed291a
using LinearAlgebra

# ╔═╡ 10d0546c-88c5-11eb-2c81-896255df6f0f
using OrdinaryDiffEq

# ╔═╡ 704e295a-88c5-11eb-358e-1f2bb8c7abeb
using Plots

# ╔═╡ 193b4592-88a4-11eb-26a6-03e17972a112
md"# Stability

The general transfer function is given as

$G(s) = \frac{Y(s)}{U(s)} = \frac{b_{0} + b_{1} ~ s + b_{2} s^2 + \cdots + b_{n} ~ s^{n} }{a_{0} + a_{1} ~ s + a_{2} s^2 + \cdots + a_{m} ~ s^{m}}$

with coefficients $a_{i}, b_{j} \in \mathbb{R}$ and $i \in \left\{0,1,\cdots, m\right\}$, $j \in \left\{0,1,\cdots, n\right\}$, and the Laplace transformed system output $Y(s)$ and input $U(s)$. The system state $X(s)$ is used often instead of the system output $Y(s)$. 

The transfer function is called 

- Proper if $ n \leq m $,
- Strictly proper if $ n < m $ and
- Not proper if $ n > m $.

Here, only (strictly) proper systems are discussed because they are relevant for technical processes. 

The roots of the numerator are called **zeros**

$b_{0} + b_{1} ~ s + b_{2} s^2 + \cdots + b_{n} ~ s^{n} = 0$

and the roots of the denominator are called **poles**

$a_{0} + a_{1} ~ s + a_{2} s^2 + \cdots + a_{m} ~ s^{m} = 0 \text{.}$

The zeros have impact on the property of the [minimum phase](https://en.wikipedia.org/wiki/Minimum_phase) and the zero dynamics of a system. The zero dynamics is the internal dynamics of a system while the output is measured as zero. However, the position of the zeros are not very important for the basics of Control Engineering. 

The poles are crucial for the stability of a dynamical system. 
- If the real part of **all poles** are **smaller** than zero, then the system is called bounded-input-bounded-output (BIBO) **stable**. 
- If the real part of **any pole** is **greater** than zero, then the system is called **unstable**. 
- If the real part of **any pole** is **equal** to zero, then the system is called **marginally stable**. 

The last case is common for open loop systems with an integrator. 

The dynamics and especially the output $y(t)$ of a **stable** system will not tend to infinity for any bounded (or limited) input - or mathematical written

$\lvert y(t) \rvert < \infty ~ \text{for all} ~ t \geq 0 \text{.}$

"

# ╔═╡ 7677debe-88a4-11eb-1746-93e1ba079cea
md"### Example: Unstable Systems

The transfer functions 

$G_{1}(s) = \frac{s + 1}{s - 8} \quad \text{and} \quad
G_{2}(s) = \frac{1}{s^2 + a~s + b}$

are assumed. 

The first system $G_{1}$ has one zero root at $s = -1$ and one pole at $s = 8$. This pole is greater than zero and thus system $G_{1}$ is unstable.

The second system $G_{2}$ has no zeros and two poles 

$s_{1,2} = \frac{-a}{2} \pm \frac{1}{2} \sqrt{a^2 - 4~b} \text{.}$

If $a < 0$ or $b < 0$ then this system has at least one pole greater than zero and thus system $G_{2}$ is **unstable**. For example, the values are chosen as $a = 2$ and $b = -3$ and therefore one yields $s_{1} = -3$ and $s = 1$.

System $G_{2}(s)$ is simulated below with the final simulation time $T = 5$ seconds."

# ╔═╡ ba9c361c-88a4-11eb-1731-01b7a1b0fd23
begin
	a = 2; # Parameters
	b = -3;
end

# ╔═╡ ba6c1b62-88a4-11eb-1f54-fdef814c6ebc
G₂ = tf(1, [1, a, b]) # Transfer function of second order system

# ╔═╡ dcd8f9c2-88a4-11eb-168f-91b910b294a3
Tf = 5; # Final simulation time

# ╔═╡ e03c60ea-88a4-11eb-2fa6-59ed22a74b95
stepplot(G₂, Tf, label="y(t)")

# ╔═╡ 14ef7f16-88a5-11eb-0b79-bd5db7d56d75
md"## Marginally stable Systems

Integrator blocks have the transfer function

$G_{I}(s) = \frac{1}{s}$ 

and are often used to minimize the (steady-state) error $e(t) = r(t) - y(t)$ between reference $r(t)$ and system output $y(t)$ - see also [PID controller](https://en.wikipedia.org/wiki/PID_controller). The series connection of an integrator and a classical n-th order system will lead to an open loop system with one pole at the origin.

### Example: Integrator

The stable system 

$G_{st}(s) = \frac{1}{s^2 + 2s + 10}$

has conjugated-complex poles at $s_{1,2} = -1.0 \pm 3.0j$."

# ╔═╡ a0dfeaf0-88a6-11eb-3118-9d31dac8e68d
Gst = tf(1, [1, 2, 10]) # Transfer function of second order system

# ╔═╡ b608180a-88a6-11eb-27b2-b71e882b8e57
stepplot(Gst, Tf, label="Gst: y(t)")

# ╔═╡ a6fe594e-88a6-11eb-2bc9-e91e069ac619
md"The series connection of $G_{st}(s)$ and $G_{I}(s)$ leads to 

$G_{m}(s) = G_{I}(s) ~ G_{st}(s) =  \frac{1}{s(s^2 + 2s + 10)} \text{.}$

The system $G_{m}(s)$ has three poles at $s_{1} = 0$ and $s_{2,3} = -1.0 \pm 3.0j$ and consequently it is *only* marginally stable. It will be shown in a subsequent section that an integrator is useful to  design feedback loops."

# ╔═╡ 69309316-88a6-11eb-323e-09b6297de049
Gm = tf(1, [1, 2, 10, 0]) # Transfer function of second order system

# ╔═╡ bcd0b44c-88a6-11eb-16a3-0325ca32125e
stepplot(Gm, Tf, label="Gm: y(t)")

# ╔═╡ f2fec428-88a6-11eb-38f7-23afad0059b4
md"**Remark:** The integrator block sums up the step response of system $G_{st}(s)$!"

# ╔═╡ 68d34cba-88a6-11eb-286e-43480beb0962
md"### Example: Conjugated-complex Pole

The second order system

$G_{cc}(s) = \frac{1}{s^2 + a~s + b}$

is marginally stable, if $a = 0$ and $b > 0$. For example, if $b = 4$ is chosen, then the pure conjugated-complex poles are $s_{1,2} = \pm 2.0j$. 

System $G_{cc}(s)$ is simulated below."

# ╔═╡ 4947ff62-88a6-11eb-173f-f759f1715008
b2 = 4; # Parameter

# ╔═╡ 5284e9aa-88a6-11eb-3f68-63985696a9ba
G_cc = tf(1, [1, 0, b2]) # Transfer function of second order system

# ╔═╡ 5ae3a53c-88a6-11eb-2037-0d89e50f9e57
stepplot(G_cc, Tf, label="y(t)")

# ╔═╡ 070aba70-88a8-11eb-2968-b3f70d39a936
md"## Eigenvalues of LTI Systems

The poles of the transfer function and the Eigenvalues of dynamical systems in state-space representation are equivalent in case of linear-time invariant (LTI) systems.

All LTI systems can be noted either as transfer functions or in matrix-vector notation (state-space representation). Consider again the oscillator

$G_{2}(s) = \frac{1}{s^2 + a~s + b}$

which is noted in [state-space representation](https://en.wikipedia.org/wiki/State-space_representation#Canonical_realizations) as

$\begin{pmatrix}
\dot{x}_{1}(t) \\
\dot{x}_{2}(t)
\end{pmatrix}
=
\begin{pmatrix}
0 & 1 \\
-b & -a
\end{pmatrix} ~
\begin{pmatrix}
{x}_{1}(t) \\
{x}_{2}(t)
\end{pmatrix}
+
\begin{pmatrix}
0 \\
1
\end{pmatrix}
~u(t)$

and

$y(t) = (1~,~ 0) ~ x(t) = x_{1}(t).$

The Eigenvalues $\lambda$ of the free system (without input) $\dot{x}(t) = A ~ x(t)$ are found with

$\det( \lambda I - a ) = 0$

and thus here one yields

$\det \left(
\begin{pmatrix}
\lambda & 0 \\
0 & \lambda
\end{pmatrix}
\right)
-\begin{pmatrix}
0 & 1 \\
-b & -a
\end{pmatrix}
=
\det
\begin{pmatrix}
\lambda & -1 \\
b & \lambda+a
\end{pmatrix}
= \lambda ~ (\lambda + a) + b = 0.$

One notes that the characteristic polynomial of the transfer function 

$s^2 + a~s + b = 0$

corresponds with the characteristic polynomial of the Eigenvalues

$\lambda^2 + a~\lambda + b = 0$

"

# ╔═╡ c779701e-88c4-11eb-21ef-cbcbd9a6da91
A = [0 1; -b -a]

# ╔═╡ 9e1b60be-88c5-11eb-3356-71ee9940e5e3
md"**Eigenvalues**"

# ╔═╡ dafb4678-88c4-11eb-1b41-57345953c612
ev = eigvals(A)

# ╔═╡ fd51ce98-88c4-11eb-246e-0d16243d0db3
dyn_sys(x,p,t) = A*x

# ╔═╡ 33182644-88c5-11eb-1c63-e1cb0b3724c0
x₀ =  [1.0; 0.0]

# ╔═╡ 403d1e38-88c5-11eb-3375-9585715d2bb5
tspan = (0.0, Tf)

# ╔═╡ 17c25de2-88c5-11eb-368a-e563ff478ae4
prob = ODEProblem(dyn_sys, x₀, tspan)

# ╔═╡ 4d0ecaee-88c5-11eb-0f5c-edbd5fdf1ef1
sol = solve(prob, Tsit5())

# ╔═╡ 6b4b8d8a-88c5-11eb-15b3-8b11743c3bd0
plot(sol, label=["x1" "x2"])

# ╔═╡ Cell order:
# ╟─193b4592-88a4-11eb-26a6-03e17972a112
# ╟─7677debe-88a4-11eb-1746-93e1ba079cea
# ╠═b6cfe490-88a4-11eb-29d9-9f8f16950a7c
# ╠═ba9c361c-88a4-11eb-1731-01b7a1b0fd23
# ╠═ba6c1b62-88a4-11eb-1f54-fdef814c6ebc
# ╠═dcd8f9c2-88a4-11eb-168f-91b910b294a3
# ╠═e03c60ea-88a4-11eb-2fa6-59ed22a74b95
# ╟─14ef7f16-88a5-11eb-0b79-bd5db7d56d75
# ╠═a0dfeaf0-88a6-11eb-3118-9d31dac8e68d
# ╠═b608180a-88a6-11eb-27b2-b71e882b8e57
# ╟─a6fe594e-88a6-11eb-2bc9-e91e069ac619
# ╠═69309316-88a6-11eb-323e-09b6297de049
# ╠═bcd0b44c-88a6-11eb-16a3-0325ca32125e
# ╟─f2fec428-88a6-11eb-38f7-23afad0059b4
# ╟─68d34cba-88a6-11eb-286e-43480beb0962
# ╠═4947ff62-88a6-11eb-173f-f759f1715008
# ╠═5284e9aa-88a6-11eb-3f68-63985696a9ba
# ╠═5ae3a53c-88a6-11eb-2037-0d89e50f9e57
# ╟─070aba70-88a8-11eb-2968-b3f70d39a936
# ╠═c779701e-88c4-11eb-21ef-cbcbd9a6da91
# ╠═de1d432e-88c4-11eb-068b-895eb2ed291a
# ╟─9e1b60be-88c5-11eb-3356-71ee9940e5e3
# ╠═dafb4678-88c4-11eb-1b41-57345953c612
# ╠═fd51ce98-88c4-11eb-246e-0d16243d0db3
# ╠═33182644-88c5-11eb-1c63-e1cb0b3724c0
# ╠═403d1e38-88c5-11eb-3375-9585715d2bb5
# ╠═10d0546c-88c5-11eb-2c81-896255df6f0f
# ╠═17c25de2-88c5-11eb-368a-e563ff478ae4
# ╠═4d0ecaee-88c5-11eb-0f5c-edbd5fdf1ef1
# ╠═704e295a-88c5-11eb-358e-1f2bb8c7abeb
# ╠═6b4b8d8a-88c5-11eb-15b3-8b11743c3bd0
