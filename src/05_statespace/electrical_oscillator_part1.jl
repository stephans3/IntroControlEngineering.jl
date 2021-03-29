### A Pluto.jl notebook ###
# v0.12.12

using Markdown
using InteractiveUtils

# ╔═╡ 93468f12-7c0c-11eb-08f0-b51375661637
using LinearAlgebra

# ╔═╡ 012acde6-7c1b-11eb-3cb8-4b5dd45c4194
using DifferentialEquations

# ╔═╡ 8e2e24ea-7c1d-11eb-2a8e-1592f75ae0b3
using Plots

# ╔═╡ 7f35a26e-7c08-11eb-1c76-1d1d7ef04e73
md"# Electrical Oscillator I: State Controller

The transfer function of a simple [RC circuit](https://en.wikipedia.org/wiki/RC_circuit) is derived as 

$G(s) = \frac{1}{T~s +1}$

and the series connection of three RC circuits is found by multiplication as

$G(s) = G_{1}(s) ~ G_{2}(s) ~ G_{3}(s) = \frac{1}{(T_{1}~s +1) (T_{2}~s +1) (T_{3}~s +1)}. \tag{1}$

Transfer function $(1)$ is noted in polynomial form as

$G(s) = \frac{1}{T_{1} T_{2} T_{3} ~ s^{3} + (T_{1} T_{2} + T_{1} T_{3} + T_{2} T_{3}) ~s^2 + (T_{1} + T_{2} + T_{3})~s +1}$

and further reformulated as 

$G(s) = \frac{n_{0}}{s^3 + d_{2} s^2 + d_{1} s + d_{0}} \tag{2}$

with 

$n_{0} = d_{0} = (T_{1} ~ T_{2} ~ T_{3})^{-1} \text{,}$

$d_{1} = \frac{T_{1} + T_{2} + T_{3}}{T_{1} ~ T_{2} ~ T_{3}}$

and

$d_{2} = \frac{T_{1} T_{2} + T_{1} T_{3} + T_{2} T_{3}}{T_{1} ~ T_{2} ~ T_{3}} \text{.}$
"

# ╔═╡ 7c9fa30a-7c18-11eb-051d-3db5c1bcf4f8
R = 10^4 # Resistivity: 10kΩ 

# ╔═╡ 7c835e3e-7c18-11eb-31ba-2d5768996c64
C = 100 * 10^(-6) # Capacity 100μF

# ╔═╡ f94fa7f8-7c18-11eb-1621-ad32e0aa3311
T₁ = T₂ = T₃ = R*C # Time constants

# ╔═╡ 3f48146e-7c19-11eb-3450-f7d2340fc4cc
d₀ = 1/(T₁ * T₂ * T₃)

# ╔═╡ 5512f8fe-7c19-11eb-3c00-6173a632e2fd
d₁ = (T₁ + T₂ + T₃)/(T₁ * T₂ * T₃)

# ╔═╡ 661d921c-7c19-11eb-2b55-eb1ad3cd9aa0
d₂ = (T₁*T₂ + T₁*T₃ + T₂*T₃)/(T₁ * T₂ * T₃)

# ╔═╡ 9382e638-7c0c-11eb-348e-e3bf484b7045
md"# State-space representation

The state-space model

$\dot{x}(t) = A x(t) + B u(t)$

$y(t) = C x(t)$

is derived from transfer function $(2)$ as

$\begin{pmatrix}
\dot{x}_{1}(t) \\
\dot{x}_{2}(t) \\
\dot{x}_{3}(t)
\end{pmatrix}
=
\begin{pmatrix}
0 & 1 & 0 \\
0 & 0 & 1 \\
-d_{0} & -d_{1} & -d_{2}
\end{pmatrix}
~
\begin{pmatrix}
{x}_{1}(t) \\
{x}_{2}(t) \\
{x}_{3}(t)
\end{pmatrix}
+
\begin{pmatrix}
0 \\
0 \\
1
\end{pmatrix}
~ u(t) \tag{3}$

with output 

$y(t) = 
\begin{pmatrix}
n_{0} & 0 & 0
\end{pmatrix} 
\begin{pmatrix}
{x}_{1}(t) \\
{x}_{2}(t) \\
{x}_{3}(t)
\end{pmatrix}
= d_{0} ~ x_{1}(t). \tag{4}$


The dynamical system $(3)$ is stable because all Eigenvalues $\lambda_{i}$ of system matrix 

$A = \begin{pmatrix}
0 & 1 & 0 \\
0 & 0 & 1 \\
-d_{0} & -d_{1} & -d_{2}
\end{pmatrix}$

are smaller than zero.

"

# ╔═╡ 9367d0be-7c0c-11eb-2902-6fd35527a7b9
A = [0 1 0; 0 0 1; -d₀ -d₁ -d₂]

# ╔═╡ 9268db92-7c0c-11eb-0215-61ffc307720b
evals = eigvals(A)

# ╔═╡ c4c4211a-7c1a-11eb-0d0a-eb9117150022
rc_free(x, p, t) = A * x

# ╔═╡ 80383d7e-7c1b-11eb-3a35-2fb2563e9cc1
tspan = (0.0, 10.0)

# ╔═╡ 0de2bfe2-7c1b-11eb-374c-77552beb35ea
prob = ODEProblem(rc_free, 10*rand(3), tspan)  # Dynamical system with random initial conditions

# ╔═╡ 96e4c58a-7c1b-11eb-21ca-bd6f5a7d9e44
sol = solve(prob)

# ╔═╡ 91af2100-7c1d-11eb-1bb6-63db4fbbacf7
plot(sol, label=["x₁" "x₂" "x₃"])

# ╔═╡ a22040bc-7c1d-11eb-2eef-61d2f2c99ee0
y_free = d₀*sol[1,:] # Output

# ╔═╡ d5147da8-7c1d-11eb-13c7-d982c589f191
plot(sol.t, y_free, label="y")

# ╔═╡ fbecd40e-7c1b-11eb-2335-7db5e56987b2
md"# Controllability and State Controller

A system is called controllable if all states can be driven from any state to the equilibrium. This means, if a system shall be steered from an arbitrary state to the stable equilibrium it has to be guaranteed that the system is controllable. The controllability property can be checked with various methods, e.g. with the Kalman rank condidtion. 

Firstly, the controllability matrix is built with

$\mathcal{C} = [B, A~B, A^2 B] = 
\begin{pmatrix}
0 & 0 & 1 \\
0 & 1 & -d_{2} \\
1 & -d_{2} & -d_{1} + d_{2}^2 
\end{pmatrix}$

and secondly it is checked whether $\mathcal{C}$ has full rank with

$\det(\mathcal{C}) = -1 \neq 0.$


"

# ╔═╡ 61c6f994-7c21-11eb-1741-85f35b125fe3
B = [0; 0; 1]

# ╔═╡ 23feba70-7c21-11eb-32d7-1f31d1a9d9d9
Cntrmat = hcat(B, A*B, A^2*B)

# ╔═╡ 23def686-7c21-11eb-332d-5b4fb312774e
rank(Cntrmat) # rank = 3 -> full rank

# ╔═╡ bcd0ce20-7c1d-11eb-1ef0-0b51461e84f8
md"Therefore, the reachability matrix $\mathcal{R}$ has full rank and the system is controllable. 

The **state feedback** is defined by

$u(t) = -K x(t)$

with contol vector 
$K = 
\begin{pmatrix}
K_{1} & K_{2} & K_{3}
\end{pmatrix}.$

The closed-loop system dynamics is noted as 

$\dot{x}(t) = A ~ x(t) - B ~ u(t) = A x(t) - B ~ K ~ x(t) = (A - B~K) ~ x(t) = A_{CL} ~ x(t)$

with the closed-loop state matrix 

$A_{CL} = A - B~K 
= 
\begin{pmatrix}
0 & 1 & 0 \\
0 & 0 & 1 \\
-d_{0} - K_{1} & -d_{1} - K_{2} & -d_{2} - K_{3}
\end{pmatrix}.$

The filter values $K_{i}$ can be either found manually by calculating the desired Eigenvalues or [Ackermann's formula](https://en.wikipedia.org/wiki/Ackermann%27s_formula) is used to calculate the gain vector directly as

$K = 
\begin{pmatrix}
0 & 0 & \cdots 0 & 1
\end{pmatrix}
~ \mathcal{C}^{-1} ~ p(A)$

in which $p(\lambda)$ is the characteristic polynomial of the desired eigenvalues.

Here, the following Eigenvalues $\Lambda = \left\{-4, -2, -2\right\}$ are desired and thus one holds the characteristic polynomial

$p(\lambda) = (4+\lambda) ~ (2 + \lambda)^2 = \lambda^3 + 8 \lambda^2 + 20 \lambda + 16$

"

# ╔═╡ add0de0a-7c34-11eb-0722-5f5c9cb914e9
p(x) = x^3 + 8*x^2 + 20*x + 16*I # Characteristic polynomial

# ╔═╡ cc9ef362-7c34-11eb-1a51-0dae6d29e81b
K = [0 0 1] * inv(Cntrmat) * p(A)

# ╔═╡ 5cdbcebe-7c35-11eb-0761-41807f408401
Acl = A - B * K # Closed-loop system matrix

# ╔═╡ ecee329e-7c35-11eb-2466-53f2652ee2b0
rc_controlled(x, p, t) = Acl * x

# ╔═╡ f1897e08-7c35-11eb-0674-bd5e4062b116
prob_cntr = ODEProblem(rc_controlled, 10*rand(3), tspan)  # Controlled system dynamics

# ╔═╡ 135336e6-7c36-11eb-015b-3113b7ade37c
sol_cntr = solve(prob_cntr);

# ╔═╡ 27c5bedc-7c36-11eb-27aa-e72f8098f537
plot(sol_cntr, label=["x₁" "x₂" "x₃"])

# ╔═╡ 18543c8a-7c36-11eb-0a93-cfa5247d6e02
y_cntr = d₀*sol_cntr[1,:] # Output

# ╔═╡ 542c90e2-7c36-11eb-2b94-8db400edbcc3
plot(sol_cntr.t, y_cntr, label="y")

# ╔═╡ 52fef61e-7c3a-11eb-0101-2b8ecd70473f
md"## Digital Control Algorithm

In many cases the controller has to be implemented as a digital algorithm. This control algorithm needs the information of all states during the computation and thus the system dynamics has to be calculated numerically.

The dynamical system $(3)$ is integrated as

$x(t_{k+1}) = \exp(A \Delta T) ~ x(t_{k}) + \int_{0}^{\Delta T} \exp\left(A~(\Delta T - \tau)\right) ~ B ~ d\tau ~ u(t_{k})$

where the integral is solved as  

$\int_{0}^{\Delta T} \exp\left(A~(\Delta T - \tau)\right) ~ B ~ d\tau
= A^{-1} ~ [ \exp(A \Delta T) - I ] ~ B$

and finally one holds

$x(t_{k+1}) = \exp(A \Delta T) ~ x(t_{k}) + A^{-1} ~ [ \exp(A \Delta T) - I ] ~ B ~ u(t_{k}) = A_{d} x(t_{k}) + B_{d} u(t_{k}).$

"

# ╔═╡ 3f9972e6-7c41-11eb-05db-675bdb1e22be
ΔT = 0.01

# ╔═╡ fc35d72e-7c40-11eb-38ce-1934b6350215
Ad = exp(A*ΔT)

# ╔═╡ 5b49333c-7c41-11eb-3f22-b748af1cf0b4
Bd = inv(A) * (exp(A*ΔT) - I)*B

# ╔═╡ d7ede7ec-7c42-11eb-199e-5fd2d73a969f
xd0 = rand(3)

# ╔═╡ 543d8cc6-7c43-11eb-0144-2b841267be5c
n_max = 1000;

# ╔═╡ eafb976c-7c42-11eb-019c-174d62afb3c4
xd = zeros(3, n_max)

# ╔═╡ 62ac622a-7c43-11eb-3619-67998b8ac4fe
xd[:,1] = xd0

# ╔═╡ e2f12e92-7c42-11eb-1349-4396c7e95fbf
for i = 1 : n_max-1
	ud = -K * xd[:,i]
	xd[:,i+1] = Ad * xd[:,i] + Bd * ud[1]
end

# ╔═╡ 9ac8e050-7c43-11eb-3d76-97568a6ac0be
plot(xd')

# ╔═╡ 7a8c009a-7c45-11eb-20b7-e7841124b238
md"## Systems with noise

If the system contains a disturance or noise signal that cannot be measured directly then the controller could fail. Consider for example the discrete system 

$x(t_{k+1}) = A_{d} ~ x(t_{k}) + B_{d} ~ u(t_{k}) + w(t_{k})$

with an unknown noise signal $w(t)$. Here, the controller only work with the computed state values $\tilde{x}(t)$ from the nominal system   

$\tilde{x}(t_{k+1}) = A_{d} ~\tilde{x}(t_{k}) + B_{d} ~ u(t_{k})$

without any noise. One notes in the figures below, that the real state values $x(t)$ cannot be steered to the equilibrium (zero). 

"

# ╔═╡ 5289b640-7c44-11eb-01ac-857dbaa83fa9
xd_real = zeros(3, n_max);

# ╔═╡ 526eaa28-7c44-11eb-064a-c5b4a8da62f4
xd_nom = zeros(3, n_max);

# ╔═╡ 51fd0756-7c44-11eb-2c10-cf3888581aa4
xd_real[:,1] = xd_nom[:,1] = xd0

# ╔═╡ 0565da8e-7c45-11eb-047d-9f966c89ba65
for i = 1 : n_max-1
	ud = -K * xd_nom[:,i]
	w = 0.01*rand(3) # noise
	xd_real[:,i+1] = Ad * xd_real[:,i] + Bd * ud[1] + w
	xd_nom[:,i+1] = Ad * xd_nom[:,i] + Bd * ud[1]
end

# ╔═╡ 52a58cb4-7c44-11eb-1279-2f092b232ac9
plot(xd_real', title="Real system")

# ╔═╡ 1ebbd4fa-7c45-11eb-357a-db388a3e85df
plot(xd_nom', title="Nominal system")

# ╔═╡ ac7fde28-7c45-11eb-3113-2330ad967b37
xd[:,1]

# ╔═╡ Cell order:
# ╟─7f35a26e-7c08-11eb-1c76-1d1d7ef04e73
# ╠═7c9fa30a-7c18-11eb-051d-3db5c1bcf4f8
# ╠═7c835e3e-7c18-11eb-31ba-2d5768996c64
# ╠═f94fa7f8-7c18-11eb-1621-ad32e0aa3311
# ╠═3f48146e-7c19-11eb-3450-f7d2340fc4cc
# ╠═5512f8fe-7c19-11eb-3c00-6173a632e2fd
# ╠═661d921c-7c19-11eb-2b55-eb1ad3cd9aa0
# ╟─9382e638-7c0c-11eb-348e-e3bf484b7045
# ╠═9367d0be-7c0c-11eb-2902-6fd35527a7b9
# ╠═93468f12-7c0c-11eb-08f0-b51375661637
# ╠═9268db92-7c0c-11eb-0215-61ffc307720b
# ╠═c4c4211a-7c1a-11eb-0d0a-eb9117150022
# ╠═012acde6-7c1b-11eb-3cb8-4b5dd45c4194
# ╠═80383d7e-7c1b-11eb-3a35-2fb2563e9cc1
# ╠═0de2bfe2-7c1b-11eb-374c-77552beb35ea
# ╠═96e4c58a-7c1b-11eb-21ca-bd6f5a7d9e44
# ╠═8e2e24ea-7c1d-11eb-2a8e-1592f75ae0b3
# ╠═91af2100-7c1d-11eb-1bb6-63db4fbbacf7
# ╠═a22040bc-7c1d-11eb-2eef-61d2f2c99ee0
# ╠═d5147da8-7c1d-11eb-13c7-d982c589f191
# ╟─fbecd40e-7c1b-11eb-2335-7db5e56987b2
# ╠═61c6f994-7c21-11eb-1741-85f35b125fe3
# ╠═23feba70-7c21-11eb-32d7-1f31d1a9d9d9
# ╠═23def686-7c21-11eb-332d-5b4fb312774e
# ╟─bcd0ce20-7c1d-11eb-1ef0-0b51461e84f8
# ╠═add0de0a-7c34-11eb-0722-5f5c9cb914e9
# ╠═cc9ef362-7c34-11eb-1a51-0dae6d29e81b
# ╠═5cdbcebe-7c35-11eb-0761-41807f408401
# ╠═ecee329e-7c35-11eb-2466-53f2652ee2b0
# ╠═f1897e08-7c35-11eb-0674-bd5e4062b116
# ╠═135336e6-7c36-11eb-015b-3113b7ade37c
# ╠═27c5bedc-7c36-11eb-27aa-e72f8098f537
# ╠═18543c8a-7c36-11eb-0a93-cfa5247d6e02
# ╠═542c90e2-7c36-11eb-2b94-8db400edbcc3
# ╟─52fef61e-7c3a-11eb-0101-2b8ecd70473f
# ╠═3f9972e6-7c41-11eb-05db-675bdb1e22be
# ╠═fc35d72e-7c40-11eb-38ce-1934b6350215
# ╠═5b49333c-7c41-11eb-3f22-b748af1cf0b4
# ╠═d7ede7ec-7c42-11eb-199e-5fd2d73a969f
# ╠═543d8cc6-7c43-11eb-0144-2b841267be5c
# ╠═eafb976c-7c42-11eb-019c-174d62afb3c4
# ╠═62ac622a-7c43-11eb-3619-67998b8ac4fe
# ╠═e2f12e92-7c42-11eb-1349-4396c7e95fbf
# ╠═9ac8e050-7c43-11eb-3d76-97568a6ac0be
# ╟─7a8c009a-7c45-11eb-20b7-e7841124b238
# ╠═5289b640-7c44-11eb-01ac-857dbaa83fa9
# ╠═526eaa28-7c44-11eb-064a-c5b4a8da62f4
# ╠═51fd0756-7c44-11eb-2c10-cf3888581aa4
# ╠═0565da8e-7c45-11eb-047d-9f966c89ba65
# ╠═52a58cb4-7c44-11eb-1279-2f092b232ac9
# ╠═1ebbd4fa-7c45-11eb-357a-db388a3e85df
# ╠═ac7fde28-7c45-11eb-3113-2330ad967b37
