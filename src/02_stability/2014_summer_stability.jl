### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 923cffb6-8e5b-11eb-281a-6dbdefd9fdc1
using ControlSystems

# ╔═╡ ec2ed6ac-8e5e-11eb-3cae-3ffebf0cd801
using LinearAlgebra

# ╔═╡ 98368060-8e5e-11eb-2df5-e16684f13eed
using OrdinaryDiffEq

# ╔═╡ fec26aae-8e5e-11eb-35f0-3da183c3bfa3
using Plots

# ╔═╡ 590d7aac-8e56-11eb-3bc1-9f582de65c84
md"# Example: BIBO stability

Consider the dynamical system 

$\dddot{y}(t) + 4 ~ \ddot{y}(t) + 5 ~ \dot{y}(t) = u(t) \tag{1}$

with initial value $y(0) = y_{0}$.

## Transfer function

It is assumed that the initial value $y(0) = 0$ and the Laplace transform of ODE $(1)$ is noted as

$s^3 Y(s) + 4~s^2~Y(s) + 5~s~Y(s) = U(s)$

and the transfer function is derived as

$G(s) = \frac{Y(s)}{U(s)} = \frac{1}{s^3 + 4 s^2 + 5s} = \frac{1}{(s^2 + 4 s + 5) s}.\tag{2}$

The poles of transfer function $(2)$ are found with the characteristic polynomial 

$p(s) = (s^2 + 4 s + 5) s = 0.$

One pole is simply found as $s_{1}=0$ and the other two poles are found with

$s_{2,3} = \frac{-4}{2} \pm \frac{1}{2} \sqrt{16 - 4*5}$
"

# ╔═╡ 31f5bdc2-8e58-11eb-2ceb-4dfdf031fd43
a₁ = 5

# ╔═╡ 45820b86-8e58-11eb-345d-61f446ff4a59
a₂ = 4

# ╔═╡ 4ae12c56-8e58-11eb-2127-378e043f40e2
s₂ = -a₂/2 - (1/2) * sqrt(Complex(a₂^2 - 4*a₁))

# ╔═╡ 6a33bdda-8e58-11eb-2a74-87321a371926
s₃ = -a₂/2 + (1/2) * sqrt(Complex(a₂^2 - 4*a₁))

# ╔═╡ 96d8a636-8e5b-11eb-2fed-173b68b5063b
G = tf(1, [1, 4, 5, 0])

# ╔═╡ 92008d16-8e5d-11eb-3d7e-57ddaba4181f
Tf = 10.0; # Simulation time

# ╔═╡ 8ce0218e-8e5d-11eb-2b0f-fbecc03cf8f4
stepplot(G, Tf)

# ╔═╡ d803f906-8e5d-11eb-0248-3bb79cb74813
impulseplot(G,Tf)

# ╔═╡ bdade814-8e58-11eb-1dea-b1aac45cddda
md"Due to the one pole at $s_{1}=0$ the system is called marginally stable.

## State-space representation

Now, the stability is proved in the state-space representation. The states of dynamical systems $(1)$ are introduced as

$x_{1}(t) = y(t), \quad x_{2}(t) = \dot{y}(t) = \dot{x}_{1}(t), \quad x_{3}(t) = \ddot{y}(t) = \dot{x}_{2}(t)$

and the system of first-order ODEs is noted as

$\begin{pmatrix}
\dot{x}_{1}(t) \\
\dot{x}_{2}(t) \\
\dot{x}_{3}(t) 
\end{pmatrix} =
\begin{pmatrix}
0 & 1 & 0 \\
0 & 0 & 1 \\
0 & -5 & -4 
\end{pmatrix} ~
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
u(t)$

with output

$y(t) = (1,~ 0,~ 0) ~ x(t) = x_{1}(t).$

The stability of dynamical systems in the state-space are proved with Eigenvalues which are found with

$\det(\lambda I - A) = \det
\begin{pmatrix}
\lambda & -1 & 0 \\
0 & \lambda & -1 \\
0 & 5 & \lambda+4 
\end{pmatrix} =\lambda^2 (\lambda+4) - 5 ~ (-1) ~ \lambda = 0$

and 

$\lambda^2 (\lambda+4) - 5 ~ (-1) ~ \lambda = \lambda^3 + 4~\lambda^2 + 5 \lambda = \lambda (\lambda^2 + 4~\lambda + 5) = 0.$

The characteristic polynomial of the Eigenvalues 

$p(\lambda) = \lambda (\lambda^2 + 4~\lambda + 5) = 0$

is equal to characteristic polynomial of the poles and thus the Eigenvalues are 

$\lambda = \{0, -2 \pm j\}.$

"

# ╔═╡ f616f770-8e5d-11eb-2166-d3a317813d93
A = [0 1 0; 0 0 1; 0 -a₁ -a₂] # System matrix

# ╔═╡ 3f5449d4-8e5f-11eb-23d0-13b701945c9e
md"**Eigenvalues**"

# ╔═╡ da5eff60-8e5e-11eb-16aa-f5a1b0902244
evals = eigvals(A)

# ╔═╡ 7ea87eec-8e5e-11eb-2c6e-1b8b79915548
B = [0, 0, 1] # Input vector

# ╔═╡ 13dc8e66-8e5e-11eb-377e-a3233b3d2454
function dyn_sys!(dx, x, p, t)
	
	u = 1.0 # Step response
	#t > 0.0 ? u=0.0 : u=1 # Impulse response
	
	dx .= A * x + B * u
	
end

# ╔═╡ ab03bf60-8e5e-11eb-11ea-492a509e8cf0
x0 = [0.0, 0.0, 0.0]

# ╔═╡ b2aa6b58-8e5e-11eb-07e7-bb72b6fcf984
tspan = (0.0, Tf)

# ╔═╡ 9f7da4d2-8e5e-11eb-3058-45a15ac77204
prob = ODEProblem(dyn_sys!, x0, tspan)

# ╔═╡ c78c3d8c-8e5e-11eb-1f89-ffffdbf52483
sol = solve(prob, Tsit5());

# ╔═╡ 01d5035a-8e5f-11eb-0341-833ea46b1047
plot(sol, legend=:topleft)

# ╔═╡ Cell order:
# ╟─590d7aac-8e56-11eb-3bc1-9f582de65c84
# ╠═31f5bdc2-8e58-11eb-2ceb-4dfdf031fd43
# ╠═45820b86-8e58-11eb-345d-61f446ff4a59
# ╠═4ae12c56-8e58-11eb-2127-378e043f40e2
# ╠═6a33bdda-8e58-11eb-2a74-87321a371926
# ╠═923cffb6-8e5b-11eb-281a-6dbdefd9fdc1
# ╠═96d8a636-8e5b-11eb-2fed-173b68b5063b
# ╠═92008d16-8e5d-11eb-3d7e-57ddaba4181f
# ╠═8ce0218e-8e5d-11eb-2b0f-fbecc03cf8f4
# ╠═d803f906-8e5d-11eb-0248-3bb79cb74813
# ╟─bdade814-8e58-11eb-1dea-b1aac45cddda
# ╠═f616f770-8e5d-11eb-2166-d3a317813d93
# ╠═ec2ed6ac-8e5e-11eb-3cae-3ffebf0cd801
# ╟─3f5449d4-8e5f-11eb-23d0-13b701945c9e
# ╠═da5eff60-8e5e-11eb-16aa-f5a1b0902244
# ╠═7ea87eec-8e5e-11eb-2c6e-1b8b79915548
# ╠═13dc8e66-8e5e-11eb-377e-a3233b3d2454
# ╠═98368060-8e5e-11eb-2df5-e16684f13eed
# ╠═ab03bf60-8e5e-11eb-11ea-492a509e8cf0
# ╠═b2aa6b58-8e5e-11eb-07e7-bb72b6fcf984
# ╠═9f7da4d2-8e5e-11eb-3058-45a15ac77204
# ╠═c78c3d8c-8e5e-11eb-1f89-ffffdbf52483
# ╠═fec26aae-8e5e-11eb-35f0-3da183c3bfa3
# ╠═01d5035a-8e5f-11eb-0341-833ea46b1047
