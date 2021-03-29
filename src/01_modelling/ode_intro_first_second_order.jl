### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 16b42ff6-7ab1-11eb-030f-975301cead26
using Plots

# ╔═╡ 2a8d86bc-7ab6-11eb-3230-d776bd8e2b60
using LinearAlgebra

# ╔═╡ 5a1c01d0-7aaf-11eb-2c49-b748d4a61780
md"# Ordinary Differential Equations


## First order system
Consider the first order differential equation 

$\dot{y}(t) - \alpha ~ y(t) ~=~ 0$
with the initial values 

$y(0) = y_{0} \text{.}$

The solution is found by integration as

$\int_{y_{0}}^{y} \frac{d\tilde{y}}{\tilde{y}} = \int_{t_{0}}^{t} \alpha ~ d\tilde{t}$

and following 

$\ln({y}) - \ln ({y_{0}}) =  \alpha ~ (t - t_{0}) \text{.}$

Therefore, the solution is given as 

$y(t) = y_{0} ~ \exp\left( \alpha ~ (t - t_{0}) \right)$ 

and simplified with $ t_{0} = 0 $ to

$y(t) = y_{0} ~ \operatorname{e}^{\alpha t} \text{.}$ 

If $\alpha < 0$ then $\lim\limits_{t \rightarrow \infty} y(t) = 0$. Otherwise for $\alpha > 0$ then $\lim\limits_{t \rightarrow \infty} y(t) = \infty$.

### Numerical solution

Firstly, the constants and the time interval with its time steps are specified.
"

# ╔═╡ 06254fba-7ab0-11eb-0775-6de23f7fac3c
y₀ = 2.0; # Initial value

# ╔═╡ 90413b1c-7ab0-11eb-0d44-cf79424ccf9a
α = -1.2;

# ╔═╡ 9726b036-7ab0-11eb-2381-595fa350110b
t = 0.0 : 0.1 : 10; # Time range

# ╔═╡ e725252e-7ab0-11eb-0c80-e30e20f1d0c3
md"Secondly, the solution of the differential equation $y(t)$ is defined."

# ╔═╡ f8dbdcae-7ab0-11eb-2d02-5b6782abd3b5
y = y₀ * exp.(α*t); # Analytical solution

# ╔═╡ 023ab072-7ab1-11eb-2df6-b1413ddf3e1b
md"Finally, the graph of the solution is plotted."

# ╔═╡ 1ec0fb66-7ab1-11eb-0a6d-4102843df0a7
plot(t, y, label="y", title="First order system")

# ╔═╡ 3fe74872-7ab1-11eb-3413-c186287b3655
md"## Second order systems

Consider the second order differential equation

$\ddot{y}(t) = \alpha y(t) + \beta \dot{y}(t)$

with the initial values

$y(0) = y_{01} \quad \text{and} \quad \dot{y}(0) = y_{02} \text{.}$

### Undamped oscillation

If $\beta = 0$, then one holds

$\ddot{y}(t) = \alpha ~ y(t)$

which can be solved with the approach

$y(t) = A ~ \sin(\omega t) + B ~ \cos(\omega t).$

From the second derivative 

$\ddot{y}(t) = -A~\omega^2 ~ \sin(\omega t) - B~ \omega^2 ~ \cos(\omega t) = -\omega^2 (A ~ \sin(\omega t) + B ~ \cos(\omega t)) = - \omega^2 ~ y(t)$

follows  $\alpha = \omega^2$.

The constants $A$ and $B$ are specified by

$y(0) = B = y_{0,1} \quad \text{and} \quad \dot{y}(0) = A \omega =  y_{0,2}$

and thus $A = \frac{y_{02}}{\omega} = \frac{y_{02}}{\sqrt{\alpha}}$.

### Numerical solution
The initial values and the constants are defined as
"

# ╔═╡ 6ecd9532-7ab2-11eb-2191-33ca1f7f0431
y01 = 2.; 	# Initial value, e.g initial position

# ╔═╡ 86c4561c-7ab2-11eb-0b16-a57f8e43ce24
y02 = 1.2; 	# Initial value, e.g initial velocity

# ╔═╡ a29f6b60-7ab2-11eb-0b70-bf6e94017fba
α₂ = 3.1; 	# Parameter

# ╔═╡ 05da3926-7ab3-11eb-3535-add547916ab1
ω = sqrt(α₂); # Natural frequency or Eigenfrequency

# ╔═╡ 091d4c36-7ab3-11eb-3d3d-53d83d5e76db
A = y02/ω;

# ╔═╡ 0c2e2efe-7ab3-11eb-0b62-6fb8ec3f47f8
B = y01;

# ╔═╡ df025c66-7ab2-11eb-2b6e-7b977f23dd88
md"Thus, the solution is given as"

# ╔═╡ ecb2e524-7ab2-11eb-33dd-21cbfc722299
y₂ = A * sin.(ω*t) + B * cos.(ω*t); # Analytical solution

# ╔═╡ 4e44dc16-7ab3-11eb-373b-ddc0b3f7bf4c
plot(t, y₂, label="y", title="Undamped oscillation")

# ╔═╡ 4dc55c52-7ab3-11eb-35bb-5750adac8e8c
md"### Damped oscillation

If $ \beta \neq 0 $, then the second order system is transformed via

$x_{1}(t) = y(t) \quad \text{and} \quad x_{2}(t) = \dot{y}(t).$

to the first order ODE

$\begin{pmatrix}
\dot{x}_{1}(t) \\
\dot{x}_{2}(t)
\end{pmatrix}
= 
\begin{pmatrix}
0 & 1 \\
\alpha & \beta
\end{pmatrix}
\begin{pmatrix}
{x}_{1}(t) \\
{x}_{2}(t)
\end{pmatrix}
= A ~ x(t) \text{.}$


#### Stability

The first order system is stable if all eigenvalues of matrix $A$ are smaller than zero. This is proved by

$\det(\lambda I - A) = \det
\begin{pmatrix}
\lambda & -1 \\
-\alpha & \lambda - \beta
\end{pmatrix}
= \lambda ~ (\lambda - \beta) - \alpha = \lambda^2 - \beta~\lambda  - \alpha = 0$
and further

$\lambda_{1,2} = \frac{\beta}{2} \pm \frac{1}{2} \sqrt{\beta^2 + 4 \alpha} \text{.}$ 

The eigenvalues $\lambda_{1,2}$ are smaller than zero if $\beta < 0$ and

$\sqrt{\beta^2 + 4 \alpha} < \beta$  and thus $\alpha < 0$.
"

# ╔═╡ 5e0f9ddc-7ab5-11eb-275a-e1d8bef0dcc9
α₃ = -0.25; # Parameters

# ╔═╡ 5df31022-7ab5-11eb-00fc-fb69baa05a10
β₃ = -2.0;

# ╔═╡ 5d5255a6-7ab5-11eb-34e2-bdfefdeaf74f
A₃ = [0 1; α₃ β₃]

# ╔═╡ a25d9386-7ab5-11eb-0094-03f733b1d020
λ₁ = β₃/2 + sqrt(β₃^2 + 4*α₃)/2 # Eigenvalues

# ╔═╡ ca3f15e6-7ab5-11eb-2a20-95812bef0ef1
λ₂ = β₃/2 - sqrt(β₃^2 + 4*α₃)/2

# ╔═╡ d4e32494-7ab5-11eb-380a-d98f4e5939e4
md"
#### Solution

The solution of the first order system is calculated with the eigenvectors $v_{1}$, $v_{2}$ and constants $C_{1}$, $C_{2}$ as

$x(t) = C_{1} ~ v_{1} ~ \exp(\lambda_{1} ~ t) + C_{2} ~ v_{2} ~ \exp(\lambda_{2} ~ t) \text{.}$

The eigenvectors $v_{1}$, $v_{2}$ are calculated with the eigenvalues as 
$(\lambda_{1} I - A) v_{1} = 0$ and $(\lambda_{2} I - A) v_{2} = 0$.

The constants $C_{1}$, $C_{2}$ are calculated via the initial values $y(0) = x_{1}(0) = y_{0,1}$ and $\dot{y}(0) = x_{2}(0) = y_{0,2}$ as


$\begin{pmatrix}
{x}_{1}(0) \\
{x}_{2}(0)
\end{pmatrix}
=
\begin{pmatrix}
y_{01} \\
y_{02}
\end{pmatrix}
= 
[v_{1}, v_{2}]
\begin{pmatrix}
C_{1} \\
C_{2}
\end{pmatrix}
\quad \text{and thus}$

and thus

$\begin{pmatrix}
C_{1} \\
C_{2}
\end{pmatrix}
= 
[v_{1}, v_{2}]^{-1}
\begin{pmatrix}
y_{01} \\
y_{02}
\end{pmatrix} \text{.}$

### Numerical solution

The initial values $y_{01}$, $y_{02}$ and the constants $\alpha$, $\beta$ are given as above. The eigenvalues and eigenvectors of matrix $A$ are calculated as 
"

# ╔═╡ 2a70f2fe-7ab6-11eb-21c8-f78ce47f59df
eival = eigvals(A₃) # Eigenvalues of A

# ╔═╡ 3448b096-7ab6-11eb-15c3-89fb20f70db0
eivec = eigvecs(A₃) # Eigenvectors of A

# ╔═╡ 5a0f4a38-7ab6-11eb-1c59-b35c939d5c62
md"The constants $C_{1}$, $C_{2}$ are calculated as discussed theoretically and the solution $y(t)$ is calculated and plotted."

# ╔═╡ 631b5f4a-7ab6-11eb-25a6-9528020d35de
C = inv(eivec)*[y01, y02];

# ╔═╡ 73b2f57a-7ab6-11eb-2d12-19d4916f8f08
t₃ = 0.0 : 0.1 : 50.0; # Time range

# ╔═╡ 77ccf7f0-7ab6-11eb-3445-b1626fae1ce6
y₃ = eivec * [C[1]*exp.(eival[1]*t₃),C[2]*exp.(eival[2]*t₃)]; # Analytical solution

# ╔═╡ 82abefaa-7ab6-11eb-193c-eba3cbfd0402
plot(t₃, y₃, label=["x1" "x2"], title="Damped oscillation")

# ╔═╡ Cell order:
# ╟─5a1c01d0-7aaf-11eb-2c49-b748d4a61780
# ╠═06254fba-7ab0-11eb-0775-6de23f7fac3c
# ╠═90413b1c-7ab0-11eb-0d44-cf79424ccf9a
# ╠═9726b036-7ab0-11eb-2381-595fa350110b
# ╟─e725252e-7ab0-11eb-0c80-e30e20f1d0c3
# ╠═f8dbdcae-7ab0-11eb-2d02-5b6782abd3b5
# ╟─023ab072-7ab1-11eb-2df6-b1413ddf3e1b
# ╠═16b42ff6-7ab1-11eb-030f-975301cead26
# ╠═1ec0fb66-7ab1-11eb-0a6d-4102843df0a7
# ╟─3fe74872-7ab1-11eb-3413-c186287b3655
# ╠═6ecd9532-7ab2-11eb-2191-33ca1f7f0431
# ╠═86c4561c-7ab2-11eb-0b16-a57f8e43ce24
# ╠═a29f6b60-7ab2-11eb-0b70-bf6e94017fba
# ╠═05da3926-7ab3-11eb-3535-add547916ab1
# ╠═091d4c36-7ab3-11eb-3d3d-53d83d5e76db
# ╠═0c2e2efe-7ab3-11eb-0b62-6fb8ec3f47f8
# ╟─df025c66-7ab2-11eb-2b6e-7b977f23dd88
# ╠═ecb2e524-7ab2-11eb-33dd-21cbfc722299
# ╠═4e44dc16-7ab3-11eb-373b-ddc0b3f7bf4c
# ╟─4dc55c52-7ab3-11eb-35bb-5750adac8e8c
# ╠═5e0f9ddc-7ab5-11eb-275a-e1d8bef0dcc9
# ╠═5df31022-7ab5-11eb-00fc-fb69baa05a10
# ╠═5d5255a6-7ab5-11eb-34e2-bdfefdeaf74f
# ╠═a25d9386-7ab5-11eb-0094-03f733b1d020
# ╠═ca3f15e6-7ab5-11eb-2a20-95812bef0ef1
# ╟─d4e32494-7ab5-11eb-380a-d98f4e5939e4
# ╠═2a8d86bc-7ab6-11eb-3230-d776bd8e2b60
# ╠═2a70f2fe-7ab6-11eb-21c8-f78ce47f59df
# ╠═3448b096-7ab6-11eb-15c3-89fb20f70db0
# ╟─5a0f4a38-7ab6-11eb-1c59-b35c939d5c62
# ╠═631b5f4a-7ab6-11eb-25a6-9528020d35de
# ╠═73b2f57a-7ab6-11eb-2d12-19d4916f8f08
# ╠═77ccf7f0-7ab6-11eb-3445-b1626fae1ce6
# ╠═82abefaa-7ab6-11eb-193c-eba3cbfd0402
