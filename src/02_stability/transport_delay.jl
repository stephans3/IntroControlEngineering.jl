### A Pluto.jl notebook ###
# v0.12.12

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 9d8d3e38-5580-11eb-2526-550f43044d61
using ControlSystems

# ╔═╡ 71e4cfa8-5585-11eb-2f0a-eba8f744d3d4
md"# Example: Transport systems

Transport systems are typical candidates for systems. They can be described by the one-dimensional partial differential equation (or simply transport equation)

$\frac{\partial x(t,z)}{\partial t} + \nu  \frac{\partial x(t,z)}{\partial x} = u(t,x) \quad \text{for} ~ z \in (0,L), t>0$ 

with velocity $\nu > 0$, left boundary condition $x(t, 0) = u_{0}(t)$ for $t > 0$ and the initial data 

$x(0, z) = x_{z}(0) \quad \text{for} ~ z \in [0,L].$

Function $u(t,x)$ is the distributed control input along the complete length and function $u_{0}(t)$ is the boundary input on the left side ($x=0$).

For simplicity, it is assumed that control input $u(t,x)$ is constant along the whole length.

## Solution in Laplace-domain

Further the output is assumed at the right side of the system ($x=L$) as $y(t) = x(t,L)$. A Laplace transform of the PDE and further calculations lead to

$\hat{y}(s) = \exp(-s ~ T_{t}) ~ \hat{u}_{0}(s) + \frac{1}{s} \left[1 - \exp(-s ~ T_{t}) \right] ~ \hat{u}(s)$

with $T_{t} = \frac{L}{\nu}$.


## Test cases

Two scenarios are considered next:
- Fluid flow in a pipe with single boundary input $u_{0}(t) > 0$, thus $u(t,x) = 0$  and
- Conveyor system with distributed input $u(t) = u(t,x) > 0$, thus $u_{0}(t) = 0$.

"

# ╔═╡ 71cb4ec0-5585-11eb-17b8-3943c40b9006
md"## Fluid flow in a pipe

The transfer function is noted as

$\hat{g}_{0} = \frac{\hat{y}(s)}{\hat{u}_{0}(s)} = \exp(-s ~ T_{t}) ~  \quad \text{with} ~ T_{t} = \frac{L}{\nu}.$
"

# ╔═╡ a85a65de-5580-11eb-14cd-b1412d15a4eb
L = 1.0 # Length of pipe

# ╔═╡ a95e5f62-5580-11eb-2047-ff24cf1662e9
ν = 0.2 # Flow velocity [m/s]

# ╔═╡ 0ddc8b80-5581-11eb-2819-696169f8fc52
Tt = L/ν # Delay

# ╔═╡ 153fc7a2-5581-11eb-19b5-87ee0368c939
g0 = delay(Tt) # Ideal delay

# ╔═╡ 72004dd2-558a-11eb-1830-79352529ab3e
md"The ideal delay is approximated by

$\hat{g}_{approx} = \frac{a}{s + a} ~ \exp(-s ~ T_{t})$

with $a > 0$.
"

# ╔═╡ 26a03de4-5584-11eb-37f1-75bb7ba6648d
@bind a html"<input type='range' start='1' step='1'>"

# ╔═╡ 2912264a-558f-11eb-12b6-61891524f28f
md"Coefficient a=$a"

# ╔═╡ 1d00f396-5584-11eb-3d5d-9587b26069d0
g_approx = tf(a, [1, a])*g0

# ╔═╡ 15796bcc-558b-11eb-0239-7728b45c7e23
md"### Step response"

# ╔═╡ 85c1bc9c-5581-11eb-2311-ab784036ad88
Tf = 10 # Final simulation time

# ╔═╡ 226347f0-5582-11eb-393b-c13858c1b6e4
stepplot(g_approx, Tf)

# ╔═╡ a02fad28-558b-11eb-01f1-4f49f23b9c69
md"### Bode plot"

# ╔═╡ dac1a49e-5583-11eb-17f9-1f6fa47e8c56
ω = 10^(-2) :  10^(-2) :  10^(2) # Frequencies for Bode and Nyquist plot

# ╔═╡ a6ef6bbc-558b-11eb-0fd8-8f20ca665a4c
bodeplot(g_approx, ω, legend=false)

# ╔═╡ 051aca9c-558c-11eb-23b1-8d7f5500a9d7
setPlotScale("dB")

# ╔═╡ 163f2a20-558c-11eb-0044-35193342b447
md"### Nyquist plot"

# ╔═╡ 213440c8-558c-11eb-228b-e55457c147d7
nyquistplot(g_approx, ω, legend=false)

# ╔═╡ 49050f24-558c-11eb-1f3b-e18b685e0c02
md"## Conveyor system
The transfer function is noted as

$\hat{g}_{0} = \frac{\hat{y}(s)}{\hat{u}_{0}(s)} = \frac{1}{s} \left[1 - \exp(-s ~ T_{t})\right]  ~  \quad \text{with} ~ T_{t} = \frac{L}{\nu}.$

Values for $L$ and $\nu$ are the same as above.
"

# ╔═╡ 35175d56-5581-11eb-0a5c-03471278b153
g1 = tf(1, [1, 0]) - tf(1, [1, 0])*delay(Tt)

# ╔═╡ 3ce70d04-558d-11eb-39a9-95328cba1821
md"### Step plot"

# ╔═╡ 505baaac-558d-11eb-3e44-e7754436d73c
stepplot(g1, Tf)

# ╔═╡ 596acf92-558d-11eb-0846-ddecd6627774
md"### Bode plot"

# ╔═╡ c3e371fe-5582-11eb-2084-55dcb05584f2
bodeplot(g1, ω, legend=false)

# ╔═╡ 63c51f38-558d-11eb-22c9-f9f0f2566db7
md"### Nyquist plot"

# ╔═╡ fccb111a-5583-11eb-0cdf-e1cd7d9b2e89
nyquistplot(100.0*g1, ω)

# ╔═╡ Cell order:
# ╟─71e4cfa8-5585-11eb-2f0a-eba8f744d3d4
# ╟─71cb4ec0-5585-11eb-17b8-3943c40b9006
# ╠═9d8d3e38-5580-11eb-2526-550f43044d61
# ╠═a85a65de-5580-11eb-14cd-b1412d15a4eb
# ╠═a95e5f62-5580-11eb-2047-ff24cf1662e9
# ╠═0ddc8b80-5581-11eb-2819-696169f8fc52
# ╠═153fc7a2-5581-11eb-19b5-87ee0368c939
# ╟─72004dd2-558a-11eb-1830-79352529ab3e
# ╠═26a03de4-5584-11eb-37f1-75bb7ba6648d
# ╠═2912264a-558f-11eb-12b6-61891524f28f
# ╠═1d00f396-5584-11eb-3d5d-9587b26069d0
# ╟─15796bcc-558b-11eb-0239-7728b45c7e23
# ╠═85c1bc9c-5581-11eb-2311-ab784036ad88
# ╠═226347f0-5582-11eb-393b-c13858c1b6e4
# ╟─a02fad28-558b-11eb-01f1-4f49f23b9c69
# ╠═dac1a49e-5583-11eb-17f9-1f6fa47e8c56
# ╠═a6ef6bbc-558b-11eb-0fd8-8f20ca665a4c
# ╠═051aca9c-558c-11eb-23b1-8d7f5500a9d7
# ╠═163f2a20-558c-11eb-0044-35193342b447
# ╠═213440c8-558c-11eb-228b-e55457c147d7
# ╠═49050f24-558c-11eb-1f3b-e18b685e0c02
# ╠═35175d56-5581-11eb-0a5c-03471278b153
# ╠═3ce70d04-558d-11eb-39a9-95328cba1821
# ╠═505baaac-558d-11eb-3e44-e7754436d73c
# ╠═596acf92-558d-11eb-0846-ddecd6627774
# ╠═c3e371fe-5582-11eb-2084-55dcb05584f2
# ╠═63c51f38-558d-11eb-22c9-f9f0f2566db7
# ╠═fccb111a-5583-11eb-0cdf-e1cd7d9b2e89
