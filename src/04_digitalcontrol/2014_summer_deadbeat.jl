### A Pluto.jl notebook ###
# v0.12.12

using Markdown
using InteractiveUtils

# ╔═╡ bdfe6152-5b3b-11eb-2ab8-9507d5ad5a60
using Plots

# ╔═╡ 1b3761c0-5b34-11eb-3a51-915a1fb66273
md"# Dead-beat controller design: Introduction

A general plant is described as discrete transfer function by

$G_{p}^{*}(z) = \frac{B(z)}{A(z)} ~ z^{-d}$

with 

$A(z) = a_{0} + a_{1} z^{-1} +  a_{2} z^{-2} + \cdots a_{m} z^{-m} \text{,}$

$B(z) = b_{0} + b_{1} z^{-1} +  b_{2} z^{-2} + \cdots b_{n} z^{-n}$

and discrete time-delay $z^{-d}$.

The Deadbeat controller for such a general plant is defined as 

$G_{c}^{*}(z) = \frac{A(z)}{B(1) - z^{-d} ~ B(z)}$

and many authors also note it in terms 

$G_{c}^{*}(z) = \frac{q_{0} ~ A(z)}{1 - z^{-d} ~ q_{0} ~  B(z)}$

with 

$q_{0} = \frac{1}{b_{0} + b_{1} + \cdots b_{n}}.$ 

**Both formulas are equivalent!**

The open-loop system of the general plant and the Deadbeat controller is derived as 

$G_{OL}^{*}(z) = G_{p}^{*}(z) ~ G_{c}^{*}(z) = \frac{z^{-d} ~ B(z)}{B(1)- z^{-d} ~ B(z)}$

and consequently the closed-loop system is noted as 

$G_{CL}^{*}(z) = \frac{G_{ol}^{*}(z)}{1 + G_{ol}^{*}(z)} = \frac{B(z)}{B(1)}.$

**Remark:** The Deadbeat controller should only be applied to stable plants because the possible instable poles (from A(z)) are only eliminated in the calculation - not in the real system!

"

# ╔═╡ 78a8935a-5b2b-11eb-0bc0-157b068a2335
md"## Example

Consider the transfer function of a plant as

$G_{p}(s) = \frac{1}{5 s^2 + s + 1}.$

This plant is stable because its poles are 

$s_{1,2} = -0.1 \pm j \frac{\sqrt{19}}{10}.$

## Z-transform

The continous plant is tranformed to the Z-domain with

$G_{p}^{*}(z) = \frac{z - 1}{z} ~ \mathcal{Z}\left\{ \left. \mathcal{L}^{-1}\left\{ \frac{G_{p}(s)}{s} \right\} \right\rvert_{t=n \Delta T} \right\}.$

Here, a sampling time of $\Delta T = 0.5$ seconds is used to determine 

$G_{p}^{*}(z) = \frac{0.0241 ~ z^{-1} + 0.0233 ~ z^{-2}}{1 - 1.8575 ~ z^{-1} + 0.9048 ~ z^{-2}}.$

The Z-transformed discrete system is stable because the discrete poles are inside the unit circle, namely

$z_{1,2} = 0.9287 \pm j0.2055 \quad \text{and} \quad \lvert z_{1,2} \rvert = 0.9512 < 1.$

The poles can also be simply found by evaluating the formula

$z_{1,2} = e^{s_{1,2} \Delta T}.$
"

# ╔═╡ e4956e92-5b32-11eb-0962-b36f35ccd3bb
s0 = -0.1 + (sqrt(19)/10)im

# ╔═╡ 08b76d6a-5b34-11eb-3489-55eedc98d28a
ΔT = 0.5

# ╔═╡ f1906638-5b32-11eb-2006-870fdd2b0f43
z0 = exp(s0*ΔT)

# ╔═╡ 337dfa6a-5b33-11eb-2586-7d378fb12025
abs(z0) < 1

# ╔═╡ 1b1e2eda-5b34-11eb-2f93-353e2d242f61
md"## Application of dead-beat controller

The discrete transfer function given as

$G_{p}^{*}(z) = \frac{0.0241 ~ z^{-1} + 0.0233 ~ z^{-2}}{1 - 1.8575 ~ z^{-1} + 0.9048 ~ z^{-2}} = \frac{B(z)}{A(z)}$

is used to derive the dead-beat controller as

$G_{c}^{*}(z) = \frac{1 - 1.8575 ~ z^{-1}  + 0.9048 ~ z^{-2}}{0.0474 - 0.0241 ~ z^{-1} - 0.0233 ~ z^{-2}}.$

Further the closed-loop system is calculated by

$G_{CL}^{*}(z) = \frac{B(z)}{B(1)} = \frac{0.0241 ~ z^{-1} + 0.0233 ~ z^{-2}}{0.0474} = 0.5084 z^{-1} + 0.4916 z^{-2}.$

## Simulation of feedback system

Consider the discrete reference $r_{n}$ and output value $y_{n}$ of closed-loop system

$G_{CL}^{*}(z) = \frac{Y(z)}{R(z)} =  0.5084 z^{-1} + 0.4916 z^{-2}.$ 

This transfer function is transformed to the discrete-time algorithm with

$Y(z) = 0.5084 z^{-1} R(z) + 0.4916 z^{-2} R(z)$

and

$y_{n} = 0.5084 ~ r_{n-1} + 0.4916  ~ r_{n-2}.$
"

# ╔═╡ bbdf95c4-5b3a-11eb-0513-71595f7699b8
N = 5 # Number of time steps

# ╔═╡ bbc45868-5b3a-11eb-3653-c7b015841dd8
y = zeros(N) # Output

# ╔═╡ bba3ba54-5b3a-11eb-212f-25d5ec2e3512
r = ones(N) # Reference

# ╔═╡ bb807cba-5b3a-11eb-2d5f-b160c1de903d
for i = 1 : N
	if i == 1
		y[i] = 0
	elseif i == 2
		y[i] = 0.5084 * r[i-1]
	else
		y[i] = 0.5084 * r[i-1] + 0.4916*r[i-2]
	end
end

# ╔═╡ bb61b212-5b3a-11eb-13fd-4727abfac6c4
plot([0,1,2,3,4], y, legend=false)

# ╔═╡ 188ebe40-5b3c-11eb-11c0-03ab24de191c
md"### Dead-beat time

The controller needs two steps to reach the reference: $n_{db} = 2$, and the dead-beat time is 

$t_{db} = n_{db} ~ \Delta T = 1 \quad \text{second.}$"

# ╔═╡ Cell order:
# ╟─1b3761c0-5b34-11eb-3a51-915a1fb66273
# ╟─78a8935a-5b2b-11eb-0bc0-157b068a2335
# ╠═e4956e92-5b32-11eb-0962-b36f35ccd3bb
# ╠═08b76d6a-5b34-11eb-3489-55eedc98d28a
# ╠═f1906638-5b32-11eb-2006-870fdd2b0f43
# ╠═337dfa6a-5b33-11eb-2586-7d378fb12025
# ╟─1b1e2eda-5b34-11eb-2f93-353e2d242f61
# ╠═bbdf95c4-5b3a-11eb-0513-71595f7699b8
# ╠═bbc45868-5b3a-11eb-3653-c7b015841dd8
# ╠═bba3ba54-5b3a-11eb-212f-25d5ec2e3512
# ╠═bb807cba-5b3a-11eb-2d5f-b160c1de903d
# ╠═bdfe6152-5b3b-11eb-2ab8-9507d5ad5a60
# ╠═bb61b212-5b3a-11eb-13fd-4727abfac6c4
# ╟─188ebe40-5b3c-11eb-11c0-03ab24de191c
