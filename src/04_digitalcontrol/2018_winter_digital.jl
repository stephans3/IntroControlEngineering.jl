### A Pluto.jl notebook ###
# v0.12.12

using Markdown
using InteractiveUtils

# ╔═╡ ccaaea7e-5b04-11eb-39a1-efa7a38b778a
using Plots

# ╔═╡ 19821148-5b04-11eb-3fa4-dd2ff13141fe
md"# Digital Control

The controller transfer function

$G_{c}(s) = \frac{U(s)}{X(s)} = \frac{s + 1}{s + 2} = \frac{s}{s + 2} + \frac{1}{s + 2}$

is transfered to the time domain with the [inverse Laplace transform](https://en.wikipedia.org/wiki/Laplace_transform#Table_of_selected_Laplace_transforms) as

$\mathcal{L}^{-1} \left\{\frac{G(s)}{s} \right\} = \mathcal{L}^{-1} \left\{ \frac{1}{s + 2} + \frac{1}{2} \frac{2}{s~(s + 2)} \right\} = \frac{1}{2} \left[1 + \exp(-2~t) \right] = h(t)$

where $h(t)$ is known as the step response.

## Simulation of the step response
"

# ╔═╡ 8e97827e-5b04-11eb-0fa9-518c196d4c81
N = 200    # Number of time steps

# ╔═╡ 8e7e4d34-5b04-11eb-0d62-0990e2e4da07
ΔT = 0.01  # Sampling time

# ╔═╡ 8df5a332-5b04-11eb-3417-556436b55f0d
tspan = 0:ΔT:ΔT*(N-1) # Time range

# ╔═╡ 8da2ee06-5b04-11eb-183b-f7337f8e1fd4
h(t) = y = 0.5*(1 .+ exp.(-2*t)) # Impulse response

# ╔═╡ d12a4c70-5b04-11eb-37b1-d7b2944c57b1
y = h.(tspan) # Output

# ╔═╡ e2ba017e-5b04-11eb-2897-93f64e37c334
plot(tspan, y)

# ╔═╡ f14eacb2-5b04-11eb-0921-e3650b9c8a02
md"## Zero-order hold and Z-Transform

The zero-order hold method for the conversion of an analog $G(s)$ to a digital systems $G^{*}(z)$ is defined as

$G^{*}(z) = \frac{z-1}{z} ~ \mathcal{Z} \left\{ \left. \mathcal{L}^{-1} \left\{\frac{G(s)}{s} \right\}\right\rvert_{t=nT} \right\}.$

The [Z-transform](https://en.wikipedia.org/wiki/Z-transform) is applyied on the time-continuous step response to yield

$\mathcal{Z} \left\{ \frac{1}{2} \left[1 + \exp(-2~t) \right]_{t = n~T}\right\} = \frac{1}{2} \left( \frac{z}{z - 1} + \frac{z}{z - \exp\left(-2T\right)} \right)$

this result is multiplied with $\frac{z-1}{z}$ and further reformulated as 

$\frac{U(z)}{X(z)} = \frac{z - 0.5~\left(1 + \exp(-2T)\right)}{z - \exp(-2T)}$

and the resulting discrete transfer function is noted as

$G_{c}^{*}(z) = \frac{U(z)}{X(z)} = \frac{1 - 0.5~\left(1 + e^{-2T}\right) z^{-1} }{1 - e^{-2T}  z^{-1}} \text{.}$

### Control algorithm

The discrete transfer function is rewritten as

$U(z) \left[ 1 - e^{-2T}  z^{-1} \right] = X(z) \left[ 1 - 0.5~\left(1 + e^{-2T}\right) z^{-1}\right]$

and the inverse Z-Transform is applied to gain the algorithm

$u_{n} = e^{-2T} ~ u_{n-1} + x_{n} -  0.5~\left(1 + e^{-2T}\right) ~ x_{n-1} \text{.}$

### Simulation of the Control algorithm
"

# ╔═╡ 543b489e-5b05-11eb-08cb-b39bf900c13a
N2 = 50 # Number of time steps

# ╔═╡ 542097b0-5b05-11eb-1e37-9b2ac0f12d93
ΔT2 = 0.1 # Sampling time

# ╔═╡ 5402014c-5b05-11eb-1615-93b11b081016
x = ones(N2) # Input of controller

# ╔═╡ 53e493be-5b05-11eb-17e6-29a7fc1e29bf
u = zeros(N2) # Output of controller

# ╔═╡ 53c76852-5b05-11eb-0593-4f7f62e8c9b1
for i = 1 : N2-1
    u[i+1] = exp(-2*ΔT2)*u[i] + x[i+1] - 0.5*(1 + exp(-2*ΔT2))*x[i]
end

# ╔═╡ f12a2612-5b04-11eb-0d5a-99b38eb7fc0e
tspan2 = 0 : ΔT2 : ΔT2 * (N2-1) # Time range for plot

# ╔═╡ f107e1ec-5b04-11eb-294f-c74672bb32df
plot(tspan2, [x, u], label=["input" "output"])

# ╔═╡ 89584bf4-5bcc-11eb-3081-6f0bb2385c0d
T5 = 0.1

# ╔═╡ 91c08ca2-5bcc-11eb-265e-b7a71146d2d8
s5 = -0.1 + 100im

# ╔═╡ 9c6fc10e-5bcc-11eb-051b-970da172e2db
z5 = exp(s5*T5)

# ╔═╡ c895c454-5bcc-11eb-169c-83b28efab43b
abs(z5)

# ╔═╡ f067da4e-5b04-11eb-3e00-275c60c5c86b
md"## Tustin's method

The relationship between the continuous Laplace and discrete Z domain is discribed by 

$z = e^{s~T}$

which can also be noted as

$s T = \ln( z ) \quad \text{or equivalent} \quad s = \frac{1}{T} \ln(z).$

The term $\ln(z)$ can be calculated with [power series](https://en.wikipedia.org/wiki/Logarithm#Power_series) as

$\ln(z) = 2 \sum\limits_{n = 0}^{\infty} \frac{1}{2~n + 1} \left(\frac{z-1}{z+1}\right)^{2~n + 1}.$

The [Tustin transform](https://en.wikipedia.org/wiki/Bilinear_transform) uses the approximation of the power series approach as

$s \approx \frac{2}{T} \frac{z - 1}{z + 1}$

to transfer the transfer function from the continuous frequency domain directly to the discrete frequency domain as 

$\tilde{G}_{c}(z) = \left. G_{c}(s) \right\rvert_{s =  \frac{2}{T} \frac{z - 1}{z + 1}} = \frac{\frac{2}{T} \frac{z - 1}{z + 1} + 1}{\frac{2}{T} \frac{z - 1}{z + 1} + 2} = \cdots = \frac{ \frac{21}{22} - \frac{19}{22} z^{-1} }{ 1 - \frac{9}{11} z^{-1} } = \frac{\tilde{U}(z)}{\tilde{X}(z)}$

with $ T = 0.1 $ seconds.

### Control algorithm

The control algorithm is derived from the discrete transfer function using the inverse Z-Transform to yield

$u_{n} = \frac{9}{11} ~ u_{n-1} + \frac{21}{22} ~ x_{n} - \frac{19}{22} ~ x_{n-1} \text{.}$

### Simulation
"

# ╔═╡ 361f0314-5b0a-11eb-05bf-d1607f0ee111
N3 = 50 # Number of time steps

# ╔═╡ 360405d2-5b0a-11eb-1d78-31c2cc2b20c5
ΔT3 = 0.1 # Sampling time

# ╔═╡ 35e1ea1a-5b0a-11eb-39bc-53c1ae54870e
x3 = ones(N3); # Input of controller

# ╔═╡ 35c70bbe-5b0a-11eb-20d1-9f2b738823d9
u3 = zeros(N3); # Output of controller

# ╔═╡ 35aaf3ac-5b0a-11eb-3ff2-a7b1ffa21387
for i = 1 : N3-1
	u3[i+1] = (9/11)*u3[i] + (21/22)*x3[i+1] - (19/22)*x3[i]
end

# ╔═╡ 359041f6-5b0a-11eb-05f2-bd3076fc58fc
tspan3 = 0 : ΔT3 : ΔT3 * (N3-1) # Time range for plot

# ╔═╡ 35767640-5b0a-11eb-2c1d-ab803f859a2b
plot(tspan3, [x3, u3], label=["input" "output"])

# ╔═╡ 3527112e-5b0a-11eb-2f86-e9baa50c3873
u - u3

# ╔═╡ 34b0aa0a-5b0a-11eb-10d2-fbf75b2f064c


# ╔═╡ Cell order:
# ╟─19821148-5b04-11eb-3fa4-dd2ff13141fe
# ╠═8e97827e-5b04-11eb-0fa9-518c196d4c81
# ╠═8e7e4d34-5b04-11eb-0d62-0990e2e4da07
# ╠═8df5a332-5b04-11eb-3417-556436b55f0d
# ╠═8da2ee06-5b04-11eb-183b-f7337f8e1fd4
# ╠═d12a4c70-5b04-11eb-37b1-d7b2944c57b1
# ╠═ccaaea7e-5b04-11eb-39a1-efa7a38b778a
# ╠═e2ba017e-5b04-11eb-2897-93f64e37c334
# ╟─f14eacb2-5b04-11eb-0921-e3650b9c8a02
# ╠═543b489e-5b05-11eb-08cb-b39bf900c13a
# ╠═542097b0-5b05-11eb-1e37-9b2ac0f12d93
# ╠═5402014c-5b05-11eb-1615-93b11b081016
# ╠═53e493be-5b05-11eb-17e6-29a7fc1e29bf
# ╠═53c76852-5b05-11eb-0593-4f7f62e8c9b1
# ╠═f12a2612-5b04-11eb-0d5a-99b38eb7fc0e
# ╠═f107e1ec-5b04-11eb-294f-c74672bb32df
# ╠═89584bf4-5bcc-11eb-3081-6f0bb2385c0d
# ╠═91c08ca2-5bcc-11eb-265e-b7a71146d2d8
# ╠═9c6fc10e-5bcc-11eb-051b-970da172e2db
# ╠═c895c454-5bcc-11eb-169c-83b28efab43b
# ╟─f067da4e-5b04-11eb-3e00-275c60c5c86b
# ╠═361f0314-5b0a-11eb-05bf-d1607f0ee111
# ╠═360405d2-5b0a-11eb-1d78-31c2cc2b20c5
# ╠═35e1ea1a-5b0a-11eb-39bc-53c1ae54870e
# ╠═35c70bbe-5b0a-11eb-20d1-9f2b738823d9
# ╠═35aaf3ac-5b0a-11eb-3ff2-a7b1ffa21387
# ╠═359041f6-5b0a-11eb-05f2-bd3076fc58fc
# ╠═35767640-5b0a-11eb-2c1d-ab803f859a2b
# ╠═3527112e-5b0a-11eb-2f86-e9baa50c3873
# ╠═34b0aa0a-5b0a-11eb-10d2-fbf75b2f064c
