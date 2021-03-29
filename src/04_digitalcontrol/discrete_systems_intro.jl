### A Pluto.jl notebook ###
# v0.12.12

using Markdown
using InteractiveUtils

# ╔═╡ 0baf0bb8-599a-11eb-2b22-4977da7a6ab0
using Plots

# ╔═╡ 03163c02-59b4-11eb-0b97-1deb1454a296
using ControlSystems

# ╔═╡ c36c5df0-5995-11eb-1d7c-75900aa19960
md"# Discrete Systems

This Pluto notebook gives an introduction to 
- discretization with Euler method: [forward](https://en.wikipedia.org/wiki/Euler_method) and [backward](https://en.wikipedia.org/wiki/Backward_Euler_method)
- [Z-transform](https://en.wikipedia.org/wiki/Z-transform) 
- [Digital control systems](https://en.wikipedia.org/wiki/Digital_control)

## Dynamical systems

Consider a dynamical system 

$\dot{x}(t) = a ~ x(t) \quad \text {with} \quad x(0) = x_{0} > 0$

and some $a \in \mathbb{R}.$ Applying the **forward Euler scheme**

$\frac{d x(t)}{dt} \approx \frac{x(t_{k+1}) - x(t_{k})}{ \Delta T} = f(t, x(t_{k}))$

with sampling time $\Delta T > 0$ leads to 

$\frac{x(t_{k+1}) - x(t_{k})}{ \Delta T} = a ~ x(t_{k})$

or equivalent

$x(t_{k+1}) = x(t_{k}) + a ~ \Delta T ~ x(t_{k}) = (1 + a~\Delta T) ~ x(t_{k}).$

This time-discrete formula can be noted as an iterative algorithm

$x(t_{0}) = x_{0} \quad \text{(initial value)}$

$x(t_{1}) = (1 + a~\Delta T) ~ x(t_{0})$

$x(t_{2}) = (1 + a~\Delta T) ~ x(t_{1}) = (1 + a~\Delta T)^{2} ~ x(t_{0})$

$x(t_{3}) = (1 + a~\Delta T) ~ x(t_{2}) = (1 + a~\Delta T)^{3} ~ x(t_{0})$

$\vdots$

$x(t_{n}) = (1 + a~\Delta T) ~ x(t_{n-1}) = (1 + a~\Delta T)^{n} ~ x(t_{0})$

"

# ╔═╡ 51ae8768-5999-11eb-072e-2f3d10b4e33a
md"### Example: Comparison of formula and algorithm"

# ╔═╡ 672100c4-5999-11eb-184f-495f49fbc0fc
a = -2.0

# ╔═╡ 6e3e07ee-5999-11eb-01ca-4b480d92ba7a
ΔT = 0.2 # Sampling time

# ╔═╡ 763379b6-5999-11eb-01b3-f3cb44a5cbf1
x₀ = 4.0 # initial value

# ╔═╡ 8677176a-5999-11eb-1685-7dfc07d07388
N = 20 # Maximum number of steps

# ╔═╡ ba08d69c-5999-11eb-3b89-b1006b8952a6
x = zeros(N+1);

# ╔═╡ b9c25788-5999-11eb-04ee-176a63d2afe5
x[1] = x₀

# ╔═╡ b9a46656-5999-11eb-2884-a567947f53b1
begin
	for i = 1 : N-1
		x[i+1] = x[i] + ΔT * a*x[i]
	end
end

# ╔═╡ 206505e6-599a-11eb-2cf0-070687e5a1ae
x̃ = zeros(N+1);

# ╔═╡ 1dd1a510-599a-11eb-2973-49ff4ed6510d
x̃[1] = x₀

# ╔═╡ 4700d62e-599a-11eb-1bea-5f81797ea86e
for i = 1 : N-1
	x̃[i+1] = (1 + a*ΔT)^(i)*x₀
end

# ╔═╡ 0c7503c2-599a-11eb-0d77-3522746efbb6
scatter(x, label="x")

# ╔═╡ 6f901cee-599a-11eb-3ef8-c169f8557a51
scatter(x̃, label="x̃")

# ╔═╡ ead0db4e-59a6-11eb-213a-c1ded424ee47
md"## Numerical Stability

The time-discrete equation 

$x(t_{k+1}) = x(t_{k}) + a ~ \Delta T ~ x(t_{k}) = (1 + a~\Delta T) ~ x(t_{k})$

is called [numerical stable](https://en.wikipedia.org/wiki/Numerical_stability) if 

$(1 + a~\Delta T) > -1.$

This condition guarantees that the term $(1 + a~\Delta T)^{n}$ in the algorithm  

$x(t_{n}) = (1 + a~\Delta T)^{n} ~ x(t_{0})$

is always greater than zero and there **does not occur** any trembling."

# ╔═╡ ebd0a2e0-59a6-11eb-0d78-15c1f496ca00
q = 1 + a*ΔT

# ╔═╡ 0a34b62e-599b-11eb-20c9-35bac2c7bbd1
md"## Backward Euler method

The backward Euler scheme uses the recent (or next) time step on both sides as 

$\frac{d x(t)}{dt} \approx \frac{x(t_{k+1}) - x(t_{k})}{ \Delta T} = f(t, x(t_{k+1})).$

For the given example one holds 

$\frac{x(t_{k+1}) - x(t_{k})}{ \Delta T} = a ~ x(t_{k+1})$

or equivalent

$x(t_{k+1}) - a ~ \Delta T ~ x(t_{k+1}) = (1 - a ~ \Delta T) ~ x(t_{k+1}) = x(t_{k}).$

Finally the time-discrete equation is noted as

$x(t_{k+1}) = (1 - a ~ \Delta T)^{-1} ~ x(t_{k})$

and the algorithm is derived analog to the forward Euler method as

$x(t_{n}) = (1 - a ~ \Delta T)^{-n} ~ x(t_{0}).$

If $a < 0$ - this means the system is bounded-input-bounded-output (BIBO) stable - then the algorithm is stable for all $\Delta T > 0$.

### Example
"

# ╔═╡ 4d96ed44-599c-11eb-0ca6-53f505e22809
a2 = -2.0

# ╔═╡ 262db5e0-59a0-11eb-0891-45b271f95b17
ΔT2 = 1.5

# ╔═╡ 2fa8406a-59a0-11eb-2537-e577b3523bd5
x2 = zeros(N + 1);

# ╔═╡ 6d08835c-59a0-11eb-385d-0b1934cb0179
x2[1] = x₀

# ╔═╡ 465fc4f6-59a0-11eb-007f-07114f74ace4
for i = 1 : N-1
	x2[i+1] = (1 - a2*ΔT2)^(-i) * x₀
end

# ╔═╡ 85901f6e-59a0-11eb-18dd-732c8dd36219
scatter(x2, legend=false)

# ╔═╡ c4241fe8-59a0-11eb-241e-0fe5f7a95da0
md"**Remark:** The effort computation of the inverse can be costly in case of large systems."

# ╔═╡ 099ed2b6-59a1-11eb-3fc9-f9592ceb4c10
md"## Z-Transform

Consider the time-discrete equation 

$x(t_{k+1}) = (1 + a~\Delta T) ~ x(t_{k})$

with $x_{0} > 0$ and $a \in \mathbb{R}$. The [Z-transform](https://en.wikipedia.org/wiki/Z-transform) is applied analog to the Laplace-transform as

$z ~ X(z) - z ~ x_{0} = (1 + a~\Delta T) ~ X(z).$

This equation is reformulated to 

$\left[z - (1 + a~\Delta T)\right] ~ X(z) = z ~ x_{0}$

or equivalent

$X(z) = \frac{z}{z - (1 + a~\Delta T)} ~ x_{0} = \frac{1}{1 - (1 + a~\Delta T) z^{-1}} ~ x_{0}.$

The **pole** of transfer function $G(z) = \frac{X(z)}{x_{0}}$ is found with 

$z - (1 + a~\Delta T) = 0 \quad \text{equivalent to} \quad z = 1 + a~\Delta T$

and the pole is called **stable** if $\lvert z \rvert < 1$. As one notes, the discrete system has two parts to be called stable:
- the discrete system $x(t_{k+1}) = f(t_{k}, x)$ shall not grow to infinity which is comparable with the **BIBO stability** and
- the discrete system has to be **numerical stable**

"

# ╔═╡ 0982c904-59a1-11eb-146d-05c7d46cbf5e
md"## Digital Control Systems

Consider the continuous control system

$\dot{x}(t) = a~x(t) + b~u(t)$
$y(t) = c~x(t)$

with state $x$, input $u$ and output $y$. The differential equation describing the variation of the state is solved by integration as 

$x(t) = e^{a~(t-t_{0})} ~ x_{0} + \int\limits_{t_{0}}^{t} e^{a~(t-s)} ~ b ~ u(s) ~ ds.$

Now, the continuous time is discretized as $t = n~\Delta T$ with sampling time $\Delta T >0$ and discrete step $n = 0, 1 , 2, \cdots N-1$. Consequently, one yields the discrete equation

$x(n \Delta T) = e^{a~n~\Delta T} ~ x_{0} + \int\limits_{0}^{n~\Delta T} e^{a~(n~\Delta T-s)} ~ b ~ u(s) ~ ds$

which is **approximated** as

$x(n \Delta T) = e^{a~n~\Delta T} ~ x_{0} + \sum\limits_{i=0}^{n-1} e^{a~(n-i)~\Delta T} ~ b ~ u(i \Delta T) ~ \Delta T.$

**Remark:** Note the equality $e^{a~n~\Delta T} = (1 + a \Delta T)^{n}$. 

### Example

"

# ╔═╡ 096441d2-59a1-11eb-1339-17efaf3e04b6
md"a = $(a)"

# ╔═╡ 9a43493a-59b0-11eb-246b-0bfc5845c05e
md"ΔT = $(ΔT)"

# ╔═╡ 738d5bb0-59b1-11eb-2d83-a561b6593506
md"Stability (numerical + BIBO): q=$(q) 

$q \in (-1,1)$"

# ╔═╡ 0945cac2-59a1-11eb-1dea-8793877ace76
md"Initial condition $x_{0}$ = $(x₀)"

# ╔═╡ 0928088e-59a1-11eb-31e8-79991e949db3
md"Input $u(t) = 1$ for $t > 0$."

# ╔═╡ 08ddedf8-59a1-11eb-2586-575a87b0e83f
u = ones(N); # input

# ╔═╡ 748ab37c-59b2-11eb-25cb-bbb341f2bb76
b = 1.0 # input vector

# ╔═╡ 08c048de-59a1-11eb-250f-abd3b136661c
x3 = zeros(N+1);

# ╔═╡ 5d38ab58-59b1-11eb-2914-25c083aa1be5
x3[1] = 0. # Initial value

# ╔═╡ 69e7c41a-59b1-11eb-264a-ab07131fc8ce
for n=1 : N
	in_sum = 0
	for i = 1 : n
		in_sum = in_sum + exp(a*(n-i)*ΔT) * b * u[i]*ΔT
	end
	
	x3[n+1] = exp(a*n*ΔT)*x3[1] + in_sum
end

# ╔═╡ 1f70d9f8-5a3d-11eb-00e1-dbd4bec0c12c


# ╔═╡ 9ffa7f1a-59b2-11eb-1f8f-7903602fa03c
scatter(x3, legend=false)

# ╔═╡ 6b2d6256-59b5-11eb-0bd2-0d10936747d6
md"### Compare with higher-order solver"

# ╔═╡ 052800d6-5a32-11eb-1395-0f6569e077a2
sys = ss(a, b, 1, 0)

# ╔═╡ 1a674c86-5a32-11eb-218b-87b9b8b1cf16
stepplot(sys, 4., legend=false)

# ╔═╡ 3663799a-5a3d-11eb-3e46-9be3224c6638
md"#### Discussion

Obviously, there is some difference between the approach using a sum and the ODE solver method. One reason of this difference is the **approximation** of the integral using the sum. The smaller the sampling time $\Delta T$ is chosen the more precise is the solution.

## Digital Control Algorithm

Now, a algorithm is introduced to solve the digital control exactly. Instead of solving the dynamical system for the complete time span as 

$x(t) = e^{a~(t-t_{0})} ~ x_{0} + \int\limits_{t_{0}}^{t} e^{a~(t-s)} ~ b ~ u(s) ~ ds$

the system only solved from one time step to the next as 

$x(t_{k+1}) = e^{a~ \Delta T} ~ x(t_{k}) + \int\limits_{t_{k}}^{t_{k+1}} e^{a~(t_{k+1}-s)} ~ b ~ u(s) ~ ds$

with $\Delta T = t_{k+1} - t_{k}$. Firstly, a new variable $\tau := s - t_{k}$ and the solution is noted as

$x(t_{k+1}) = e^{a~ \Delta T} ~ x(t_{k}) + \int\limits_{0}^{\Delta T} e^{a~(\Delta T - \tau)} ~ b ~ d\tau ~ u(t_{k}).$

The input $u(t_{k})$ can be noted outside the integral because it is constant during the sampling time ([Zero-order hold](https://en.wikipedia.org/wiki/Zero-order_hold)). The integral is solved as 

$\int\limits_{0}^{\Delta T} e^{a~(\Delta T - \tau)} ~ b ~ d\tau = \frac{b}{a}\left(e^{a\Delta T} - 1\right)$

and the new coefficients of the digital state-space representation are
1. state coefficient: $a_{d} = e^{a~ \Delta T}$
2. input coefficient: $b_{d} = \frac{b}{a}\left(e^{a\Delta T} - 1\right)$
3. output coefficient: $c_{d} = c$

and the digital control algorithm is derived as

$x_{k+1} = a_{d} ~ x_{k} + b_{d} ~ u_{k} \text{,}$
$y_{k} = c_{d} ~ x_{k} \text{.}$

### Implementation
"

# ╔═╡ 9768a30c-5a3b-11eb-0fa3-9d567fcd19cc
x4 = zeros(N+1);

# ╔═╡ 8bb6a22e-5a41-11eb-182b-e9f53551d675
a4 = -2.0;

# ╔═╡ 976e8f32-5a41-11eb-3ff7-c7a498eaaa1f
b4 = 1.0;

# ╔═╡ 6a88fce6-5a41-11eb-033a-7d1f857316b2
ad4 = exp(a*ΔT) # State coefficient

# ╔═╡ 974d4cb0-5a3b-11eb-14a4-cb8d91a8ada9
bd4 = b/a * (exp(a*ΔT) - 1) # Input coefficient

# ╔═╡ 972e8af0-5a3b-11eb-0cc2-d721dfaafff9
for n = 1 : N
	x4[n+1] = ad4*x4[n] + bd4*u[n]
end

# ╔═╡ 9714a4a0-5a3b-11eb-37bc-e3b037c71637
scatter(x4, legend=false)

# ╔═╡ Cell order:
# ╟─c36c5df0-5995-11eb-1d7c-75900aa19960
# ╟─51ae8768-5999-11eb-072e-2f3d10b4e33a
# ╠═672100c4-5999-11eb-184f-495f49fbc0fc
# ╠═6e3e07ee-5999-11eb-01ca-4b480d92ba7a
# ╠═763379b6-5999-11eb-01b3-f3cb44a5cbf1
# ╠═8677176a-5999-11eb-1685-7dfc07d07388
# ╠═ba08d69c-5999-11eb-3b89-b1006b8952a6
# ╠═b9c25788-5999-11eb-04ee-176a63d2afe5
# ╠═b9a46656-5999-11eb-2884-a567947f53b1
# ╠═206505e6-599a-11eb-2cf0-070687e5a1ae
# ╠═1dd1a510-599a-11eb-2973-49ff4ed6510d
# ╠═4700d62e-599a-11eb-1bea-5f81797ea86e
# ╠═0baf0bb8-599a-11eb-2b22-4977da7a6ab0
# ╠═0c7503c2-599a-11eb-0d77-3522746efbb6
# ╠═6f901cee-599a-11eb-3ef8-c169f8557a51
# ╟─ead0db4e-59a6-11eb-213a-c1ded424ee47
# ╠═ebd0a2e0-59a6-11eb-0d78-15c1f496ca00
# ╟─0a34b62e-599b-11eb-20c9-35bac2c7bbd1
# ╠═4d96ed44-599c-11eb-0ca6-53f505e22809
# ╠═262db5e0-59a0-11eb-0891-45b271f95b17
# ╠═2fa8406a-59a0-11eb-2537-e577b3523bd5
# ╠═6d08835c-59a0-11eb-385d-0b1934cb0179
# ╠═465fc4f6-59a0-11eb-007f-07114f74ace4
# ╠═85901f6e-59a0-11eb-18dd-732c8dd36219
# ╟─c4241fe8-59a0-11eb-241e-0fe5f7a95da0
# ╟─099ed2b6-59a1-11eb-3fc9-f9592ceb4c10
# ╟─0982c904-59a1-11eb-146d-05c7d46cbf5e
# ╠═096441d2-59a1-11eb-1339-17efaf3e04b6
# ╠═9a43493a-59b0-11eb-246b-0bfc5845c05e
# ╠═738d5bb0-59b1-11eb-2d83-a561b6593506
# ╟─0945cac2-59a1-11eb-1dea-8793877ace76
# ╟─0928088e-59a1-11eb-31e8-79991e949db3
# ╠═08ddedf8-59a1-11eb-2586-575a87b0e83f
# ╠═748ab37c-59b2-11eb-25cb-bbb341f2bb76
# ╠═08c048de-59a1-11eb-250f-abd3b136661c
# ╠═5d38ab58-59b1-11eb-2914-25c083aa1be5
# ╠═69e7c41a-59b1-11eb-264a-ab07131fc8ce
# ╠═1f70d9f8-5a3d-11eb-00e1-dbd4bec0c12c
# ╠═9ffa7f1a-59b2-11eb-1f8f-7903602fa03c
# ╟─6b2d6256-59b5-11eb-0bd2-0d10936747d6
# ╠═03163c02-59b4-11eb-0b97-1deb1454a296
# ╠═052800d6-5a32-11eb-1395-0f6569e077a2
# ╠═1a674c86-5a32-11eb-218b-87b9b8b1cf16
# ╟─3663799a-5a3d-11eb-3e46-9be3224c6638
# ╠═9768a30c-5a3b-11eb-0fa3-9d567fcd19cc
# ╠═8bb6a22e-5a41-11eb-182b-e9f53551d675
# ╠═976e8f32-5a41-11eb-3ff7-c7a498eaaa1f
# ╠═6a88fce6-5a41-11eb-033a-7d1f857316b2
# ╠═974d4cb0-5a3b-11eb-14a4-cb8d91a8ada9
# ╠═972e8af0-5a3b-11eb-0cc2-d721dfaafff9
# ╠═9714a4a0-5a3b-11eb-37bc-e3b037c71637
