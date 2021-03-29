### A Pluto.jl notebook ###
# v0.12.12

using Markdown
using InteractiveUtils

# ╔═╡ bc2d7c20-6539-11eb-03d9-6967ca513380
using ControlSystems, Plots

# ╔═╡ dc795b48-6630-11eb-2d8e-71c92ca5cc94
md"# Root-locus: Analysis

Consider a dynamical system with transfer function

$G_{OL}(s) = K \frac{s^3 + 9 s^2 + 31 s + 39}{s^3 + 3 s^2 + 7 s + 5}.$

## Root-locus diagram
"

# ╔═╡ e502b996-6539-11eb-2d57-9b7a8c4eb7c9
Kp = 1

# ╔═╡ d744ff62-6539-11eb-0965-41ad6960220e
Gol = Kp * tf([1, 9, 31, 39],[1, 3, 7, 5])

# ╔═╡ d7444676-6539-11eb-3e81-af80999770fb
rlocus(Gol,yticks = -3:0.5:4, ylims=(-3.0, 4.), xlims=(-3.5, 0.5))

# ╔═╡ 64a6974c-6631-11eb-0fd8-9982cbba6a4d
md"The poles of the open-loop system are $\{-1., -1 \pm 2j\}$ and the zeros are $\{-3., -3 \pm 2j\}$. Using these poles and zeros, one finds that

1. for $K=0$ the poles of the closed-loop system are $\{-1., -1 \pm 2j\}$ and
2. for $K\rightarrow \infty$ the poles of the closed-loop system are $\{-3., -3 \pm 2j\}$.

## Finding gain values

If one seeks for the gain of some pole, e.g. $s_{1} = -2 + 2.65j$, then the closed-loop transfer function has to be determined as  

$G_{CL}(s) = \frac{G_{OL}(s)}{1 + G_{OL}(s)} =  \frac{K (s^3 + 9s^2 + 31s + 39)}{(s^3 + 3s^2 + 7s + 5) + K (s^3 + 9s^2 + 31s + 39)}.$

The characteristic polynomial is set to zero 

$p(s) = (s^3 + 3s^2 + 7s + 5) + K (s^3 + 9s^2 + 31s + 39) = 0$

and the resulting gain is found with

$K = \frac{-(s^3 + 3s^2 + 7s + 5)}{s^3 + 9s^2 + 31s + 39}.$"

# ╔═╡ 15c4ce78-653c-11eb-15c7-e77e9af3ea9f
K(s) = -(s^3 + 3*s^2 + 7*s + 5)/(s^3 + 9*s^2 + 31*s + 39)

# ╔═╡ 6eb8e884-653c-11eb-3318-1fcd1e924cd1
s₁ = -2 + 2.65im # Pole on the root-locus

# ╔═╡ 081d5400-6633-11eb-0e0b-8334058e42d8
md"**Note:** Only the real part of the gain is of interest."

# ╔═╡ 6842856e-653c-11eb-0f4d-c505d0e9d7e0
K₁ = real(K(s₁))

# ╔═╡ d24f4532-653e-11eb-00d1-95f610b7b878
md"## Damping / Attenuation

The poles of an arbitrary second-order system (e.g. oscillator)

$G(s) = \frac{K_{p}}{1 + 2 \zeta T s + T^2 s^2}$

with attenuation $\zeta \geq 0$ and time constant $T > 0$ is found with

$s_{1,2} = - \frac{\zeta}{T} \pm \frac{j}{T} \sqrt{1 - \zeta^2} = \rho \pm j \omega.$

The damping of an arbitrary system (not necessarily of second order) can be calculated as

$\zeta = \frac{-\rho}{\sqrt{\rho^2 + \omega^2}} = \sin(\alpha)$

where $\alpha$ represents the angle between the imaginary axis and a line connecting the origin and desired pole.  

In this example, the attenuation of the pole $s_{1} = 2.0 + 2.65j$ shall be found where
- the real part of $s_{1}: \rho = 2.0$ and
- the imaginary part of $ s_{1}: \omega = 2.65$.

"

# ╔═╡ 18176586-653d-11eb-087a-e77b641ce663
ζ = 2/sqrt(2^2 + 2.65^2) # Attenuation

# ╔═╡ 4689266c-6541-11eb-2fa9-7310869b57ce
md"If the attenuation is given and the controller gain $K$ is desired, then the angle $\alpha$ is calculated and the corresponding pole on the root locus is found graphically.

Assume the attenuation 

$\zeta_{2} = 0.7 = \sin(\alpha)$

which corresponds to the angle

$\alpha[°] = \arcsin(0.7) ~ \frac{180°}{\pi} \approx 45°.$
"

# ╔═╡ eabe261e-653d-11eb-3f8e-2577be1f5e22
α₂ = asin(0.7)*180/pi

# ╔═╡ 77176110-6639-11eb-10e0-bf5ad0206816
md"Drawing a line from the origin to the root locus with an angle $\alpha \approx 45°$  results in a pole at $s_{2} \approx -2.5 \pm 2.5j$."

# ╔═╡ 99a59a1a-6541-11eb-193e-0ba6724c38d8
s₂ =-2.5 + 2.5im

# ╔═╡ bd075c3e-6639-11eb-0348-e9f51523b899
md"Finally, the corresponding pole is found as above with 

$K = \frac{-(s^3 + 3s^2 + 7s + 5)}{s^3 + 9s^2 + 31s + 39}.$"

# ╔═╡ dc7cdf46-6542-11eb-2857-97c3f7177481
K₂ = real(K(s₂))

# ╔═╡ e580391c-6542-11eb-0e14-4b468d69cdf3


# ╔═╡ Cell order:
# ╟─dc795b48-6630-11eb-2d8e-71c92ca5cc94
# ╠═bc2d7c20-6539-11eb-03d9-6967ca513380
# ╠═e502b996-6539-11eb-2d57-9b7a8c4eb7c9
# ╠═d744ff62-6539-11eb-0965-41ad6960220e
# ╠═d7444676-6539-11eb-3e81-af80999770fb
# ╟─64a6974c-6631-11eb-0fd8-9982cbba6a4d
# ╠═15c4ce78-653c-11eb-15c7-e77e9af3ea9f
# ╠═6eb8e884-653c-11eb-3318-1fcd1e924cd1
# ╟─081d5400-6633-11eb-0e0b-8334058e42d8
# ╠═6842856e-653c-11eb-0f4d-c505d0e9d7e0
# ╟─d24f4532-653e-11eb-00d1-95f610b7b878
# ╠═18176586-653d-11eb-087a-e77b641ce663
# ╟─4689266c-6541-11eb-2fa9-7310869b57ce
# ╠═eabe261e-653d-11eb-3f8e-2577be1f5e22
# ╟─77176110-6639-11eb-10e0-bf5ad0206816
# ╠═99a59a1a-6541-11eb-193e-0ba6724c38d8
# ╟─bd075c3e-6639-11eb-0348-e9f51523b899
# ╠═dc7cdf46-6542-11eb-2857-97c3f7177481
# ╠═e580391c-6542-11eb-0e14-4b468d69cdf3
