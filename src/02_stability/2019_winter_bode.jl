### A Pluto.jl notebook ###
# v0.12.12

using Markdown
using InteractiveUtils

# ╔═╡ 9502a1d8-3fa6-11eb-187d-2dff3e492226
using Plots

# ╔═╡ 261fff78-3fac-11eb-1e00-bb85b363d664
using ControlSystems

# ╔═╡ 95230452-3fa6-11eb-114b-e300526fdffc
md"# Example 1: Bode plot

Assume the transfer function 

$G_{OL}(s) = K \frac{(s + 10)^2}{(s + 1)^3}$

with $K > 0$. The open loop system $G_{OL}$ is BIBO stable because it has only poles smaller than zero: $s_{i} = -1$ for $i = 1, 2, 3$.

The transfer function $G_{OL}$ consists of three types of transfer function:
- Static gain $K > 0$,
- PT$_{1}$ system $G_{1}(s) = \frac{1}{s + 1}$ and
- DT$_{1}$ system $G_{2}(s) = (s + 10)$.

The latter two systems in the frequency domain yield the form

$G_{1}(j \omega) = \frac{1}{1 + j \omega} = \frac{1}{1 + \omega^2} - j \frac{\omega}{1 + \omega^2}$

and

$G_{2}(j \omega) = 10 + j \omega \text{.}$"

# ╔═╡ 2c3f0c66-3efb-11eb-1c40-9555f67871b6
md"

## Amplitude

The amplitude of a transfer function is found by calculating

$\lvert G(j\omega) \rvert = \sqrt{ \mathcal{Re}\{G(j\omega)\}^2 + \mathcal{Im}\{G(j\omega)\}^2}$

and thus one notes 

$\lvert G_{1}(j\omega) \rvert = \frac{\sqrt{1 + \omega^2}}{1 + \omega^2} = \frac{1}{\sqrt{1 + \omega^2}} \quad \text{and}$

$\lvert G_{2}(j\omega) \rvert = \sqrt{100 + \omega^2}$

The amplitude of the open loop system at $\omega=0$ is simply calculated by 

$\lvert G_{OL}(0) \rvert = K ~ \lvert G_{1}(0) \rvert^3 ~ \lvert G_{2}(0) \rvert^2 = 100~K.$

The amplitude in Decibel [dB] is found with

$A[dB](\omega) = 20 ~ \log( \lvert G(\omega) \rvert )$

and thus the amplitude in [dB] for $K=1$ is 

$A[dB](0) = 20 ~ \log( 100 ) = 40.$

"

# ╔═╡ 6d87b94c-3fa6-11eb-0565-57aa96804bc3
G₁(ω) = 1/(sqrt(1 + ω^2))

# ╔═╡ 47bc2fbc-3fa7-11eb-2364-73007a336964
G₂(ω) = sqrt(100 + ω^2)

# ╔═╡ fcfa6c1e-3fa6-11eb-2e45-f5ed54428077
wspan = 0.1 : 0.1 : 100

# ╔═╡ dcbdc25c-3fa8-11eb-014a-25c0bd2679bc
K = 1

# ╔═╡ 74795452-3fa8-11eb-3fb8-891e82128d2b
Gol(ω) = K * G₁(ω)^3 * G₂(ω)^2

# ╔═╡ 73def6c8-3fa8-11eb-2049-a9239256c5e9
A(ω) = 20*log10(Gol(ω))

# ╔═╡ f3570e86-3fa8-11eb-3131-fdafb40ca2a5
begin
	mag = A(0);
	md"A(0) = $mag"
end

# ╔═╡ 3f610888-3fa7-11eb-001c-cfc5c0e9858c
plot(wspan, A.(wspan), xlabel="ω", ylabel="A [dB]", xaxis=:log, legend=false)

# ╔═╡ f007da02-3f02-11eb-3f4f-85c953ca9c02
md"## Phase

The phase of transfer function can be calculated with

$\phi(\omega) = \arctan\left( \frac{ \mathcal{Im}\{G(j\omega)\}}{\mathcal{Re}\{G(j\omega)\}}\right)$

and thus one holds for both subsystems

$\phi_{1}(\omega) = \arctan\left( -\omega \right) = - \arctan\left( \omega \right) \quad \text{and}$

$\phi_{2}(\omega) = \arctan\left( 0.1 ~ \omega \right).$

The phase of the open loop system is the sum of all phases as
$\phi_{OL}(\omega) = -3 ~ \arctan\left( \omega \right) + 2 ~ \arctan\left( 0.1 ~ \omega \right) \text{.}$

The phase at $\omega \rightarrow \infty$ can be simply be found with

$\lim\limits_{\omega \rightarrow \infty} \arctan(\omega) = \frac{\pi}{2}  = 90^{\circ}$

and thus the phase of the open loop at infinity is calculated with

$\lim\limits_{\omega \rightarrow \infty} \phi_{OL}(\omega) = -3 ~ \frac{\pi}{2} + 2 ~ \frac{\pi}{2}  = - \frac{\pi}{2} = -90^{\circ}.$
"

# ╔═╡ 79b10d9c-3f04-11eb-18b5-55316fc05dad
ϕ₁(ω) = -1*atan(ω)

# ╔═╡ 50fb5eac-3fa9-11eb-2f99-77e1de0e3bf9
ϕ₂(ω) = atan(0.1 * ω)

# ╔═╡ 639f9f1e-3fa9-11eb-0ccd-e1d989640d69
ϕol(ω) = 3*ϕ₁(ω) + 2*ϕ₂(ω)

# ╔═╡ 7f479da2-3fa9-11eb-3804-153dffb98355
begin
	w = 10000
	phase = ϕol(w);
	md"The phase is $\phi$($w) = $phase"
end

# ╔═╡ e29861e8-3fa9-11eb-1305-2fe747b1a3b4
plot(wspan, ϕol.(wspan), xlabel="ω", ylabel="ϕ", xaxis=:log, yticks=round.([-pi, -0.5*pi, 0], digits=2), legend=false)

# ╔═╡ 7074821e-413d-11eb-0e23-93dcfb7beb41
md"## Simulation with ControlSystems

Firstly, the transfer function is built as

$G(s) = \frac{s^2 + 20 s + 100}{s^3 + 3 s^2 + 3 s + 1}.$
"

# ╔═╡ 27bc0fe8-413c-11eb-313a-cf1309372fd3
begin
	import Pkg
	Pkg.add("ControlSystems")
end

# ╔═╡ 26095de0-3fac-11eb-3c8f-f37a9a0ef5fd
# Open loop transfer function 
G = tf([1, 20, 100], [1, 3, 3, 1])

# ╔═╡ 5a6b6ce8-413e-11eb-16e5-fb17f239638a
md"### Bode diagram

The magnitude starts $40 dB$ and the phase approaches $-90°$ for $\omega \rightarrow \infty$."

# ╔═╡ 08e2bfb2-3fb0-11eb-3625-0d489716b987
bodeplot(G, wspan, legend=false)

# ╔═╡ 122072e2-3fb0-11eb-2791-1105c9b8f8dc
setPlotScale("dB")

# ╔═╡ 05c0133a-4141-11eb-1b3f-4f72f08d78ad
md"## Phase margin

The magnitude $A[dB]$ is almost zero at $\omega=5$. The phase is given as $\phi(5) \approx -183° < -180°.$ Therefore, the phase margin $\Delta \phi = \phi(5) - (-180°) < 0°$  which means the system is **unstable**!"

# ╔═╡ 442abd22-5687-11eb-101a-17bbc4ce8e15
ϕ = -360.0 + bode(G, [5.0])[2][1,1,1] # A[dB] ≈ 0 -> ω=5 [1/s]

# ╔═╡ 307b4300-4142-11eb-124f-61b566dee570
md"## Stability

Due to the fact that the phase margin is almost zero, the closed-loop system

$G_{CL}(s = \frac{G_{OL}(s)}{1 + G_{OL}(s)}$

can be assumed as unstable. The impulse response of the closed-loop system unveils an unstable behaviour."

# ╔═╡ 259bfd0e-3fac-11eb-04e3-776e1e86940e
Gcl = G/(1 + G)

# ╔═╡ c6a79c30-5686-11eb-241a-1dc7c7d6cc75
Tf = 4.0 # Final simulation time

# ╔═╡ 6723d1ca-3fb6-11eb-0b03-65734298211a
impulseplot(Gcl, Tf, legend=false)

# ╔═╡ 2451f22c-4142-11eb-1eed-3d52a0501acb
md"## Stabilization

The closed-loop system can stabilized by varying the gain $K > 0$. If one desires a phase margin $\Delta \phi = 40°$, one has to find the frequency $\tilde{\omega}$ of the corresponding phase $\phi = -140°$."

# ╔═╡ 9358a562-3fb7-11eb-2ef0-978c3b3c785a
bode(G, [17])

# ╔═╡ de59a224-4146-11eb-13dc-31b2ed3964b2
md"The phase $\phi=-140°$ is found for $\omega=17$. The magnitude at this frequency is $A[db] = -20$. Therefore, one has a gain margin of $\Delta A[dB] = 20$. The gain in absolute values is calculated with

$V = 10^\frac{\Delta A[dB]}{20} = 10^\frac{20}{20} = 10.$"

# ╔═╡ d1c96570-3fb7-11eb-31ed-a10d6abcd120
V = 10^(20/20)

# ╔═╡ ea01876c-4147-11eb-2468-5bf379985237
md"The Bode plot for the dynamical system with gain $K=10$ is shown below."

# ╔═╡ 257fbfea-3fac-11eb-2cb6-2f4856d703eb
# Open loop transfer function 
G1 = V*tf([1, 20, 100], [1, 3, 3, 1])

# ╔═╡ 2568bea8-3fac-11eb-1f9f-6121fa90bac3
bodeplot(G1, wspan)

# ╔═╡ 56dd90ec-4148-11eb-00ad-5d56f856ae5a
md"At frequency $\omega=17$ the amplitude is $A[dB](17) \approx 0$ and the phase is $\phi(17) \approx -140°$. Therefore, the phase margin of $\Delta \phi = 40°$ is guaranteed."

# ╔═╡ b76c67d2-3fb6-11eb-3dab-97f754ade630
bode(G1, [17])

# ╔═╡ eb0e05a8-4148-11eb-0fee-b10f1422bde0
md"The stability of the closed-loop system is guaranteed because of the phase margin $\Delta \phi > 0$. The simulation of the impulse response below shows a stable behaviour."

# ╔═╡ 2550442c-3fac-11eb-0a09-bda85235be5a
Gcl1 = G1/(1 + G1)

# ╔═╡ 252501fe-3fac-11eb-315e-0bd6986dbec5
impulseplot(Gcl1, 3.0)

# ╔═╡ 4d98c0be-4149-11eb-29ea-ed09bc92feed
md"## Steady-state error

There exist a (small or tiny) steady-state error because the open-loop system does **not** contain an integrator element $\frac{1}{s}$. This can be either noticed in the transfer function or in the phase plot. If a system contains one integrator element then the phase will always start at $\phi(0) = -90°$.

Due to the large gain $V=10 > 1$ the steady-state error is neglectable in this case." 

# ╔═╡ f4bf5516-4148-11eb-1e4c-5fbb7df849cc
stepplot(Gcl1, Tf, legend=false)

# ╔═╡ Cell order:
# ╟─95230452-3fa6-11eb-114b-e300526fdffc
# ╟─2c3f0c66-3efb-11eb-1c40-9555f67871b6
# ╠═6d87b94c-3fa6-11eb-0565-57aa96804bc3
# ╠═47bc2fbc-3fa7-11eb-2364-73007a336964
# ╠═9502a1d8-3fa6-11eb-187d-2dff3e492226
# ╠═fcfa6c1e-3fa6-11eb-2e45-f5ed54428077
# ╠═dcbdc25c-3fa8-11eb-014a-25c0bd2679bc
# ╠═74795452-3fa8-11eb-3fb8-891e82128d2b
# ╠═73def6c8-3fa8-11eb-2049-a9239256c5e9
# ╠═f3570e86-3fa8-11eb-3131-fdafb40ca2a5
# ╠═3f610888-3fa7-11eb-001c-cfc5c0e9858c
# ╟─f007da02-3f02-11eb-3f4f-85c953ca9c02
# ╠═79b10d9c-3f04-11eb-18b5-55316fc05dad
# ╠═50fb5eac-3fa9-11eb-2f99-77e1de0e3bf9
# ╠═639f9f1e-3fa9-11eb-0ccd-e1d989640d69
# ╠═7f479da2-3fa9-11eb-3804-153dffb98355
# ╠═e29861e8-3fa9-11eb-1305-2fe747b1a3b4
# ╟─7074821e-413d-11eb-0e23-93dcfb7beb41
# ╠═27bc0fe8-413c-11eb-313a-cf1309372fd3
# ╠═261fff78-3fac-11eb-1e00-bb85b363d664
# ╠═26095de0-3fac-11eb-3c8f-f37a9a0ef5fd
# ╟─5a6b6ce8-413e-11eb-16e5-fb17f239638a
# ╠═08e2bfb2-3fb0-11eb-3625-0d489716b987
# ╠═122072e2-3fb0-11eb-2791-1105c9b8f8dc
# ╟─05c0133a-4141-11eb-1b3f-4f72f08d78ad
# ╠═442abd22-5687-11eb-101a-17bbc4ce8e15
# ╟─307b4300-4142-11eb-124f-61b566dee570
# ╠═259bfd0e-3fac-11eb-04e3-776e1e86940e
# ╠═c6a79c30-5686-11eb-241a-1dc7c7d6cc75
# ╠═6723d1ca-3fb6-11eb-0b03-65734298211a
# ╟─2451f22c-4142-11eb-1eed-3d52a0501acb
# ╠═9358a562-3fb7-11eb-2ef0-978c3b3c785a
# ╟─de59a224-4146-11eb-13dc-31b2ed3964b2
# ╠═d1c96570-3fb7-11eb-31ed-a10d6abcd120
# ╟─ea01876c-4147-11eb-2468-5bf379985237
# ╠═257fbfea-3fac-11eb-2cb6-2f4856d703eb
# ╠═2568bea8-3fac-11eb-1f9f-6121fa90bac3
# ╟─56dd90ec-4148-11eb-00ad-5d56f856ae5a
# ╠═b76c67d2-3fb6-11eb-3dab-97f754ade630
# ╟─eb0e05a8-4148-11eb-0fee-b10f1422bde0
# ╠═2550442c-3fac-11eb-0a09-bda85235be5a
# ╠═252501fe-3fac-11eb-315e-0bd6986dbec5
# ╟─4d98c0be-4149-11eb-29ea-ed09bc92feed
# ╠═f4bf5516-4148-11eb-1e4c-5fbb7df849cc
