### A Pluto.jl notebook ###
# v0.12.12

using Markdown
using InteractiveUtils

# ╔═╡ ff6d01d2-5fc2-11eb-1ea1-f7bd3c807b18
using ControlSystems

# ╔═╡ fa48723e-5fe4-11eb-3ab3-03111ca9933f
using Plots

# ╔═╡ 04ac48e2-5ff0-11eb-0b70-498de97ecd94
md"# Ziegler-Nichols method

## Original System
"

# ╔═╡ 02f23862-5fe2-11eb-2ced-c527d05c2f04
Gp = tf(0.125, [1, 1.5, 0.75, 0.125])

# ╔═╡ 02d3cc60-5fe2-11eb-111a-e3c0010f96a4
Tf = 20.0 # Simulation time

# ╔═╡ 01a6d71c-5fe2-11eb-141f-2d761892bc6f
stepplot(Gp, Tf, legend=false)

# ╔═╡ 47ced8c4-5fe6-11eb-21af-21c560fab8af
begin
	amp = step(Gp, Tf)[1]
	tsteps = step(Gp, Tf)[2]
end

# ╔═╡ 589cb5b8-5fe6-11eb-2e5a-75022f54d759
begin
	approx_grad = (amp[75] - amp[65])/(tsteps[75] - tsteps[65])
	approx_offset = amp[70] -approx_grad * tsteps[70]  
	approx = approx_grad * tsteps .+ approx_offset
end

# ╔═╡ 8e8d41d8-5fe6-11eb-1758-29fa986af146
plot(tsteps[1:150], [approx[1:150], amp[1:150]], legend=false)

# ╔═╡ 89cfeb04-5fe7-11eb-07fa-f76673f99ac0
t0 = tsteps[23]

# ╔═╡ 8f9d3f00-5fe7-11eb-1d54-d1466c68a861
f0 = approx[23]

# ╔═╡ 53323364-5fe8-11eb-2b1f-e5feab0943e6
t1 = tsteps[133]

# ╔═╡ 3a7ca3d4-5fe8-11eb-2f9c-458c5b2d9fb4
f1 = approx[133]

# ╔═╡ ae331d88-5fef-11eb-2c4c-bd6b90d0730d
Td = t0;

# ╔═╡ d90eefa2-5fe8-11eb-1f00-917b4589e812
Tg = t1 - t0

# ╔═╡ bf848976-5fe7-11eb-31d6-adb3eedb1715
md"The tangiant intersects
- the zero line at t = $t0 : f( $t0 ) = $f0
- the reference at t = $t1 : f( $t1 ) = $f1

Thus, the time constant are given as 
- delay $T_{d} =$ $t0 seconds and 
- gradient $T_{g} =$ $t1 - $t0 = $Tg seconds.

## Approximated first-order system

The general approximated first-order system with delay is noted as

$\tilde{G}_{p}(s) = \frac{K}{1 + T_{g} ~ s} e^{-T_{d}~s}$

and in this case one holds

$\tilde{G}_{p}(s) = \frac{1}{1 + 7.7 ~ s} e^{-1.54~s}.$

"

# ╔═╡ 74e34678-5fe7-11eb-26fd-296a0b7b6c34
Gapprox = tf(1, [Tg, 1])*delay(Td)

# ╔═╡ aa2bd1a0-5fe6-11eb-326f-a9f513e86d44
stepplot([Gp, Gapprox], 3*Tf, label=["original" "approx"])

# ╔═╡ c739abfc-5fe6-11eb-3469-35988cb11800
md"## Controller design with Ziegler-Nichols approach

The coefficients for each PID controller type are defined as 

- Proportional controller $G_{c}(s) = K_{p}$

$K_{p} = \frac{T_{g}}{T_{d}}$

- PI controller $G_{c}(s) = K_{p} + \frac{1}{T_{I} s} = K_{p} \left(1 + \frac{1}{T_{n} s} \right)$

$K_{p} = 0.9 ~ \frac{T_{g}}{T_{d}} \quad \text{and} \quad T_{n} = 3.3 ~ T_{d}$

- PID controller $G_{c}(s) = K_{p} + \frac{1}{T_{I} s} + T_{D} s = K_{p} \left(1 + \frac{1}{T_{n} s} + T_{v} s \right)$

$K_{p} = 1.2 ~ \frac{T_{g}}{T_{d}} \quad \text{and} \quad T_{n} = 2 ~ T_{d} \quad \text{and}  \quad T_{v} = 0.5 ~ T_{d}$

Here, a PI controller is desired and its coefficients are

$K_{p} = 0.9 \cdot \frac{7.7}{1.54} = 4.5 \quad \text{and} \quad T_{n} = 3.3 \cdot 1.54 \approx 5.1$

and thus one holds

$G_{c}(s) = 4.5 \left(1 + \frac{1}{5.1 ~ s} \right) \approx 4.5 \left(1 + \frac{0.2}{s}\right) = \frac{4.5~s + 0.9}{s}.$

"

# ╔═╡ e493e940-5fe3-11eb-2e4f-bd08a44647b3
Kp = 0.9*(Tg/Td) # Proportional gain

# ╔═╡ 1989140e-5fee-11eb-0b15-b9bc0e967efc
Tn = 3.3 * t0 # Integrator time constant

# ╔═╡ 24f4cc66-5fee-11eb-14f1-419ae3ee8a51
Ki = Kp/Tn # Integrator gain

# ╔═╡ 364b1042-5fe4-11eb-1c71-837c01e7a808
Gc = tf([Kp, Ki],[1, 0]) # PI controller

# ╔═╡ b5afad5c-5fee-11eb-01cb-69b7f416b1b3
md"### Feedback simulation of approximated system"

# ╔═╡ 134d6d3c-5fe5-11eb-3228-cf8ef8daa608
Gol_approx = Gc * Gapprox; # Open-Loop for approximated system 

# ╔═╡ 8a4c7ffa-5fee-11eb-0f50-8910d73a53cd
Tf2 = 50; # Simulation time

# ╔═╡ 6990b04e-5fe5-11eb-0697-9b8c13d6cf80
stepplot(feedback(Gol_approx), Tf2, legend=false)

# ╔═╡ c335920c-5fee-11eb-1c6a-4366d8eda15f
md"### Feedback simulation of original system"

# ╔═╡ cb14bcdc-5fee-11eb-0789-230f4ff81ef5
Gol = Gc * Gp; # Open-Loop for approximated system 

# ╔═╡ e99aa446-5fee-11eb-3173-33ae4a5a1bc3
stepplot(feedback(Gol), Tf2, legend=false)

# ╔═╡ Cell order:
# ╟─04ac48e2-5ff0-11eb-0b70-498de97ecd94
# ╠═ff6d01d2-5fc2-11eb-1ea1-f7bd3c807b18
# ╠═02f23862-5fe2-11eb-2ced-c527d05c2f04
# ╠═02d3cc60-5fe2-11eb-111a-e3c0010f96a4
# ╠═01a6d71c-5fe2-11eb-141f-2d761892bc6f
# ╠═47ced8c4-5fe6-11eb-21af-21c560fab8af
# ╠═589cb5b8-5fe6-11eb-2e5a-75022f54d759
# ╠═fa48723e-5fe4-11eb-3ab3-03111ca9933f
# ╠═8e8d41d8-5fe6-11eb-1758-29fa986af146
# ╠═89cfeb04-5fe7-11eb-07fa-f76673f99ac0
# ╠═8f9d3f00-5fe7-11eb-1d54-d1466c68a861
# ╠═53323364-5fe8-11eb-2b1f-e5feab0943e6
# ╠═3a7ca3d4-5fe8-11eb-2f9c-458c5b2d9fb4
# ╠═ae331d88-5fef-11eb-2c4c-bd6b90d0730d
# ╠═d90eefa2-5fe8-11eb-1f00-917b4589e812
# ╟─bf848976-5fe7-11eb-31d6-adb3eedb1715
# ╠═74e34678-5fe7-11eb-26fd-296a0b7b6c34
# ╠═aa2bd1a0-5fe6-11eb-326f-a9f513e86d44
# ╟─c739abfc-5fe6-11eb-3469-35988cb11800
# ╠═e493e940-5fe3-11eb-2e4f-bd08a44647b3
# ╠═1989140e-5fee-11eb-0b15-b9bc0e967efc
# ╠═24f4cc66-5fee-11eb-14f1-419ae3ee8a51
# ╠═364b1042-5fe4-11eb-1c71-837c01e7a808
# ╟─b5afad5c-5fee-11eb-01cb-69b7f416b1b3
# ╠═134d6d3c-5fe5-11eb-3228-cf8ef8daa608
# ╠═8a4c7ffa-5fee-11eb-0f50-8910d73a53cd
# ╠═6990b04e-5fe5-11eb-0697-9b8c13d6cf80
# ╟─c335920c-5fee-11eb-1c6a-4366d8eda15f
# ╠═cb14bcdc-5fee-11eb-0789-230f4ff81ef5
# ╠═e99aa446-5fee-11eb-3173-33ae4a5a1bc3
