### A Pluto.jl notebook ###
# v0.12.12

using Markdown
using InteractiveUtils

# ╔═╡ 58dca28c-615c-11eb-21f2-87d0713cf56c
using ControlSystems

# ╔═╡ cbf89b1a-615d-11eb-1578-53ea1f1e6841
md"# Ziegler-Nichols method: Example

## System approximation
"

# ╔═╡ 4626b6b4-615c-11eb-3a2e-a90663dc91e3
Td = 2.5 # Time delay

# ╔═╡ 4fb8fc78-615c-11eb-16df-8f1db543a53b
Tg = 13.5 - Td # Time constant

# ╔═╡ 620274fe-615c-11eb-3742-573a7d402e09
Gapprox = tf(1, [Tg, 1])*delay(Td)

# ╔═╡ 86c8d2e2-615c-11eb-14e7-7f6e98b514f7
Tf = 60 # Simulation time

# ╔═╡ 84f783dc-615c-11eb-17e0-d1735e6548a2
stepplot(Gapprox, Tf, legend=false)

# ╔═╡ 11156da2-615d-11eb-0073-89c995c2e330
md"## Controller design

PI controller 

$G_{c}(s) = K_{p} + \frac{1}{T_{I} s} = K_{p} \left(1 + \frac{1}{T_{n} s} \right)$

Here, the coefficients are calculated as

$K_{p} = 0.9 \cdot \frac{11}{2.5} \approx 4.0 \quad \text{and} \quad T_{n} = 3.3 \cdot 2.5 \approx 8.25$

and thus one holds

$G_{c}(s) = 4.0 \left(1 + \frac{1}{8.25 ~ s} \right) \approx 4.0 \left(1 + \frac{0.13}{s}\right) = \frac{4.0~s + 0.5}{s}.$

"

# ╔═╡ 9c97be76-615c-11eb-0ed2-3b04fb02c701
Kp = 0.9*(Tg/Td) # Proportional gain

# ╔═╡ a6bb72ee-615c-11eb-3a88-3379f0d259ff
Tn = 3.3 * Td # Integrator time constant

# ╔═╡ ac60052a-615c-11eb-0dfe-3d4423340791
Ki = Kp/Tn # Integrator gain

# ╔═╡ ac2afdd0-615c-11eb-18d2-17c7e7dfbfe2
Gc = tf([Kp, Ki],[1, 0]) # PI controller

# ╔═╡ e9a912d2-615c-11eb-165f-f399e656be58
Gol_approx = Gc * Gapprox; # Open-Loop for approximated system 

# ╔═╡ ac08deb2-615c-11eb-1322-cbf4ddd9fc58
stepplot(feedback(Gol_approx), Tf, legend=false)

# ╔═╡ abe01ab8-615c-11eb-02a4-a90f75e1d389


# ╔═╡ ab9685b2-615c-11eb-1e6a-61a2e16e6aa4


# ╔═╡ ab2561aa-615c-11eb-3359-71adfd4cc04b


# ╔═╡ Cell order:
# ╟─cbf89b1a-615d-11eb-1578-53ea1f1e6841
# ╠═4626b6b4-615c-11eb-3a2e-a90663dc91e3
# ╠═4fb8fc78-615c-11eb-16df-8f1db543a53b
# ╠═58dca28c-615c-11eb-21f2-87d0713cf56c
# ╠═620274fe-615c-11eb-3742-573a7d402e09
# ╠═86c8d2e2-615c-11eb-14e7-7f6e98b514f7
# ╠═84f783dc-615c-11eb-17e0-d1735e6548a2
# ╟─11156da2-615d-11eb-0073-89c995c2e330
# ╠═9c97be76-615c-11eb-0ed2-3b04fb02c701
# ╠═a6bb72ee-615c-11eb-3a88-3379f0d259ff
# ╠═ac60052a-615c-11eb-0dfe-3d4423340791
# ╠═ac2afdd0-615c-11eb-18d2-17c7e7dfbfe2
# ╠═e9a912d2-615c-11eb-165f-f399e656be58
# ╠═ac08deb2-615c-11eb-1322-cbf4ddd9fc58
# ╠═abe01ab8-615c-11eb-02a4-a90f75e1d389
# ╠═ab9685b2-615c-11eb-1e6a-61a2e16e6aa4
# ╠═ab2561aa-615c-11eb-3359-71adfd4cc04b
