### A Pluto.jl notebook ###
# v0.12.12

using Markdown
using InteractiveUtils

# ╔═╡ 1d620648-6544-11eb-1dd9-732308a21a8d
using ControlSystems

# ╔═╡ a50903e6-6559-11eb-1eaa-e1955e182775
md"# Root locus design

Consider the transfer function of the open-loop system

$G_{OL}(s) = \frac{K}{s(s + 7.5)(s + 12.5)}.$

This system is **marginally stable** because it has the poles $\{-12.5, -7.5, 0\}$.

"

# ╔═╡ 23309a44-6544-11eb-2c25-d1f2ef6f4cb1
Gol = tf(1, [1, 20, 7.5*12.5, 0])

# ╔═╡ 4e6b1ff4-6544-11eb-21bb-8d3c0b937518
rlocusplot(Gol, K=10000, legend=:topleft)

# ╔═╡ 7327160a-6570-11eb-2984-518ab9baf102
md"## Construction of the root locus

1. All poles of the closed-loop system start ($K=0$) at the poles of open-loop system:

$\{-12.5, -7.5, 0\}$

2. The poles of the closed-loop system proceed to the zeros of the open-loop system: here none

3. The roots leave the real axis at the bifurcation point which is calculated by

$\frac{d}{ds} G_{OL}(s) = \frac{d}{ds} \left( \frac{N(s)}{D(s)} \right) = \frac{N'(s)~D(s) + N(s)~D'(s)}{D(s)^2} = 0.$

This is calculated here with $N(s) = K$, $N'(s) = 0$,

$D(s) = s(s + 7.5)(s + 12.5) = s^3 + 20~s^2 + 93.75~s$

and 

$D'(s) = 3~s^2 + 40~s + 93.75.$

Consequently, one yields

$N(s) D'(s) = K (3~s^2 + 40~s + 93.75) = 0$

and further

$s_{1,2} = \frac{-40}{6} \pm \frac{1}{6} \sqrt{40^2 - 12 \cdot 93.75} \approx \{-10, -3\}$"

# ╔═╡ a5084c96-6559-11eb-27a1-03bee116d2de
s₁ = -(40/6) + sqrt(40^2 - 12*93.75)/6

# ╔═╡ a4c4084c-6559-11eb-1b9f-f76e5ae691a5
s₂ = -(40/6) - sqrt(40^2 - 12*93.75)/6

# ╔═╡ 3e7bc014-6560-11eb-3946-2b06bfc2ee30
md"Only $s \approx -3.0$ is a valid bifurcation point because it has to be a pole of the closed-loop system. This is proved with

$G_{CL}(s) = \frac{G_{OL}(s)}{1 + G_{OL}(s)}$

and

$1 + G_{OL}(s) = 1 + \frac{K}{s(s + 7.5)(s + 12.5)} = 0 \quad \text{or} \quad K = -s~(s + 7.5)~(s + 12.5).$"

# ╔═╡ 3ed2a3d8-6561-11eb-2b0a-d33ff63d959b
K(s) = -s * (s + 7.5) * (s + 12.5)

# ╔═╡ 4c4fdbac-6561-11eb-21c7-259119f2babb
K₁ = K(s₁) # Gain of bifurcation 

# ╔═╡ 5bd2d55c-6561-11eb-0bba-677d3a064225
K₂ = K(s₂) 

# ╔═╡ 716583a6-6561-11eb-1285-2312faed82e3
md"The gain has to be greater zero. Therefore, only s=-3 is a bifurcation point.

4. The angle at the bifurcation point is found with 

$\phi = \frac{180°}{r} = \frac{180°}{2} = 90°$

where $r$ is the number of branches.


5. The root locus branches intersect the imaginary axis at $(\rho, \omega)$. The intersection point is found with

$1 + G_{OL}(s) = 1 + \left.\frac{N(s)}{D(s)}\right\rvert_{s=j\omega} = 0 \quad \text{or} \quad N(j~\omega) + D(j~\omega) = 0.$

Here, this point is calculated with

$N(j~\omega) + D(j~\omega) = K + (-j \omega^3) - 20~\omega^2 + 93.75~j \omega = K - 20 \omega^2 + j \omega (93.75 - \omega^2)$

and separated in real and imaginary part as
- Imaginary: $\omega (93.75 - \omega^2) = 0$ thus $\omega = \sqrt{93.75} \approx 9.7$ and
- Real: $K - 20~\omega^2 = 0$ which is equivalent $K = 20 \cdot 93.75 = 1875.$

"

# ╔═╡ 3a230be0-656e-11eb-0750-4983952ef4c2
ω = sqrt(93.75)

# ╔═╡ b807246a-656e-11eb-3720-51366e643b95
Kc = 20*93.75 # Critical gain 

# ╔═╡ 3a01739c-656e-11eb-2844-3762d6366bfd
md"6. The root locus branches converge to linear asymptotes. These asymptotes start at point

$s_{a} = \frac{\sum\limits_{i=1}^{n} s_{P,i} - \sum\limits_{j=1}^{m} s_{Z,j}}{n - m}$

where $n$ is the number of poles $s_{P,i}$ and $m$ is the number of zeros $s_{Z,j}$.

Here, the starting point is calculated as

$s_{a} = \frac{-12.5 - 7.5 - 0}{3 - 0} = \frac{-20}{3} \approx -6.7.$

7. The angular between the real axis and the asymptotes are calculated with

$\phi_{A,i} = \frac{(2 ~ i + 1) \cdot 180°}{n - m} = \frac{2 ~ i + 1}{3} \cdot 180° \quad \text{for} \quad i=0, 1, \cdots, n-m-1.$

Therefore, the angles are  

$\phi_{A,0} = \frac{180°}{3} = 60°  \quad \text{,} \quad  \phi_{A,1} = \frac{3 \cdot 180°}{3} = 180° \quad \text{and} \quad \phi_{A,2} = \frac{5 \cdot 180°}{3} = 300°.$"

# ╔═╡ d8cbf566-6570-11eb-3f5a-ef301bf410a6
rlocusplot(Gol, K=1000,xticks=-4:0.5:0.5,yticks = -6:0.5:6,xlims=(-4.0, 0.5),legend=:topleft)

# ╔═╡ dae9a748-6571-11eb-11ac-dbdc99f5a4fa
md"### Damping / Attenuation

If an attenuation 

$\zeta = 0.7 = \sin(\alpha)$

is assumed, then the corresponding angle is calculated with

$\alpha[°] = \arcsin(0.7) ~ \frac{180°}{\pi} \approx 45°.$

Drawing a line from the origin with an angle of $45°$ the point $s = -2.7 + 2.7j$ on the branch is found. 

"

# ╔═╡ 01b90156-6573-11eb-038f-0bcbbc723ea5
s₀ = -2.7 + 2.7im

# ╔═╡ 11d4e71e-6573-11eb-29fb-7b9c419c9c9c
K₀ = real(K(s₀))

# ╔═╡ Cell order:
# ╟─a50903e6-6559-11eb-1eaa-e1955e182775
# ╠═1d620648-6544-11eb-1dd9-732308a21a8d
# ╠═23309a44-6544-11eb-2c25-d1f2ef6f4cb1
# ╠═4e6b1ff4-6544-11eb-21bb-8d3c0b937518
# ╟─7327160a-6570-11eb-2984-518ab9baf102
# ╠═a5084c96-6559-11eb-27a1-03bee116d2de
# ╠═a4c4084c-6559-11eb-1b9f-f76e5ae691a5
# ╟─3e7bc014-6560-11eb-3946-2b06bfc2ee30
# ╠═3ed2a3d8-6561-11eb-2b0a-d33ff63d959b
# ╠═4c4fdbac-6561-11eb-21c7-259119f2babb
# ╠═5bd2d55c-6561-11eb-0bba-677d3a064225
# ╟─716583a6-6561-11eb-1285-2312faed82e3
# ╠═3a230be0-656e-11eb-0750-4983952ef4c2
# ╠═b807246a-656e-11eb-3720-51366e643b95
# ╟─3a01739c-656e-11eb-2844-3762d6366bfd
# ╠═d8cbf566-6570-11eb-3f5a-ef301bf410a6
# ╟─dae9a748-6571-11eb-11ac-dbdc99f5a4fa
# ╠═01b90156-6573-11eb-038f-0bcbbc723ea5
# ╠═11d4e71e-6573-11eb-29fb-7b9c419c9c9c
