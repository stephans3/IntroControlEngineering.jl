### A Pluto.jl notebook ###
# v0.12.12

using Markdown
using InteractiveUtils

# ╔═╡ f0649f3c-7b70-11eb-241f-b3dbc878e572
using Luxor

# ╔═╡ 88b19370-7b72-11eb-333c-336da378ac9e
using Plots

# ╔═╡ 5e4b0804-7b73-11eb-30dd-53f5fe5e80c0
using ControlSystems

# ╔═╡ 39a4a238-7b3e-11eb-0f9b-71dae147e0bb
md"# Electrical circuit

Consider a circuit with a voltage source $u_{in}(t)$, a resistor

$u_{R}(t) = R ~ i(t)$

with $ R > 0$, an inductor 

$u_{L}(t) = L ~ \frac{d i(t)}{dt}$

with $L > 0$, and a capacitor 

$u_{C}(t) = \frac{1}{C} ~ \int i(\tau) ~ d\tau$

with $C > 0$ at the output. The components are connected serial as portrayed in the figure below and one holds 

$u_{in}(t) = u_{R}(t) + u_{L}(t) + u_{C}(t) \text{.} \tag{1}$

Equation $(1)$ is noted in terms of output $u_{C}$ with $i(t) = C ~ \dot{u}_{C}(t)$
as

$u_{in}(t) = R ~ C ~ \dot{u}_{C}(t) + L ~ C ~ \ddot{u}_{C}(t) + u_{C}(t)  \text{.}  \tag{2}$

The second order differential equation $(2)$ can be noted as a first order system 

$\begin{pmatrix}
\dot{x}_{1}(t) \\
\dot{x}_{2}(t)
\end{pmatrix}
=
\begin{pmatrix}
0 & 1 \\
- \frac{1}{L C} & - \frac{R}{L}
\end{pmatrix}
\begin{pmatrix}
x_{1}(t) \\
x_{2}(t)
\end{pmatrix}
+
\begin{pmatrix}
0 \\
\frac{K}{L C}
\end{pmatrix}
u(t)$

with $x_{1}(t) = u_{C}(t)$, $x_{2}(t) = \dot{u}_{C}(t)$ and $u_{in}(t) = K ~ u(t)$. This electrical oscillator works analog to the mechanical oscillator and can be simulated as discussed before.

Otherwise, Equation $(2)$ can be noted with the Laplace transform (see [Wikipedia](https://en.wikipedia.org/wiki/Laplace_transform)) in the frequency domain as

$U_{in}(s) = s ~ R ~ C ~ U_{c}(s) + s^2 ~ L ~ C ~ U_{c}(s) + U_{c}(s) = \left[ s^2 ~ L C + s ~ R C + 1 \right] U_{C}(s)$

and rewritten as the transfer function

$G(s) = \frac{Y(s)}{U(s)} = \frac{K}{s^2 ~ L C + s ~ R C + 1} = \frac{K}{s^2 ~ T_{1}^2 + s ~ T_{2} + 1}$

with $U_{in}(s) = K ~ U(s)$ (U is called input), $U_{C}(s) = Y(s)$ (Y is called output), $T_{1} = \sqrt{L C}$ and $T_{2} = R C$. 

Using the transfer function $G$ the behaviour of output $Y$ can be tested for steps and impulses at the input $U$ - so called step and impulse responses."


# ╔═╡ 7b313956-7b42-11eb-0c28-97adc530d776
@svg begin
	
	# Nodes
	circle( Point(-200,-70), 5, :stroke)
	circle( Point(200,-70), 5, :stroke)
	circle( Point(-200,70), 5, :stroke)
	circle( Point(200,70), 5, :stroke)
	rect(Point(-150, -85), 80, 30, :stroke)
	
	rad = 12
	# Inductor
	for i = 1 : 4
	Luxor.arrow(Point(-30+2*i*rad, -70), rad, π, 0, arrowheadlength=0, arrowheadangle=pi/12, linewidth=2.0)	
	end
	
	# Capacitor
	line( Point(125, -10), Point(155, -10), :stroke)
	line( Point(125, 10), Point(155, 10), :stroke)
	
	# Edges
	line( Point(-195,-70), Point(-150,-70), :stroke)
	line( Point(-70,-70), Point(-30+rad+1,-70), :stroke)
	line( Point(9*rad-30,-70), Point(195,-70), :stroke)
	line( Point(-195, 70), Point(195, 70), :stroke)
	line( Point(140, -70), Point(140,-10), :stroke )
	line( Point(140, 10), Point(140, 70), :stroke )
	
	# Voltages
	Luxor.arrow( Point(-150, -100), Point(-70, -100), linewidth=2.0)
	Luxor.arrow( Point(-30+rad+1, -100), Point(9*rad-30, -100), linewidth=2.0)
	Luxor.arrow( Point(220, -65), Point(220, 65), linewidth=2.0 )
	Luxor.arrow( Point(-220, -65), Point(-220, 65), linewidth=2.0 )
	
	# Labels
	fontsize(30)
	Luxor.text("uᵢₙ", Point(-265, 10))
	Luxor.text("u", Point(-120, -110))
	Luxor.text("u", Point(20, -110))
	Luxor.text("u", Point(235, 10))
	Luxor.text("R", Point(-120, -25))
	Luxor.text("L", Point(20, -25))
	Luxor.text("C", Point(90, 10))
	
	
	fontsize(20)
	Luxor.text("R", Point(-100, -105))
	Luxor.text("L", Point(40, -105))
	Luxor.text("C", Point(255, 20))
end

# ╔═╡ 91416dc8-7b3e-11eb-0d17-9bf55cfca0b5
md"## Step and impulse response

Both responses can be explained by a simple electrical circuit with a resistor and capacitor (RC oscillator) without inductor. So, induction $L = 0$ and one holds 

$u_{in}(t) = R C ~ \dot{u}_{C}(t) + u_{C}$

or rewritten with $y(t) = u_{C}$, $u(t) = u_{in}$ and $a = 1/(R C)$ as a linear ODE

$\dot{y}(t) = -a ~ y(t) + a ~ u(t) \text{.} \tag{3}$

### Impulse response

If input $u(t)$ is an impulse, for example 

$u(t) = \delta(t) = \begin{cases} 1 \quad \text{for} \quad t = 0 \\ 0 \quad \text{else}\end{cases}$ 

then one yields an impulse response $g(t)$ as the behaviour of the output $y(t) = g(t)$. One can assume for the simple ODE $(3)$ that $u(t) = 0$ and $y(0) = 1 \cdot a$ and thus 

$y(t) = g(t) = a ~ exp(-a ~t) \text{.}$

The resulting graph is plotted below."

# ╔═╡ 2748ae28-7b58-11eb-273a-c3060406b035
a = 0.8 

# ╔═╡ 88cd908e-7b72-11eb-0852-75d7ec7f129b
t = 0.0 : 0.01 : 10; # Time range

# ╔═╡ 8896476e-7b72-11eb-362b-a9575e01485a
plot(t, t -> a * exp(-a*t), label="y= g(t)")

# ╔═╡ d9f5d3e0-7b72-11eb-0054-eb985c9c1faa
md"### Step response

If input $u(t)$ is a step, for example $u(t) = 1 ~ \text{for} ~ t \geq 0$, then one yields an step response $h(t)$ as the behaviour of the output $y(t) = h(t)$. 

The linear ODE $\dot{y}(t) = -a ~ \left[ y(t) - 1 \right] $ with the initial value $y(t = 0) = x_{0} = 0$ is solved via integration as

$\ln\left( \frac{y - 1}{y_{0} - 1} \right) = -a ~ t$

and further

$y(t) = 1 + [x_{0} - 1] ~ \exp(-a ~ t) = 1 - \exp(-a ~ t) = h(t) \text{.}$

One observes, that $g(t) = a ~ \exp(-a ~ t) = \frac{d}{dt} h(t)$. This holds in general, too, and leads with the Laplace transform in the frequency domain to $G(s) = s ~ H(s)$.

The resulting graphs for the impulse $y(t) = g(t)$ and step response $y(t) = h(t)$ are plotted below."

# ╔═╡ f21f22de-7b72-11eb-1519-d7080c3cf5f1
plot(t, [t -> a*exp(-a*t), t-> 1-exp(-a*t)], label=["y = g(t)" "y = h(t)"])

# ╔═╡ 2c6b86d8-7b73-11eb-3e8c-2daf70acbee4
md"## Dynamics of the electrical oscillator

The transfer function of the electrical oscillator is given as

$G(s) = \frac{Y(s)}{U(s)} = \frac{K}{s^2 ~ T_{1}^2 + s ~ T_{2} + 1} = \frac{K \omega_{1}^2 }{s^2  +  \omega_{2} ~ s + \omega_{1}^2}$

with frequencies $\omega_{1} = T_{1}^{-1}$, $\omega_{2} = T_{2}/T_{1}^{2}$ and the output can be calculated with $Y(s) = G(s) ~ U(s)$.

### Stability

The transfer function $G(s)$ is called bounded-input-bounded-output (BIBO) stable if and only if the real part of all poles are smaller than zero. 

Thus, the denominator is studied for

$s^2  + \omega_{2} ~ s + \omega_{1}^2 = 0$

and this leads to 
$s_{1,2} = \frac{-\omega_{2}}{2} \pm \frac{1}{2} \sqrt{\omega_{2}^2 - 4 \omega_{1}^2} \text{.}$

Therefore, if $\omega_{1} > 0$ and $\omega_{2} > 0$ than the transfer function $G(s)$ is BIBO stable. This condition is equal to the stability proof with eigenvalues in the state space representation. 

### Impulse response

The impulse $u(t) = \delta(t)$ is Laplace transformed to $U(s) = 1$ and thus $Y(s) = G(s) ~ U(s) = G(s)$. The output in the time domain can be either derived by partial fraction decomposition and inverse Laplace transform or simulated with specific software libraries.
"

# ╔═╡ 2c35a19e-7b73-11eb-0397-e5b698d004c7
T₁ = 1.0;

# ╔═╡ 5e97b884-7b73-11eb-0f54-1d088b792897
T₂ = 0.6;

# ╔═╡ 5e7f4b14-7b73-11eb-1a73-bd7c68ec5807
K = 2.3; # Gain or amplification

# ╔═╡ 5e6590fc-7b73-11eb-1815-636463bb5f9a
Tf = 20; # Final time for simulation

# ╔═╡ 2bedb2bc-7b73-11eb-1f98-6d62cdcc853b
G = tf([K], [T₁^2, T₂, 1])

# ╔═╡ 5e2fdb88-7b73-11eb-3bb1-1ff54b515847
impulseplot(G, Tf, label="y(t) = g(t)")

# ╔═╡ 835a03f0-7b73-11eb-1706-d98a60bba904
md"### Step response

The step response $y(t) = h(t) = \int g(\tau) d\tau$ is Laplace transformed to $Y(s) = H(s) = \frac{1}{s} G(s)$. Thus, the input in the frequency domain is given as $U(s) = \frac{1}{s}$ and one yields 

$Y(s) = G(s) ~ U(s) ~=~ \frac{K \omega_{1}^2 }{s^3  +  \omega_{2} ~ s^2 + \omega_{1}^2 ~ s} \text{.}$

The step response can be calculated in the time domain either analytically or numerically analog to the impulse response."

# ╔═╡ a3d0e1a0-7b73-11eb-0a20-19d2921060c0
stepplot(G, Tf, label="y(t) = h(t)")

# ╔═╡ Cell order:
# ╟─39a4a238-7b3e-11eb-0f9b-71dae147e0bb
# ╠═f0649f3c-7b70-11eb-241f-b3dbc878e572
# ╠═7b313956-7b42-11eb-0c28-97adc530d776
# ╟─91416dc8-7b3e-11eb-0d17-9bf55cfca0b5
# ╠═2748ae28-7b58-11eb-273a-c3060406b035
# ╠═88cd908e-7b72-11eb-0852-75d7ec7f129b
# ╠═88b19370-7b72-11eb-333c-336da378ac9e
# ╠═8896476e-7b72-11eb-362b-a9575e01485a
# ╟─d9f5d3e0-7b72-11eb-0054-eb985c9c1faa
# ╠═f21f22de-7b72-11eb-1519-d7080c3cf5f1
# ╟─2c6b86d8-7b73-11eb-3e8c-2daf70acbee4
# ╠═2c35a19e-7b73-11eb-0397-e5b698d004c7
# ╠═5e97b884-7b73-11eb-0f54-1d088b792897
# ╠═5e7f4b14-7b73-11eb-1a73-bd7c68ec5807
# ╠═5e6590fc-7b73-11eb-1815-636463bb5f9a
# ╠═5e4b0804-7b73-11eb-30dd-53f5fe5e80c0
# ╠═2bedb2bc-7b73-11eb-1f98-6d62cdcc853b
# ╠═5e2fdb88-7b73-11eb-3bb1-1ff54b515847
# ╟─835a03f0-7b73-11eb-1706-d98a60bba904
# ╠═a3d0e1a0-7b73-11eb-0a20-19d2921060c0
