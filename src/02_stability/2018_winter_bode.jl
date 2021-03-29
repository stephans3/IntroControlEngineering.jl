### A Pluto.jl notebook ###
# v0.12.12

using Markdown
using InteractiveUtils

# ╔═╡ 3c53f232-559f-11eb-248a-1be39821e914
using ControlSystems

# ╔═╡ 3273ac9c-559c-11eb-1433-191625b806eb
md"# Bode plot of a delay system

Consider the transfer function of a plant

$G_{p}(s) = \frac{0.5}{s + 0.5}.$

A proportional-integral (PI) controller is assumed as 

$G_{c}(s) = K ~ \frac{s + 0.5}{s} = K ~ \left( 1 + \frac{1}{2~s}\right).$

Thus, the open-loop system is found by 

$G_{ol}(s) = G_{p}(s) \cdot G_{c}(s) = K ~ \frac{1}{2~s}.$

The closed-loop system is calculated with 

$G_{cl}(s) = \frac{G_{ol}(s)}{1 + G_{ol}(s)} = \frac{K}{2s + K}$

and its pole is determinded with 

$2~s + K = 0 \quad \text{or equivalent} \quad  s = -\frac{K}{2}.$

The closed-loop system is bounded-input-bounded-output (BIBO) stable if $K>0$ which implies that $s<0$. For simplicity, the gain is assumed as $K=1$.

**This result shall be proved using the Bode plot!**
"

# ╔═╡ 3c36448a-559f-11eb-23f5-13523a76140b
K = 1 # Gain

# ╔═╡ 3c184f86-559f-11eb-23bc-01eae2cf6854
Gol = tf(K, [2, 0]) # Open-loop system

# ╔═╡ 3b8c7eaa-559f-11eb-3a67-37c76ce8e1ba
ω = 10^(-2) : 10^(-2)  : 10^2 # Frequencies

# ╔═╡ 17388f3e-55a0-11eb-2725-b7d259f6cce2
md"## Bode plot"

# ╔═╡ 3b72683a-559f-11eb-2888-716bc1a61a69
bodeplot(Gol, ω, legend=false)

# ╔═╡ 3649a270-55a0-11eb-1f73-f1083089be84
setPlotScale("dB")

# ╔═╡ 4f191b6c-55a0-11eb-1488-e316052b4a16
md"**The phase does not reach $-180°$ and so it is stable for all positive values of $K$.**"

# ╔═╡ fee8d440-55a7-11eb-2c99-a318a4f62b77
md"### Step plot of closed-loop system"

# ╔═╡ ce10a302-55a7-11eb-1103-2b0d203ae308
Gcl = feedback(Gol) # Closed-loop system

# ╔═╡ 6d5d6c90-55a5-11eb-21f1-ebf693fb4034
Tf = 20.0 # Final simulation time

# ╔═╡ e9420e5e-55a7-11eb-2781-612fb9a094c9
stepplot(Gcl, Tf)

# ╔═╡ 9657d2de-55a0-11eb-18fa-9bd5e1c774e3
md"## Additional delay

From here on an additional delay $T_{t}=0.7$ seconds is assumed for the plant. The new transfer function of the plant is given by

$G_{p}(s) = \frac{0.5}{s + 0.5} ~ e^{-0.7 ~s}$

and the open-loop system is formulated by

$G_{ol}(s) = G_{p}(s) \cdot G_{c}(s) = K ~ \frac{1}{2~s} ~ e^{-0.7 ~s}.$

Again, gain $K=1$ is assumed and the stability is checked with the Bode plot.
"

# ╔═╡ c8bce4c4-55a2-11eb-06a1-6185c826d1d7
Tt = 0.7 # Delay

# ╔═╡ 963a829c-55a0-11eb-136d-2f39bbd91bfb
Gol_delay = tf(K, [2.0, 0.0]) * delay(Tt) # Open-loop with delay

# ╔═╡ f6457082-55a2-11eb-036a-c71c681dc8fd
md"## Bode plot"

# ╔═╡ 2d91fa60-55a3-11eb-0796-6b069af25e37
ω_delay = 10^(-2) : 10^(-2) : 0.5*10^(1)

# ╔═╡ 15eb3854-55a3-11eb-2bd8-5b41b6657826
bodeplot(Gol_delay, ω_delay, legend=false)

# ╔═╡ aace512c-55a3-11eb-35e0-ef9c2b1d31c5
ϕ = bode(Gol_delay, [0.5])[2][1, 1, 1] # Phase at ω=0.5

# ╔═╡ 55724a70-55a4-11eb-0e03-d9b87d1c5d35
md"## Stability

At frequency $\omega = 0.5$ the magnitude is almost zero dB. The phase at $\omega = 0.5$ is 

$\phi(0.5) \approx -110°$

The closed-loop system with delay is BIBO stable because its phase margin is greater than zero as

$\Delta \phi = \phi(0.5) - (-180°) \approx 70°.$"


# ╔═╡ 935f1862-55aa-11eb-3e71-bdcfbb109542
md"### Step plot of closed-loop system"

# ╔═╡ 9001252a-55a5-11eb-030e-c57a5458fb59
Gcl_delay = feedback(Gol_delay)

# ╔═╡ c829fc38-55a5-11eb-1431-9310dcc80b15
stepplot(Gcl_delay,Tf)

# ╔═╡ f6f764b2-55a5-11eb-3b9b-0b6361bc29f4
md"## Adjustment

If one wants to decrease the phase margin to $30°$ which corresponds with phase $\phi = -150°$ then the magnitude at for this frequency has to be found. Here, the phase $\phi = -150°$ corresponds almost with the magnitude $A[dB] = -10$ and the gain margin is $\Delta A[dB] = 0 - A[dB] = 10$. 

In absolute values the magnitude is calculated with

$V = 10^\frac{\Delta A[dB]}{20} = 10^{\frac{10}{20}} = \sqrt{10} \approx 3.16.$

The new system is noted as 

$G_{2}(s) = V~ K ~ \frac{1}{2~s} ~ e^{-0.7 ~s} = \sqrt{10} ~ \frac{1}{2~s} ~ e^{-0.7 ~s}.$

"

# ╔═╡ de2b33c8-55a3-11eb-0632-ffd841c9754f
bode(Gol_delay, [1.5])

# ╔═╡ c364bafa-55a8-11eb-3cdf-977eb42abc24
V = sqrt(10)

# ╔═╡ bb8fdd46-55a8-11eb-294d-630c6f4c89a2
Gol_delay2 = V * Gol_delay;

# ╔═╡ e64c49e6-55aa-11eb-182d-3bd624874c85
Gcl_delay2 = feedback(Gol_delay2)

# ╔═╡ f279d012-55aa-11eb-1ed5-8da58f1347cf
stepplot(Gcl_delay2,Tf)

# ╔═╡ 68d9e8c8-55ab-11eb-3e60-7b43d176c059
md"### Comparison of Nyquist plots"

# ╔═╡ 1ed7f8d6-55a7-11eb-01a2-4722ebb829b1
nyquistplot(Gol_delay, 0.4 : 0.1 : 100.0, legend=false)

# ╔═╡ a41984c8-55a8-11eb-32b4-3bbe60da8528
nyquistplot(Gol_delay2, legend=false)

# ╔═╡ Cell order:
# ╟─3273ac9c-559c-11eb-1433-191625b806eb
# ╠═3c53f232-559f-11eb-248a-1be39821e914
# ╠═3c36448a-559f-11eb-23f5-13523a76140b
# ╠═3c184f86-559f-11eb-23bc-01eae2cf6854
# ╠═3b8c7eaa-559f-11eb-3a67-37c76ce8e1ba
# ╟─17388f3e-55a0-11eb-2725-b7d259f6cce2
# ╠═3b72683a-559f-11eb-2888-716bc1a61a69
# ╟─3649a270-55a0-11eb-1f73-f1083089be84
# ╟─4f191b6c-55a0-11eb-1488-e316052b4a16
# ╟─fee8d440-55a7-11eb-2c99-a318a4f62b77
# ╠═ce10a302-55a7-11eb-1103-2b0d203ae308
# ╠═6d5d6c90-55a5-11eb-21f1-ebf693fb4034
# ╠═e9420e5e-55a7-11eb-2781-612fb9a094c9
# ╟─9657d2de-55a0-11eb-18fa-9bd5e1c774e3
# ╠═c8bce4c4-55a2-11eb-06a1-6185c826d1d7
# ╠═963a829c-55a0-11eb-136d-2f39bbd91bfb
# ╟─f6457082-55a2-11eb-036a-c71c681dc8fd
# ╠═2d91fa60-55a3-11eb-0796-6b069af25e37
# ╠═15eb3854-55a3-11eb-2bd8-5b41b6657826
# ╠═aace512c-55a3-11eb-35e0-ef9c2b1d31c5
# ╟─55724a70-55a4-11eb-0e03-d9b87d1c5d35
# ╟─935f1862-55aa-11eb-3e71-bdcfbb109542
# ╠═9001252a-55a5-11eb-030e-c57a5458fb59
# ╠═c829fc38-55a5-11eb-1431-9310dcc80b15
# ╟─f6f764b2-55a5-11eb-3b9b-0b6361bc29f4
# ╠═de2b33c8-55a3-11eb-0632-ffd841c9754f
# ╠═c364bafa-55a8-11eb-3cdf-977eb42abc24
# ╠═bb8fdd46-55a8-11eb-294d-630c6f4c89a2
# ╠═e64c49e6-55aa-11eb-182d-3bd624874c85
# ╠═f279d012-55aa-11eb-1ed5-8da58f1347cf
# ╠═68d9e8c8-55ab-11eb-3e60-7b43d176c059
# ╠═1ed7f8d6-55a7-11eb-01a2-4722ebb829b1
# ╠═a41984c8-55a8-11eb-32b4-3bbe60da8528
