### A Pluto.jl notebook ###
# v0.12.12

using Markdown
using InteractiveUtils

# ╔═╡ ca748432-5f36-11eb-1368-edc4653d0abe
using ControlSystems

# ╔═╡ 81b85936-5f30-11eb-1a81-6f2131859292
md"# Controller design: Naslin approach

**Remark:** The Naslin approach is designed for systems with proportional and integral behaviour. 

Consider the transfer function 

$G_{CL}(s) = \frac{B(s)}{A(s)}$

of closed-loop system with polynomial

$A(s) = a_{0} + a_{1} s + a_{2} s^2 + \cdots a_{m} s^{m}.$

Firstly, the Naslin coefficients 

$\alpha_{i} = \frac{a_{i}^2}{a_{i-1} ~ a_{i+1}} \quad \text{with} \quad i = 1,2, \cdots, m-1$

have to be determined. This coefficients have to meet the requirement

$1.5 < \alpha_{i} < 2.5  \quad \text{with} \quad i = 1,2, \cdots, m-1.$

**Small coefficients** lead to a fast closed-loop system and

**Large coefficients** lead to a damped behaviour.


## Example

Consider a plant with transfer function

$G_{p}(s) = \frac{0.5}{(s+2)~(s+1)~(s+0.5)} = \frac{0.5}{s^3 + 3.5 s^2 + 3.5 s + 1}$

which shall be controlled by a PID controller

$G_{c}(s) = \frac{K_{P} s + K_{I} + K_{D} s^2}{s}.$

The closed-loop system is noted as

$G_{CL}(s) = \frac{G_{p}(s) G_{c}(s)}{1 + G_{p}(s) G_{c}(s)} = \frac{0.5 ~ (K_{P} s + K_{I} + K_{D} s^2)}{s^4 + 3.5 s^3 + 3.5 s^2 + s + 0.5 (K_{P} s + K_{I} + K_{D} s^2)}$

and further

$G_{CL}(s) = \frac{0.5 ~ (K_{P} s + K_{I} + K_{D} s^2)}{s^4 + 3.5 s^3 + (3.5 + 0.5 K_{D}) s^2 + (1 + 0.5 K_{P}) s + 0.5 K_{I}}$

### Naslin coefficients

The Naslin coefficients are calculated with

$\alpha_{1} = \frac{a_{1}^2}{a_{0} ~ a_{2}} = \frac{(1 + 0.5 K_{P})^2}{0.5 K_{I} ~ (3.5 + 0.5 K_{D})},$

$\alpha_{2} = \frac{a_{2}^2}{a_{1} ~ a_{3}} = \frac{(3.5 + 0.5 K_{D})^2}{(1 + 0.5 K_{P}) ~ 3.5} \quad \text{and}$

$\alpha_{3} = \frac{a_{3}^2}{a_{2} ~ a_{4}} = \frac{3.5^2}{3.5 + 0.5 K_{D}}.$

The PID coefficients shall be found such that $\alpha_{1} = \alpha_{2} = \alpha_{3} = 2$.

"

# ╔═╡ 74a6d8c8-614f-11eb-2733-2381d427ebbb
md"$2 \cdot (3.5 + 0.5~Kd) = 7 + Kd = 3.5^2$"

# ╔═╡ 149f7b58-5f36-11eb-1fe7-b12722d27c45
Kd = 3.5^2 - 7

# ╔═╡ 13de38c6-5f36-11eb-0a80-6b03868ac810
Kp = ((3.5 + 0.5*Kd)^2 - 7)/3.5

# ╔═╡ 9b4acda6-5f36-11eb-2449-29e43a81b1ca
Ki = (1 + 0.5*Kp)^2 / (3.5 + 0.5*Kd)

# ╔═╡ d5621b70-5f36-11eb-1005-d591f136aeef
Gcl = 0.5*tf([Kd, Kp, Ki], [1, 3.5, (3.5 + 0.5*Kd), (1 + 0.5*Kp), 0.5*Ki])

# ╔═╡ 53220e8a-5f37-11eb-28da-3df53dcb1f8f
Tf = 10.0; #Final simulation time

# ╔═╡ 5dd084d8-5f37-11eb-2af3-e5020b2b9b5b
stepplot(Gcl, Tf, legend=false)

# ╔═╡ Cell order:
# ╟─81b85936-5f30-11eb-1a81-6f2131859292
# ╠═74a6d8c8-614f-11eb-2733-2381d427ebbb
# ╠═149f7b58-5f36-11eb-1fe7-b12722d27c45
# ╠═13de38c6-5f36-11eb-0a80-6b03868ac810
# ╠═9b4acda6-5f36-11eb-2449-29e43a81b1ca
# ╠═ca748432-5f36-11eb-1368-edc4653d0abe
# ╠═d5621b70-5f36-11eb-1005-d591f136aeef
# ╠═53220e8a-5f37-11eb-28da-3df53dcb1f8f
# ╠═5dd084d8-5f37-11eb-2af3-e5020b2b9b5b
