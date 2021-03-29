### A Pluto.jl notebook ###
# v0.12.12

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 487cf5ba-5f1b-11eb-1ae4-f9853ae7161d
using ControlSystems

# ╔═╡ 0660cbfa-5f1d-11eb-25da-a3cb2a8a4c07
using PlutoUI

# ╔═╡ 4f496bf0-5f19-11eb-020c-d3139161e0e0
md"# Routh-Hurwitz Criterion: Example

Consider an open-loop system with the transfer function

$G_{OL}(s) = K ~ \frac{s^2 + s + 1}{s^4 + 3 s^3 + s^2 + s + 1}$

with $K>0$.
"

# ╔═╡ 47e59032-5f1b-11eb-1b42-e969d8c2528e
K = 1 # gain

# ╔═╡ a9d55c8a-5f1b-11eb-3da7-eb3d0ad7f636
G_ol = K * tf([1, 1, 1], [1, 3, 1, 1, 1])

# ╔═╡ 125fb99e-5f1c-11eb-311e-21ab8bdcc145
pole(G_ol) # Unstable poles

# ╔═╡ d29b20e8-5f1b-11eb-366e-7f377929f0fd
Tf = 20.0; # Final simulation time

# ╔═╡ a9bd0da6-5f1b-11eb-19b2-cd2036ffa8fe
stepplot(G_ol, Tf)

# ╔═╡ a9a0ccfe-5f1b-11eb-3874-7ff114111922
md"The closed-loop system is derived by

$G_{CL}(s) = \frac{G_{OL}(s)}{1 + G_{OL}(s)}$

and thus one yields

$G_{CL}(s) = \frac{K ~ \frac{s^2 + s + 1}{s^4 + 3 s^3 + s^2 + s + 1}}{1 + K ~ \frac{s^2 + s + 1}{s^4 + 3 s^3 + s^2 + s + 1}} 
= \frac{K \left( s^2 + s + 1 \right) }{s^4 + 3 s^3 + (K+1) s^2  + (K+1) s + 1 + K}.$
" 



# ╔═╡ 0a47ff48-5f1d-11eb-111b-594c8fbc6286
@bind K2 Slider(0:1:100, default=1)

# ╔═╡ 8d4536ce-5f1d-11eb-21b7-e18da00cea3a
md"Gain K = $K2"

# ╔═╡ a9850c30-5f1b-11eb-1ff3-b356ae3b6f84
G_cl = K2 * tf([1, 1, 1], [1, 3, (K2 +1), (K2 + 1), (K2 +1)])

# ╔═╡ a96a18b2-5f1b-11eb-0ad2-a99843cb0e1d
stepplot(G_cl, Tf, legend=false)

# ╔═╡ a933d40a-5f1b-11eb-3e43-f90737f5ce8c
md"## Characteristic polynom

The general characteristic polynom for a fourth-order system is given as 

$p(s) = a_{4} s^4 + a_{3} s^3 + a_{2} s^2  + a_{1} s + a_{0}$

and here it is noted as

$p(s) = s^4 + 3 s^3 + (K+1) s^2  + (K+1) s + 1 + K.$ 

## Routh-Hurwitz algorithm

**Step 1: Coefficients**

All coefficients have to be greater than zero: $a_{i} > 0$ for $i = 0, 1, \cdots, 4$. This is guaranteed if $K + 1 > 0$ or $K > -1$. 

**Step 2: 2-dimensional Matrix**

$H_{2} =
\begin{pmatrix}
a_{1} & a_{3} \\
a_{0} & a_{2}
\end{pmatrix}
=
\begin{pmatrix}
K+1 & 3 \\
K+1 & K+1
\end{pmatrix}$

The determinant of $H_{2}$ has to be greater than zero, as

$\det{H_{2}} = (K+1)^2 - 3~(K+1) = (K + 1) \left[K - 2 \right] > 0.$

This requirement is fulfilled if $K > 2$.

**Step 2: 3-dimensional Matrix**

$H_{3} =
\begin{pmatrix}
a_{1} & a_{3} & a_{5} \\
a_{0} & a_{2} & a_{4} \\
0 & a_{1} & a_{3}
\end{pmatrix}
=
\begin{pmatrix}
K+1 & 3    & 0 \\
K+1 & K+1  & 1 \\
0   & K+1  & 3 
\end{pmatrix}$

The determinant of $H_{3}$ has to be greater than zero, as

$\det{H_{3}} = 3 ~ (K+1)^2 -  (K+1)^2 - 9 (K +1)  = (K + 1) (2K - 7) > 0.$

This requirement is fulfilled if $K > 3.5$.

**Step 2: 4-dimensional Matrix**

Is not necessary here because $a_{4} = 1$ and thus the determinants of $H_{4}$ and  $H_{3}$ would be equal.

### Conclusion

The closed-loop system is stable for

$K \in (0, 3.5) \quad \text{or} \quad 0 < K < 3.5.$

## Critical gain

The gain $K = 3.5$ leads to a marginally stable system.
"

# ╔═╡ a7f30548-5f1b-11eb-2459-79b26b54fa2f
K3 = 3.5

# ╔═╡ 14ecd308-5f22-11eb-224c-29e06a289d3e
G_cl3 = K3 * tf([1, 1, 1], [1, 3, (K3 +1), (K3 + 1), (K3 +1)])

# ╔═╡ 19629418-5f22-11eb-1d72-fd2296daed64
stepplot(G_cl3, Tf, legend=false)

# ╔═╡ Cell order:
# ╟─4f496bf0-5f19-11eb-020c-d3139161e0e0
# ╠═487cf5ba-5f1b-11eb-1ae4-f9853ae7161d
# ╠═47e59032-5f1b-11eb-1b42-e969d8c2528e
# ╠═a9d55c8a-5f1b-11eb-3da7-eb3d0ad7f636
# ╠═125fb99e-5f1c-11eb-311e-21ab8bdcc145
# ╠═d29b20e8-5f1b-11eb-366e-7f377929f0fd
# ╠═a9bd0da6-5f1b-11eb-19b2-cd2036ffa8fe
# ╟─a9a0ccfe-5f1b-11eb-3874-7ff114111922
# ╠═0660cbfa-5f1d-11eb-25da-a3cb2a8a4c07
# ╠═0a47ff48-5f1d-11eb-111b-594c8fbc6286
# ╠═8d4536ce-5f1d-11eb-21b7-e18da00cea3a
# ╠═a9850c30-5f1b-11eb-1ff3-b356ae3b6f84
# ╠═a96a18b2-5f1b-11eb-0ad2-a99843cb0e1d
# ╟─a933d40a-5f1b-11eb-3e43-f90737f5ce8c
# ╠═a7f30548-5f1b-11eb-2459-79b26b54fa2f
# ╠═14ecd308-5f22-11eb-224c-29e06a289d3e
# ╠═19629418-5f22-11eb-1d72-fd2296daed64
