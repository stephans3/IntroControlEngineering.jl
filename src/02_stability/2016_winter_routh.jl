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

# ╔═╡ 73ec54ce-5f26-11eb-3d8f-61ca5e27222b
using PlutoUI, ControlSystems

# ╔═╡ 6dcce080-5f25-11eb-1c7c-4b9e2a6a44c8
md"# Routh-Hurwitz Criterion: Example

Consider a plant with transfer function

$G_{p}(s) = \frac{1}{s(s + 1)(s + 2)}.$

## First Controller approach

A proportional controller $G_{c}(s) = K$ is assumed and the open-loop system is noted as

$G_{OL}(s) = G_{c}(s) ~ G_{p}(s) = \frac{K}{s(s + 1)(s + 2)}$

and the closed-loop system is derived as

$G_{CL}(s) = \frac{G_{OL}(s)}{1 + G_{OL}(s)} = \frac{K}{s^3 + 3 s^2 + 2 s + K}.$

"

# ╔═╡ 2e67152a-5f2a-11eb-1328-e3f6a7d58814
@bind K Slider(0:0.1:8, default=1)

# ╔═╡ 82fc23d6-5f26-11eb-3716-2b0b3ab4fcff
md"Gain K = $K"

# ╔═╡ 89001864-5f26-11eb-3c12-29ec3a44b3e6
G_cl = K * tf([1], [1, 3, 2, K])

# ╔═╡ bee58aae-5f26-11eb-0cb5-ab467800faf4
Tf = 20.0; # Final simulation time

# ╔═╡ cb98199c-5f26-11eb-33bb-6b209e361e16
stepplot(G_cl, Tf, legend=false)

# ╔═╡ f0dc344a-5f26-11eb-03f3-915f38511246
md"## Characteristic polynom

The general characteristic polynom for a third-order system is given as 

$p(s) = a_{3} s^3 + a_{2} s^2  + a_{1} s + a_{0}$

and here it is noted as

$p(s) = s^3 + 3 s^2 + 2 s + K.$ 

## Routh-Hurwitz algorithm

**Step 1: Coefficients**

All coefficients have to be greater than zero: $a_{i} > 0$ for $i = 0, 1, \cdots, 3$. This is guaranteed if $K + 1 > 0$ or $K > 0$. 

**Step 2: 2-dimensional Matrix**

$H_{2} =
\begin{pmatrix}
a_{1} & a_{3} \\
a_{0} & a_{2}
\end{pmatrix}
=
\begin{pmatrix}
2 & 1 \\
K & 3
\end{pmatrix}$

The determinant of $H_{2}$ has to be greater than zero, as

$\det{H_{2}} = 6 - K > 0.$

This requirement is fulfilled if $K < 6$.

**Step 2: 3-dimensional Matrix**

$H_{3} =
\begin{pmatrix}
a_{1} & a_{3} & a_{5} \\
a_{0} & a_{2} & a_{4} \\
0 & a_{1} & a_{3}
\end{pmatrix}
=
\begin{pmatrix}
2 & 1 & 0 \\
K & 3 & 0 \\
0 & 2 & 1
\end{pmatrix}$

The determinant of $H_{3}$ has to be greater than zero, as

$\det{H_{3}} = 6 - K > 0.$

This requirement is fulfilled if $K < 6$.

### Conclusion

The closed-loop system is stable for 

$K \in (0,6) \quad \text{or} \quad 0 < K < 6.$


## Second Controller approach

A new [lead controller](https://en.wikipedia.org/wiki/Lead%E2%80%93lag_compensator) controller 

$G_{c} = K ~ \frac{s + 2}{s + 2.5}$

is assumed and the open-loop system is noted as

$G_{OL}(s) = G_{c} ~ G_{p} = \frac{K}{s(s + 1)(s + 2.5)}$

and the closed-loop system is derived as

$G_{CL}(s) = \frac{G_{OL}(s)}{1 + G_{OL}(s)} = \frac{K}{s^3 + 3.5 s^2 + 2.5 s + K}.$


"

# ╔═╡ 3674d14e-5f2a-11eb-0f5f-e1a4b856bda4
@bind K2 Slider(0:0.1:12, default=1)

# ╔═╡ 36512406-5f2a-11eb-20cd-11cf52e75fdb
md"Gain K = $K2"

# ╔═╡ 3630373c-5f2a-11eb-3f0b-01316e62fa8f
G_cl2 = K2 * tf([1], [1, 3.5, 2.5, K2])

# ╔═╡ 3618ad74-5f2a-11eb-3e44-65068875b562
stepplot(G_cl2, Tf, legend=false)

# ╔═╡ 35db5ae4-5f2a-11eb-2b07-93efcfb2d3d4
md"## Characteristic polynom

The general characteristic polynom for a third-order system is given as 

$p(s) = a_{3} s^3 + a_{2} s^2  + a_{1} s + a_{0}$

and here it is noted as

$p(s) = s^3 + 3.5 s^2 + 2.5 s + K.$ 

## Routh-Hurwitz algorithm

**Step 1: Coefficients**

All coefficients have to be greater than zero: $a_{i} > 0$ for $i = 0, 1, \cdots, 3$. This is guaranteed if $K + 1 > 0$ or $K > 0$. 

**Step 2: 2-dimensional Matrix**

$H_{2} =
\begin{pmatrix}
a_{1} & a_{3} \\
a_{0} & a_{2}
\end{pmatrix}
=
\begin{pmatrix}
2.5 & 1 \\
K & 3.5
\end{pmatrix}$

The determinant of $H_{2}$ has to be greater than zero, as

$\det{H_{2}} = 8.75 - K > 0.$

This requirement is fulfilled if $K < 8.75$.

**Step 3: 3-dimensional Matrix**

Is not necessary here because $a_{3} = 1$ and thus the determinants of $H_{2}$ and  $H_{3}$ would be equal.

### Conclusion

The closed-loop system is stable for 

$K \in (0,8.75) \quad \text{or} \quad 0 < K < 8.75.$

"




# ╔═╡ 35b2912c-5f2a-11eb-224c-15559cc36eac


# ╔═╡ Cell order:
# ╟─6dcce080-5f25-11eb-1c7c-4b9e2a6a44c8
# ╠═73ec54ce-5f26-11eb-3d8f-61ca5e27222b
# ╠═2e67152a-5f2a-11eb-1328-e3f6a7d58814
# ╠═82fc23d6-5f26-11eb-3716-2b0b3ab4fcff
# ╠═89001864-5f26-11eb-3c12-29ec3a44b3e6
# ╠═bee58aae-5f26-11eb-0cb5-ab467800faf4
# ╠═cb98199c-5f26-11eb-33bb-6b209e361e16
# ╟─f0dc344a-5f26-11eb-03f3-915f38511246
# ╠═3674d14e-5f2a-11eb-0f5f-e1a4b856bda4
# ╠═36512406-5f2a-11eb-20cd-11cf52e75fdb
# ╠═3630373c-5f2a-11eb-3f0b-01316e62fa8f
# ╠═3618ad74-5f2a-11eb-3e44-65068875b562
# ╟─35db5ae4-5f2a-11eb-2b07-93efcfb2d3d4
# ╠═35b2912c-5f2a-11eb-224c-15559cc36eac
