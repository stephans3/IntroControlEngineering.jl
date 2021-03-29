### A Pluto.jl notebook ###
# v0.12.12

using Markdown
using InteractiveUtils

# ╔═╡ 1db03632-802f-11eb-3ecf-7f4afb12854d
using LinearAlgebra

# ╔═╡ 4518d7f2-80c4-11eb-1081-c785fed963bf
using DifferentialEquations

# ╔═╡ 77b9bb72-80c4-11eb-10a4-7b7a254c7469
using Plots

# ╔═╡ e57510ba-802d-11eb-3ebb-db07872a1197
md"# Electrical Oscillator II: State Observer

Consider again the electrial oscillator in state-space representation

$\dot{x}(t) = A x(t) + B u(t)$

$y(t) = C x(t)$

with $x \in \mathbb{R}^3$ and

$A = 
\begin{pmatrix}
0 & 1 & 0 \\
0 & 0 & 1 \\
-d_{0} & -d_{1} & -d_{2}
\end{pmatrix}
\quad \text{,} \quad
B =
\begin{pmatrix}
0 \\
0 \\
1
\end{pmatrix}
\quad \text{,} \quad
C = 
\begin{pmatrix}
d_{0} & 0 & 0
\end{pmatrix} \text{.}$


"

# ╔═╡ f11fcb00-802e-11eb-3cf5-89740c671213
R = 10^4 # Resistivity: 10kΩ 

# ╔═╡ 03009a2a-802f-11eb-1201-7dad5001b9a6
Cap = 100 * 10^(-6) # Capacity 100μF

# ╔═╡ 0908a1a6-802f-11eb-17e0-a5a373cfb741
T₁ = T₂ = T₃ = R*Cap # Time constants

# ╔═╡ 0d348baa-802f-11eb-1f6a-6bd3f3cef719
d₀ = 1/(T₁ * T₂ * T₃)

# ╔═╡ 10bc058c-802f-11eb-13c3-39e4c666b0d9
d₁ = (T₁ + T₂ + T₃)/(T₁ * T₂ * T₃)

# ╔═╡ 14757b5c-802f-11eb-0a11-61792cd1ee69
d₂ = (T₁*T₂ + T₁*T₃ + T₂*T₃)/(T₁ * T₂ * T₃)

# ╔═╡ 19c30432-802f-11eb-34e7-7133b9a99547
A = [0 1 0; 0 0 1; -d₀ -d₁ -d₂]

# ╔═╡ 1f6b0dee-802f-11eb-1e5b-f13a3a2f91cb
md"## Observability and State Observer

If the initial state of a system can be computed from the knowledge of all input $u(t)$ and output signals $y(t)$ from $t=0$ until some (final) time $t=T_{f}$ then the system is called observable. This means, the real state can be calculated from the knowledge of the system input and output. The observability property can be checked with various methods, e.g. with the Kalman rank condidtion. 

Firstly, the observability matrix is built with

$\mathcal{O} = 
\begin{pmatrix}
C \\
C ~ A \\
C ~ A^2 
\end{pmatrix}
=
\begin{pmatrix}
d_{0} & 0 & 0 \\
0 & d_{0} & 0 \\
0 & 0 &  d_{0}
\end{pmatrix}$

and secondly it is checked whether $\mathcal{O}$ has full rank with

$\det(\mathcal{O}) = d_{0}^3 \neq 0 \quad \text{if} \quad d_{0} \neq 0.$


"

# ╔═╡ e00a98aa-8032-11eb-3b06-c3e114f45586
C = [d₀ 0 0]

# ╔═╡ 06c7bd3c-8034-11eb-320a-93da3fc398f2
C*A

# ╔═╡ 29e25480-8034-11eb-381f-07868e859375
Obsvmat = vcat(C, C*A, C*A^2)

# ╔═╡ 44f520ce-8034-11eb-144a-579ffa4a0067
rank(Obsvmat)

# ╔═╡ 775a4f06-80bb-11eb-1c31-b1b2efb0a364
md"The [state observer](https://en.wikipedia.org/wiki/State_observer) is a model of the original dynamical system which is updated by input $u(t)$ and real output $y(t)$. The difference between the real $y(t)$ and the observer output $\hat{y}(t)$ is noted as the observer error 

$e_{obs}(t) = y(t) - \hat{y}(t).$ 

The observer dynamics is stated as

$\dot{\hat{x}}(t) = A ~ \hat{x} + B ~ u(t) + L ~ e_{obs}(t)$

$\hat{y}(t) = C ~ \hat{x}(t)$

with observer gain vector (or matrix) $L$. The observer dynamics is further formed as 

$\dot{\hat{x}}(t) = A ~ \hat{x} + B ~ u(t) + L ~ (y(t) - \hat{y}(t)) = (A - L~C)~ \hat{x} + B u(t) + L y(t)$

and the observer gain $L$ has to be found such that the observer error dynamics 

$\frac{d}{dt} {e}_{obs}(t) = (A - L~C)~ e(t)$

is asymptotical stable. That means $A_{obs} = (A - L~C)$ has stable eigenvalues and so it is guaranteed that 

$e_{obs}(t) \rightarrow 0 \quad \text{for} \quad t \rightarrow 0.$

Here, the observer state matrix yields the form 

$A_{obs} = (A - L~C) =
\begin{pmatrix}
-l_{1} d_{0}  & 1 & 0 \\
-l_{2} d_{0} & 0 & 1 \\
-(l_{3} + 1) ~d_{0} & -d_{1} & -d_{2}
\end{pmatrix}$

where $L = \begin{pmatrix} l_{0}, l_{1}, l_{2} \end{pmatrix}^{\top}$. The observer has to be faster than the original dynamical system and therefore the eigenvalues of $A_{obs}$ have to be smaller than the eigenvalues of $A$.

The observer gain vector is calculated using [Ackermann's formula](https://en.wikipedia.org/wiki/Ackermann%27s_formula) as

$L^{\top} = 
\begin{pmatrix}
0 & 0 & \cdots 0 & 1
\end{pmatrix}
~ \mathcal{O}^{-1} ~ p(A^{\top})$

in which $p(\lambda)$ is the characteristic polynomial of the desired eigenvalues.

Here, the following Eigenvalues $\Lambda = \left\{-8, -8, -8\right\}$ are desired to guarantee a fast observer a and thus one holds the characteristic polynomial

$p(\lambda) = (8+\lambda)^3 = \lambda^3 + 24 \lambda^2 + 192 \lambda + 512.$


"

# ╔═╡ bff2905e-80dc-11eb-20a1-ed018c449524
λdes = (-8, -8, -8) # Desired observer Eigenvalues

# ╔═╡ 6970983c-80dc-11eb-244a-634406f79797
p(x, λ) = (x - λ[1]*I) * (x - λ[2]*I) * (x - λ[3]*I)   

# ╔═╡ b67cb5f4-80dc-11eb-144d-d91d56d2a3df


# ╔═╡ 6012195a-80c0-11eb-21fa-e908b8753736
L1 = [0 0 1] * inv(Obsvmat) * p(A', λdes)

# ╔═╡ 6604e282-80db-11eb-1bbf-7b22dabd5e76
md"**Observer gain vector**"

# ╔═╡ 11af2eca-80db-11eb-0044-0d3900a52be2
L = transpose(L1)

# ╔═╡ 8bc82be6-80db-11eb-04df-7589e2a85705


# ╔═╡ 7274019c-80db-11eb-3ce3-ff25243ac526
md"**Observer state matrix**"

# ╔═╡ 52aaffec-80c1-11eb-3bb3-afbdb286dd9f
Aobs = A - L*C

# ╔═╡ 6e69ae68-80c1-11eb-0fae-fbf235d6094a
evals_obs = eigvals(Aobs)

# ╔═╡ 9e3d9a18-80db-11eb-0806-b5fbb7531259
md"### Simulation of the observer"

# ╔═╡ ca80630e-80c1-11eb-269e-4319ef335279
function rc_observed(dx, x, p, t) # without input
	
	x_real = @view x[1:3]	# real states
	dx_real = @view dx[1:3]
	
	x_obs = @view x[4:6]	# observed states
	dx_obs = @view dx[4:6]
	
	y_real = C * x_real
	
	dx_real .= A * x_real 	# real state dynamics
	dx_obs .= Aobs * x_obs + L * y_real # observer dynamics

end

# ╔═╡ 0a6467b8-80c4-11eb-1475-0f8e1c7662e7
x_init = 10*rand(6) # initial real + observer states 

# ╔═╡ 1939a5ec-80c4-11eb-0048-25224bc9481c
tspan = (0.0, 10.0)

# ╔═╡ 2df13272-80c4-11eb-379d-85afb5f46a6f
prob_obs = ODEProblem(rc_observed, x_init, tspan)  # Dynamical system (real + observer) with random initial conditions

# ╔═╡ 4e3943a8-80c4-11eb-2f01-fb314c8c336d
sol_obs = solve(prob_obs)

# ╔═╡ 779ed7e4-80c4-11eb-3fda-6f55cd1e1f14
plot(sol_obs, label=["x1" "x2" "x3" "x̂1" "x̂2" "x̂3"], title="Open-loop System with Observer")

# ╔═╡ 4b504692-80c5-11eb-217f-3f06c95da132
sol_obs[1:3,:]

# ╔═╡ 7662cba6-80c4-11eb-3cbd-f328eedf0bc2
err_obs = sol_obs[1:3,:] - sol_obs[4:6,:]

# ╔═╡ 7e7dbdb8-80c5-11eb-0b9f-a3743f09e304
plot(sol_obs.t, err_obs', title="Observation Error")

# ╔═╡ 182b4ebc-80c6-11eb-395a-318c8191e822
md"## State Observer and Controller

The new closed-loop system consists of the original dynamical system, a state observer and a state feedback controller. The state controller is now defined as

$u(t) = -K ~ \hat{x}(t)$

because the internal state can not be measured directly. Here, the feedback law from the previous part is used.  

"

# ╔═╡ 2f86b7ae-80c6-11eb-134d-6384fd3223fd
B = [0; 0; 1]

# ╔═╡ 66df100c-80c6-11eb-2c34-47962f3c44c6
Cntrmat = hcat(B, A*B, A^2*B)

# ╔═╡ 56f602e8-80e1-11eb-23c4-efbe0fbedcd5
λcntr = (-4, -2, -2) # Eigenvalues of the controller

# ╔═╡ 4f6a215a-80c6-11eb-22bc-8f699c240f61
K = [0 0 1] * inv(Cntrmat) * p(A, λcntr)

# ╔═╡ e2fd4fc8-80df-11eb-30b0-137f4d69ad9b
md"The observer dynamics is noted in terms of states $(x, \hat{x})$ as

$\dot{\hat{x}}(t) =  (A_{obs} - B~K)~ \hat{x} + L ~ C ~ x(t)$

and the complete system is noted as

$\begin{pmatrix}
\dot{x}(t) \\
\dot{\hat{x}}(t)
\end{pmatrix}
=
\begin{pmatrix}
A & -B K \\
L~C & A_{obs} - B~K
\end{pmatrix}
\begin{pmatrix}
{x}(t) \\
{\hat{x}}(t)
\end{pmatrix} \tag{1}$

with output

$y(t) = C ~ x(t).$"

# ╔═╡ 2f6cce84-80c6-11eb-28e4-2125fb4cecdb
function rc_cntr_obs(dx, x, p, t) # Controller and Observer
	
	x_real = @view x[1:3]	# real states
	dx_real = @view dx[1:3]
	
	x_obs = @view x[4:6]	# observed states
	dx_obs = @view dx[4:6]
	
	u = -K * x_obs # State feedback
	u = u[1]
	y_real = C * x_real
	
	dx_real .= A * x_real + B*u 	# real state dynamics
	dx_obs .= Aobs * x_obs + B*u  + L * y_real # observer dynamics

end

# ╔═╡ bdc51770-80c6-11eb-1a4a-6fe5b5ac5fc6
prob_co = ODEProblem(rc_cntr_obs, x_init, tspan)  # Dynamical system (real + observer) with random initial conditions

# ╔═╡ 438e47d4-80c7-11eb-369e-713c9ca592cd
sol_co = solve(prob_co);

# ╔═╡ ccaa2bf6-80c6-11eb-1399-0fad937143e8
plot(sol_co, label=["x1" "x2" "x3" "x̂1" "x̂2" "x̂3"], title="Closed-loop system")

# ╔═╡ 813bfd14-80d2-11eb-1d84-f5f4f739dfef
err_co = sol_co[1:3,:] - sol_co[4:6,:]

# ╔═╡ 8bbdbea8-80d2-11eb-0137-8100ffe88034
plot(sol_obs.t, err_obs', title="Observation Error")

# ╔═╡ 4e99453c-80e2-11eb-10f0-d5db9ac18ada
md"## Digital Control Algorithm

The state matrix 

$A_{co} = 
\begin{pmatrix}
A & -B K \\
L~C & A_{obs} - B~K
\end{pmatrix}$ 

of the controller-observer system $(1)$ is transfered to a discrete state matrix

$A_{d,co} = e^{A_ {co} \Delta T}$

with sampling time $\Delta T=0.01$ seconds. The resulting algorithm is

$x_{co}(t_{k+1}) = A_{d,co} ~ x_{co}(t_{k})$

with united state vector
$x_{co} = 
\begin{pmatrix}
x(t) \\
\hat{x}(t)
\end{pmatrix}$.

Additionally, a noise signal is added to the real states as

$\begin{pmatrix}
x(t_{k+1}) \\
\hat{x}(t_{k+1})
\end{pmatrix}
= A_{d,co} ~ 
\begin{pmatrix}
x(t_{k+1}) \\
\hat{x}(t_{k+1})
\end{pmatrix}
+
\begin{pmatrix}
w(t_{k}) \\
0
\end{pmatrix}.$


"

# ╔═╡ 4e5298d0-80e2-11eb-32e4-f9c7abcfb893
Aco = zeros(6,6); # Building state matrix of controller and observer part

# ╔═╡ 4e336ad4-80e2-11eb-0321-29c9f75f653b
Aco[1:3,1:3] = A

# ╔═╡ 76ba4068-80f4-11eb-18b0-d79b01c0f478
Aco[1:3,4:6] = -B*K

# ╔═╡ b281dcf0-80f4-11eb-39e3-752784ef95ed
Aco[4:6,1:3] = L*C

# ╔═╡ c3867380-80f4-11eb-1a1e-61892538c262
Aco[4:6,4:6] = Aobs - B*K

# ╔═╡ d78fbce2-80f4-11eb-3a15-09eb96356fec
eigvals(Aco) # Eigenvalues of the state matrix

# ╔═╡ 0f647d7e-80f5-11eb-13a0-9b1608eb924a
ΔT = 0.01 # Sampling time

# ╔═╡ 07d2deac-80f5-11eb-21a9-f3d7da20e7b6
Adco = exp(Aco*ΔT) # Discrete state matrix

# ╔═╡ 386472f6-80f5-11eb-01e2-1b894c25ce58
xd0 = 10*rand(6) # Initial values

# ╔═╡ 55d89790-80f5-11eb-2dac-b1a2a371b578
tsteps = 1000; # time steps

# ╔═╡ 4b82bb72-80f5-11eb-0c31-c37caf088c8e
xd = zeros(6, tsteps);

# ╔═╡ 689cca2c-80f5-11eb-0f58-6bc0f7ab6acb
xd[:,1] = xd0

# ╔═╡ 79179e2c-80f5-11eb-0254-c7d9086d7a31
for i = 1 : tsteps-1
	w = 0.01*rand(3) # noise
	xd[1:3,i] = xd[1:3,i] + w # additive noise
	xd[:,i+1] = Adco * xd[:,i]
end

# ╔═╡ 9e73856e-80f5-11eb-2a09-eb3499ab0bc5
plot(xd', label=["x1" "x2" "x3" "x̂1" "x̂2" "x̂3"], title="Digital closed-loop system with noise")

# ╔═╡ f2c0a914-80f5-11eb-144c-99c87d037490
error_co = xd[1:3,:] - xd[4:6,:] 

# ╔═╡ 12a86daa-80f6-11eb-3ff4-a7d9bbf0f8a6
plot(error_co', title="Observation error")

# ╔═╡ Cell order:
# ╟─e57510ba-802d-11eb-3ebb-db07872a1197
# ╠═f11fcb00-802e-11eb-3cf5-89740c671213
# ╠═03009a2a-802f-11eb-1201-7dad5001b9a6
# ╠═0908a1a6-802f-11eb-17e0-a5a373cfb741
# ╠═0d348baa-802f-11eb-1f6a-6bd3f3cef719
# ╠═10bc058c-802f-11eb-13c3-39e4c666b0d9
# ╠═14757b5c-802f-11eb-0a11-61792cd1ee69
# ╠═19c30432-802f-11eb-34e7-7133b9a99547
# ╠═1db03632-802f-11eb-3ecf-7f4afb12854d
# ╟─1f6b0dee-802f-11eb-1e5b-f13a3a2f91cb
# ╠═e00a98aa-8032-11eb-3b06-c3e114f45586
# ╠═06c7bd3c-8034-11eb-320a-93da3fc398f2
# ╠═29e25480-8034-11eb-381f-07868e859375
# ╠═44f520ce-8034-11eb-144a-579ffa4a0067
# ╟─775a4f06-80bb-11eb-1c31-b1b2efb0a364
# ╠═bff2905e-80dc-11eb-20a1-ed018c449524
# ╠═6970983c-80dc-11eb-244a-634406f79797
# ╠═b67cb5f4-80dc-11eb-144d-d91d56d2a3df
# ╠═6012195a-80c0-11eb-21fa-e908b8753736
# ╟─6604e282-80db-11eb-1bbf-7b22dabd5e76
# ╠═11af2eca-80db-11eb-0044-0d3900a52be2
# ╠═8bc82be6-80db-11eb-04df-7589e2a85705
# ╟─7274019c-80db-11eb-3ce3-ff25243ac526
# ╠═52aaffec-80c1-11eb-3bb3-afbdb286dd9f
# ╠═6e69ae68-80c1-11eb-0fae-fbf235d6094a
# ╟─9e3d9a18-80db-11eb-0806-b5fbb7531259
# ╠═ca80630e-80c1-11eb-269e-4319ef335279
# ╠═0a6467b8-80c4-11eb-1475-0f8e1c7662e7
# ╠═1939a5ec-80c4-11eb-0048-25224bc9481c
# ╠═4518d7f2-80c4-11eb-1081-c785fed963bf
# ╠═2df13272-80c4-11eb-379d-85afb5f46a6f
# ╠═4e3943a8-80c4-11eb-2f01-fb314c8c336d
# ╠═77b9bb72-80c4-11eb-10a4-7b7a254c7469
# ╠═779ed7e4-80c4-11eb-3fda-6f55cd1e1f14
# ╠═4b504692-80c5-11eb-217f-3f06c95da132
# ╠═7662cba6-80c4-11eb-3cbd-f328eedf0bc2
# ╠═7e7dbdb8-80c5-11eb-0b9f-a3743f09e304
# ╟─182b4ebc-80c6-11eb-395a-318c8191e822
# ╠═2f86b7ae-80c6-11eb-134d-6384fd3223fd
# ╠═66df100c-80c6-11eb-2c34-47962f3c44c6
# ╠═56f602e8-80e1-11eb-23c4-efbe0fbedcd5
# ╠═4f6a215a-80c6-11eb-22bc-8f699c240f61
# ╠═e2fd4fc8-80df-11eb-30b0-137f4d69ad9b
# ╠═2f6cce84-80c6-11eb-28e4-2125fb4cecdb
# ╠═bdc51770-80c6-11eb-1a4a-6fe5b5ac5fc6
# ╠═438e47d4-80c7-11eb-369e-713c9ca592cd
# ╠═ccaa2bf6-80c6-11eb-1399-0fad937143e8
# ╠═813bfd14-80d2-11eb-1d84-f5f4f739dfef
# ╠═8bbdbea8-80d2-11eb-0137-8100ffe88034
# ╟─4e99453c-80e2-11eb-10f0-d5db9ac18ada
# ╠═4e5298d0-80e2-11eb-32e4-f9c7abcfb893
# ╠═4e336ad4-80e2-11eb-0321-29c9f75f653b
# ╠═76ba4068-80f4-11eb-18b0-d79b01c0f478
# ╠═b281dcf0-80f4-11eb-39e3-752784ef95ed
# ╠═c3867380-80f4-11eb-1a1e-61892538c262
# ╠═d78fbce2-80f4-11eb-3a15-09eb96356fec
# ╠═0f647d7e-80f5-11eb-13a0-9b1608eb924a
# ╠═07d2deac-80f5-11eb-21a9-f3d7da20e7b6
# ╠═386472f6-80f5-11eb-01e2-1b894c25ce58
# ╠═55d89790-80f5-11eb-2dac-b1a2a371b578
# ╠═4b82bb72-80f5-11eb-0c31-c37caf088c8e
# ╠═689cca2c-80f5-11eb-0f58-6bc0f7ab6acb
# ╠═79179e2c-80f5-11eb-0254-c7d9086d7a31
# ╠═9e73856e-80f5-11eb-2a09-eb3499ab0bc5
# ╠═f2c0a914-80f5-11eb-144c-99c87d037490
# ╠═12a86daa-80f6-11eb-3ff4-a7d9bbf0f8a6
