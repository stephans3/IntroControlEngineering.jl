### A Pluto.jl notebook ###
# v0.12.12

using Markdown
using InteractiveUtils

# ╔═╡ e7ed06c2-8265-11eb-3811-d55107701acd
using Luxor

# ╔═╡ 7f7b85be-8256-11eb-0234-4bc7abf5cb43
using DifferentialEquations

# ╔═╡ 8dfd45b4-8256-11eb-3c44-ff4a2b515be6
using Plots

# ╔═╡ d9b2de44-8254-11eb-1df3-4bba4c0cc785
md"# Mechanical oscillator

Consider a mechanical oscillator including a spring 

$F_{c}(t) = c ~ x(t)$

with $ c > 0 $, a damping (or attenuation)

$F_{D}(t) = D ~ v(t) = D ~ \dot{x}(t)$

with D > 0, and a mass

$F_{m}(t) = m ~ a(t) = m ~ \ddot{x}(t)$

with m > 0. This system is excited by an external force $ F_{ex}(t) $ and all forces are summarized to

$F_{m}(t) + F_{D}(t) + F_{c}(t) = m ~ \ddot{x}(t) + D ~ \dot{x}(t) + c ~ x(t) = F_{ex}(t) \text{.}$

The complete mechanical system is portrayed in the figure below. "

# ╔═╡ 50d0f1a8-8266-11eb-1981-b16cf8dcde7d
@svg begin
	
	# Wall
	line( Point(-250, -180) , Point(-250,0), :stroke ) 
	
	for i = 0 : 8
		dx = i*20
		line( Point(-252, -170+dx) , Point(-275,-160+dx), :stroke ) 
	end
	
	# Attenuator
	line( Point(-250, -140) , Point( -200, -140), :stroke) 
	line( Point(-200, -160) , Point( -200, -120), :stroke) 
	line( Point(-200, -160) , Point( -120, -160), :stroke) 
	line( Point(-200, -120) , Point( -120, -120), :stroke)
	line( Point(-140, -155) , Point( -140, -125), :stroke) 
	line( Point(-140, -140) , Point(  -50, -140), :stroke) 
	
	
	# Spring
	line( Point(-250, -40) , Point( -200, -40), :stroke) 
	line( Point( -200, -40), Point( -190, -60), :stroke)	
	
	for i = 0 : 4
	dx = i*20;
	line(  Point( -190+dx, -60), Point( -180+dx, -20), :stroke)	
	line(  Point( -180+dx, -20), Point( -170+dx, -60), :stroke)
	end
	line( Point( -90, -60), Point( -80, -40), :stroke)
	line( Point( -80, -40), Point( -50, -40), :stroke)
	
	# Mass
	rect(Point(-50, -170), 60, 160, :stroke)
	
	# External Force
	Luxor.arrow( Point(15, -90), Point(100, -90), linewidth=3.0)
	
	# Sum of forces
	
	line( Point(-20, 100), Point(-20, 220), :stroke) 
	Luxor.arrow( Point(-20, 110), Point(-80, 110), linewidth=2.0)
	Luxor.arrow( Point(-20, 160), Point(-80, 160), linewidth=2.0)
	Luxor.arrow( Point(-20, 210), Point(-80, 210), linewidth=2.0)
	Luxor.arrow( Point(-20, 170), Point( 40, 170), linewidth=2.0)
	
	fontsize(20)
	Luxor.text("Attenuator d", Point(-210, -180))
	Luxor.text("Mass m", Point(-60, -180))
	Luxor.text("Spring c", Point(-210, -70))
	Luxor.text("external ", Point(40, -120))
	Luxor.text("Force F", Point(40, -100))
	fontsize(15)
	Luxor.text("ex", Point(110, -95))
	
	fontsize(25)
	Luxor.text("Sum of forces", Point(-100, 70))
	
	fontsize(20)
	Luxor.text("F", Point(-60, 100))
	Luxor.text("F", Point(-60, 150))
	Luxor.text("F", Point(-60, 200))
	Luxor.text("F", Point( 0, 160))
	
	fontsize(15)
	Luxor.text("d", Point(-50, 105))
	Luxor.text("c", Point(-50, 155))
	Luxor.text("m", Point(-50, 205))
	Luxor.text("ex", Point(10, 165))
end

# ╔═╡ 5cf7202e-8255-11eb-02ea-d5d3a82e5532
md"The second order differential equation 

$\ddot{x}(t) + \frac{D}{m} \dot{x}(t) + \frac{c}{m} x(t) ~=~ \frac{1}{m} F_{ex}(t)$

is transfered with

$\frac{D}{m} = 2 ~ d ~ \omega_{0} \quad \text{,} \quad \frac{c}{m} = \omega_{0}^{2} \quad \text{,} \quad \frac{1}{m} = K ~ \omega_{0}^{2} \quad \text{and} \quad F_{ex}(t) = u(t)$

to the general oscillation equation

$\begin{align}
\ddot{y}(t) + 2 d ~ \omega_{0} ~ \dot{y}(t) + \omega_{0}^{2} ~ y(t) ~=~ K ~ \omega_{0}^{2} ~ u(t) \text{.} \tag{1}
\end{align}$


Equation $(1)$ is noted with $x_{1}(t) = y(t)$ and $x_{2}(t) = \dot{y}(t)$ as the first order differential equation

$\begin{align}
\begin{pmatrix}
\dot{x}_{1}(t) \\
\dot{x}_{2}(t)
\end{pmatrix}
=
\begin{pmatrix}
0 & 1 \\
- \omega_{0}^{2} & - 2 d ~ \omega_{0}
\end{pmatrix}
\begin{pmatrix}
x_{1}(t) \\
x_{2}(t)
\end{pmatrix}
+
\begin{pmatrix}
0 \\
K ~ \omega_{0}^{2}
\end{pmatrix}
u(t)
\end{align}$

## Stability

The uncontrolled system (for $u(t) = 0$) is stable if matrix 

$A = \begin{pmatrix}
0 & 1 \\
- \omega_{0}^{2} & - 2 d ~ \omega_{0}
\end{pmatrix}$

has only eigenvalues in the left complex space. The eigenvalues are calculated with

$det(\lambda I - A) = \lambda ~ (\lambda + 2 d \omega_{0}) + \omega_{0}^2 = \lambda^2 + 2 d ~ \omega_{0} \lambda + \omega_{0}^{2} = 0$

and thus one holds
$\lambda = -d ~ \omega_{0} \pm j \omega_{0} ~ \sqrt{1 - d^{2}} \text{.}$

The uncontrolled system is always stable if $d > 0$ and $\omega_{0} > 0$, which is both guaranteed for usual mechanical systems. Furthermore, a damping $d < 0$ implies complex eigenvalues and leads to an oscillating behaviour of the solution trajectory $y(t)$."

# ╔═╡ c954f7b6-8255-11eb-3698-d1542be8171d
md"## Numerical solution

The controlled mechancial oscillator $ \dot{x}(t) = A ~ x(t) + b ~ u(t) $ is simulated with a numerical integration method instead of calculating the solution via eigenvalues and eigenvectors. For further information about numerical integration see [Wikipedia](https://en.wikipedia.org/wiki/Numerical_methods_for_ordinary_differential_equations) and the [DifferentialEquations documentation](http://docs.juliadiffeq.org/latest/index.html).

The initial values are defined with"

# ╔═╡ e6338ea6-8255-11eb-39bf-971e0f65efc5
x₀ = [2.0, 1.0]; # Inital values

# ╔═╡ 33ea52ce-8256-11eb-0b76-4359ccdf13d6
tspan = (0.0, 50.0); # Time range

# ╔═╡ 33cc4432-8256-11eb-1e55-4fbbaf4ed145
d = 0.5; # Damping

# ╔═╡ 33b32e7a-8256-11eb-050e-8d563af9ca35
ω₀ = 0.7; # Eigenfrequency

# ╔═╡ f50a8f72-8257-11eb-1b0f-610bf448009f
K = 0.3;

# ╔═╡ 3fc66b46-8256-11eb-3907-4bdbe54de2a9
param = [d, ω₀,K];

# ╔═╡ 55e70020-8256-11eb-0032-0bfbc025ab61
md"The right-hand side $ A ~ x(t) + b ~ u(t) $  with an arbitrary input - here $u(t) = sign( sin(t))$ is chosen as a periodic excitation - and gain $K$ is defined as a function. The input is generated by the external force in the sense of a mechanical oscillator."

# ╔═╡ 76ce0264-8258-11eb-09f1-015ac3a62369
input(t) = sign(sin(t)) # Input signal: periodic excitation 

# ╔═╡ 6599bad0-8256-11eb-11dd-cd286378c68a
function mechanical_oscillator(dx, x, p, t)
    
	d = p[1] # Damping
	w = p[2] # Eigenfrequency
	K = p[3] # Gain or amplification
	
	u = input(t) # Input
    
	A = [0 1; -w^2 -2*d*w] 	# State matrix
	B = [0 ; K * w^2]		# Input matrix
	
    dx .= A*x + B*u
end

# ╔═╡ 6f57d4e4-8256-11eb-1101-db772f1d076e
md"Differential equations are solved in two steps in Julia. Firstly, the mathematical problem is formulated and secondly the mathematical problem is solved with a *numerical integration* method. Tolerances are set to specify the precision of the solution.  "

# ╔═╡ 924419a2-8256-11eb-2469-61a450600ea3
mech_osc_problem = ODEProblem(mechanical_oscillator,x₀,tspan, param); # Build of the ODE problem

# ╔═╡ 9844b674-8256-11eb-3978-971c1085e7eb
mech_osc_solution = solve(mech_osc_problem,Tsit5(),reltol=1e-8,abstol=1e-8); # Solution of the ODE problem

# ╔═╡ 88dd11fe-8256-11eb-1085-850b60adee63
md"The calculated solution is figured in a plot."

# ╔═╡ a54218e2-8258-11eb-2fb8-73f072b0800a
u_signals = input.(mech_osc_solution.t)

# ╔═╡ a2a5520e-8256-11eb-2eb1-aff9cee99ea3
plot(mech_osc_solution.t, [mech_osc_solution[:,:]', u_signals], title="Mechanical Oscillator", label=["x1" "x2" "u"])

# ╔═╡ Cell order:
# ╟─d9b2de44-8254-11eb-1df3-4bba4c0cc785
# ╠═e7ed06c2-8265-11eb-3811-d55107701acd
# ╠═50d0f1a8-8266-11eb-1981-b16cf8dcde7d
# ╟─5cf7202e-8255-11eb-02ea-d5d3a82e5532
# ╟─c954f7b6-8255-11eb-3698-d1542be8171d
# ╠═e6338ea6-8255-11eb-39bf-971e0f65efc5
# ╠═33ea52ce-8256-11eb-0b76-4359ccdf13d6
# ╠═33cc4432-8256-11eb-1e55-4fbbaf4ed145
# ╠═33b32e7a-8256-11eb-050e-8d563af9ca35
# ╠═f50a8f72-8257-11eb-1b0f-610bf448009f
# ╠═3fc66b46-8256-11eb-3907-4bdbe54de2a9
# ╟─55e70020-8256-11eb-0032-0bfbc025ab61
# ╠═76ce0264-8258-11eb-09f1-015ac3a62369
# ╠═6599bad0-8256-11eb-11dd-cd286378c68a
# ╟─6f57d4e4-8256-11eb-1101-db772f1d076e
# ╠═7f7b85be-8256-11eb-0234-4bc7abf5cb43
# ╠═924419a2-8256-11eb-2469-61a450600ea3
# ╠═9844b674-8256-11eb-3978-971c1085e7eb
# ╟─88dd11fe-8256-11eb-1085-850b60adee63
# ╠═8dfd45b4-8256-11eb-3c44-ff4a2b515be6
# ╠═a54218e2-8258-11eb-2fb8-73f072b0800a
# ╠═a2a5520e-8256-11eb-2eb1-aff9cee99ea3
