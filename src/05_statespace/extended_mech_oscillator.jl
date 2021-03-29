### A Pluto.jl notebook ###
# v0.12.12

using Markdown
using InteractiveUtils

# ╔═╡ 8acd485a-88d0-11eb-3c08-630c88acf4c6
using Luxor

# ╔═╡ 69119f16-88cd-11eb-3c64-23d024729ffa
using LinearAlgebra, DifferentialEquations, Plots;

# ╔═╡ 128e3a06-88cc-11eb-0925-a79101048c01
md"# State-Space Modeling: Extended Mechanical Oscillator

A mechanical oscillator with two spring-damping-mass series connected subsystems (see Figure below) are examined. 

At both subsystems $i \in \left\{1,2\right\}$ operate the forces of the springs

$F_{c,i}(t) = c_{i} \left[ y_{i}(t) - y_{i-1}(t) \right] \text{,}$

of the attenuators (damping)

$F_{c,i}(t) = d_{i} \left[ \dot{y}_{i}(t) - \dot{y}_{i-1}(t) \right]$

and of the mass

$F_{m,i}(t) = m_{i} \ddot{y}_{i}(t) \text{.}$

Here, it is assumed that the left boundary is not movable and thus

$y_{0}(t) = \dot{y}_{0}(t) \equiv 0 \text{.}$

Two external forces $F_{ex,i}$ are applied at the masses and the resulting sum of forces are

$F_{ex,1}(t) ~=~ F_{m,1}(t) + F_{c,1}(t) + F_{d_1}(t) - F_{c,2}(t) - F_{d,2}(t)$

$=~ m_{1} \ddot{y}_{1}(t) + c_{i} y_{1}(t) + d_{1} \dot{y}_{1}(t) -  c_{2} \left[ y_{2}(t) - y_{1}(t) \right] - d_{2} \left[ \dot{y}_{2}(t) - \dot{y}_{1}(t) \right]  \tag{1}$

and

$F_{ex,2}(t) ~=~ F_{m,2}(t) + F_{c,2}(t) + F_{d_2}(t)$

$=~ m_{2} \ddot{y}_{2}(t) + c_{2} \left[ y_{2}(t) - y_{1}(t) \right] + d_{2} \left[ \dot{y}_{2}(t) - \dot{y}_{1}(t) \right] \text{.} \tag{2}$

The external forces impact a change of positions $y_{i}(t)$ which are measured by a sensor. That means, the positions $y_{i}(t)$ are the outputs of the system."

# ╔═╡ 8da2d48c-88d0-11eb-3d0c-d56d38a43dca
@svg begin
	
	# Wall
	line( Point(-250, -180) , Point(-250,0), :stroke ) 
	
	for i = 0 : 8
		dx = i*20
		line( Point(-252, -170+dx) , Point(-275,-160+dx), :stroke ) 
	end
	
	
	##### First System ########
	
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
	# Luxor.arrow( Point(15, -90), Point(100, -90), linewidth=3.0)
	
	# Sum of forces
	
	line( Point(-20, 100), Point(-20, 230), :stroke) 
	Luxor.arrow( Point(-20, 110), Point(-80, 110), linewidth=2.0)
	Luxor.arrow( Point(-20, 160), Point(-80, 160), linewidth=2.0)
	Luxor.arrow( Point(-20, 210), Point(-80, 210), linewidth=2.0)
	Luxor.arrow( Point(-20, 120), Point( 40, 120), linewidth=2.0)
	Luxor.arrow( Point(-20, 170), Point( 40, 170), linewidth=2.0)
	Luxor.arrow( Point(-20, 220), Point( 40, 220), linewidth=2.0)
	
	
	
	
	##### Second System ########

	DX = 260;
	
	# Attenuator
	line( Point(-250+DX, -140) , Point( -200+DX, -140), :stroke) 
	line( Point(-200+DX, -160) , Point( -200+DX, -120), :stroke) 
	line( Point(-200+DX, -160) , Point( -120+DX, -160), :stroke) 
	line( Point(-200+DX, -120) , Point( -120+DX, -120), :stroke)
	line( Point(-140+DX, -155) , Point( -140+DX, -125), :stroke) 
	line( Point(-140+DX, -140) , Point(  -50+DX, -140), :stroke) 
	
	
	# Spring
	line( Point( -250+DX, -40) , Point( -200+DX, -40), :stroke) 
	line( Point( -200+DX, -40) , Point( -190+DX, -60), :stroke)	
	
	for i = 0 : 4
	dx = i*20 +DX;
	line(  Point( -190+dx, -60), Point( -180+dx, -20), :stroke)	
	line(  Point( -180+dx, -20), Point( -170+dx, -60), :stroke)
	end
	line( Point( -90+DX, -60), Point( -80+DX, -40), :stroke)
	line( Point( -80+DX, -40), Point( -50+DX, -40), :stroke)
	
	# Mass
	rect(Point(-50+DX, -170), 60, 160, :stroke)

	
	#Sum of forces
	line( Point(-20+DX, 100), Point(-20+DX, 220), :stroke) 
	Luxor.arrow( Point(-20+DX, 110), Point(-80+DX, 110), linewidth=2.0)
	Luxor.arrow( Point(-20+DX, 160), Point(-80+DX, 160), linewidth=2.0)
	Luxor.arrow( Point(-20+DX, 210), Point(-80+DX, 210), linewidth=2.0)
	Luxor.arrow( Point(-20+DX, 170), Point( 40+DX, 170), linewidth=2.0)
	
	
	
	# Coordinates
	line(Point(-250,10), Point(-250,40), :stroke)
	Luxor.arrow(Point(-250,25), Point(-220,25), linewidth=2.0) 
	
	line(Point(-20,10), Point(-20,40), :stroke)
	Luxor.arrow(Point(-20,25), Point(10,25), linewidth=2.0) 
	
	line(Point(-20+DX,10), Point(-20+DX,40), :stroke)
	Luxor.arrow(Point(-20+DX,25), Point(10+DX,25), linewidth=2.0) 

	fontsize(20)
	Luxor.text("y", Point(-240, 45))	
	Luxor.text("y", Point(-15, 45))
	Luxor.text("y", Point(-15+DX, 45))

	fontsize(15)
	Luxor.text("1", Point(-5, 52))
	Luxor.text("2", Point(-5+DX, 52))
	
	
	# Text
	
	fontsize(20)
	Luxor.text("Attenuator d", Point(-210, -180))
	Luxor.text("Mass", Point(-50, -200))
	Luxor.text("m", Point(-40, -180))
	Luxor.text("Spring c", Point(-210, -70))
	
	
	Luxor.text("Attenuator d", Point(-210+DX, -180))
	Luxor.text("Mass", Point(-50+DX, -200))
	Luxor.text("m", Point(-40+DX, -180))
	Luxor.text("Spring c", Point(-210+DX, -70))
	
	fontsize(15)
	Luxor.text("1", Point(-85, -175))
	Luxor.text("1", Point(-20, -175))
	Luxor.text("1", Point(-130, -65))
	
	
	Luxor.text("2", Point(-85+DX, -175))
	Luxor.text("2", Point(-20+DX, -175))
	Luxor.text("2", Point(-130+DX, -65))
	
	#Luxor.text("external ", Point(40, -120))
	#Luxor.text("Force F", Point(40, -100))
	#fontsize(15)
	#Luxor.text("ex", Point(110, -95))
	
	#fontsize(25)
	#Luxor.text("Sum of forces", Point(-100, 70))
	
	fontsize(20)
	Luxor.text("F", Point(-60, 100))
	Luxor.text("F", Point(-60, 150))
	Luxor.text("F", Point(-60, 200))
	Luxor.text("F", Point( 0, 110))
	Luxor.text("F", Point( 0, 160))
	Luxor.text("F", Point( 0, 210))
	
	Luxor.text("F", Point(-60+DX, 100))
	Luxor.text("F", Point(-60+DX, 150))
	Luxor.text("F", Point(-60+DX, 200))
	Luxor.text("F", Point( 0+DX, 160))
	
	
	fontsize(15)
	Luxor.text("d,1", Point(-50, 105))
	Luxor.text("c,1", Point(-50, 155))
	Luxor.text("m,1", Point(-50, 205))
	Luxor.text("ex,1", Point(10, 165))
	Luxor.text("d,2", Point(10, 115))
	Luxor.text("c,2", Point(10, 215))
	
	
	Luxor.text("d,2", Point(-50+DX, 105))
	Luxor.text("c,2", Point(-50+DX, 155))
	Luxor.text("m,2", Point(-50+DX, 205))
	Luxor.text("ex,2", Point(10+DX, 165))
	
end

# ╔═╡ 5f55af18-88cc-11eb-2d8a-efe679337085
md"## Differential Equation

The differential equation is built as a state-space model. To do so, a n-th order differential equation is transfered to a system of first order differential equations with n states. Here, the two second-order differential equations are reshaped as one system of first order differential equations with four states. 

The external forces are renamed as inputs
$u_{1}(t) ~=~ F_{ex,1}(t) \quad \text{and} \quad  u_{2}(t) ~=~ F_{ex,2}(t) \text{.}$

The outputs or measured positions $y_{i}(t)$ are renamed as

$x_{1}(t) := y_{1}(t) \qquad \text{and} \qquad x_{2}(t) := y_{2}(t)$

and new variables are introduced with

$x_{3}(t) ~=~ \dot{y}_{1}(t) ~=~ \dot{x}_{1}(t) \text{,}$

$x_{4}(t) ~=~ \dot{y}_{2}(t) ~=~ \dot{x}_{2}(t) \text{.}$

The sum of forces $(1)$ and $(2)$ are rewritten with the new variables as

$\dot{x}_{3}(t) ~=~ -\frac{c_{1}}{m_{1}} x_{1}(t) - \frac{d_{1}}{m_{1}} x_{3}(t) +  \frac{c_{2}}{m_{1}} \left[ x_{2}(t) - x_{1}(t) \right] +  \frac{d_{2}}{m_{1}} \left[ x_{4}(t) - x_{3}(t) \right] + \frac{1}{m_{1}} u_{1}(t)$

and 

$\dot{x}_{4}(t) ~=~ - \frac{c_{2}}{m_{2}} \left[ x_{2}(t) - x_{1}(t) \right] -  \frac{d_{2}}{m_{2}} \left[ x_{4}(t) - x_{3}(t) \right] + \frac{1}{m_{2}} u_{2}(t) \text{.}$

Finally, the variables are reorganized in a matrix-vector notation to gain the state-space model

$\begin{pmatrix}
\dot{x}_{1} \\
\dot{x}_{2} \\
\dot{x}_{3} \\
\dot{x}_{4}
\end{pmatrix}
=
\begin{pmatrix}
0 & 0 & 1 & 0 \\
0 & 0 & 0 & 1 \\
\frac{-1}{m_{1}}\left[c_{1} + c_{2}\right] & \frac{c_{2}}{m_{1}} & \frac{-1}{m_{1}}\left[d_{1} + d_{2}\right] & \frac{d_{2}}{m_{1}} \\
\frac{c_{2}}{m_{2}} & -\frac{c_{2}}{m_{2}} & \frac{d_{2}}{m_{2}} & -\frac{d_{2}}{m_{2}}
\end{pmatrix}
\begin{pmatrix}
x_{1} \\
x_{2} \\
x_{3} \\
x_{4}
\end{pmatrix}
+
\begin{pmatrix}
0 & 0 \\
0 & 0 \\
\frac{1}{m_{1}} & 0 \\
0 &\frac{1}{m_{2}}
\end{pmatrix}
\begin{pmatrix}
u_{1}(t) \\
u_{2}(t)
\end{pmatrix}$

with states $x_{i}$ and output

$\begin{pmatrix}
y_{1}(t) \\
y_{2}(t)
\end{pmatrix}
=
\begin{pmatrix}
1 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 
\end{pmatrix}
\begin{pmatrix}
x_{1}(t) \\
x_{2}(t) \\
x_{3}(t) \\
x_{4}(t)
\end{pmatrix}
\text{.}$

Linear time-invariant (LTI) systems are noted in the standard form as

$\dot{x}(t) = A ~ x(t) + B ~ u(t)$ 

$y(t) = C ~ x(t) + D ~ u(t)$

in which the matrices correspond here to 

$A = 
\begin{pmatrix}
0 & 0 & 1 & 0 \\
0 & 0 & 0 & 1 \\
\frac{-1}{m_{1}}\left[c_{1} + c_{2}\right] & \frac{c_{2}}{m_{1}} & \frac{-1}{m_{1}}\left[d_{1} + d_{2}\right] & \frac{d_{2}}{m_{1}} \\
\frac{c_{2}}{m_{2}} & -\frac{c_{2}}{m_{2}} & \frac{d_{2}}{m_{2}} & -\frac{d_{2}}{m_{2}}
\end{pmatrix} 
\qquad \text{,} \qquad
B = 
\begin{pmatrix}
0 & 0 \\
0 & 0 \\
\frac{1}{m_{1}} & 0 \\
0 &\frac{1}{m_{2}}
\end{pmatrix}$

and

$C = 
\begin{pmatrix}
1 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 
\end{pmatrix}
\qquad \text{,} \qquad
D = 0_{2 \times 2}
\text{.}$


## Stability

The system without control 

$\dot{x}(t) = A x(t)$

is stable if **all** Eigenvalues $\lambda_{i}$ of A are smaller than zero. The Eigenvalues are found by calculating

$\lambda ~ x(t) = A ~ x(t) \quad \Rightarrow \quad (\lambda I - A) ~ x(t) = 0$

and solving the determinant

$\det(\lambda I - A) = 0 \text{.}$

## Simulation

The mechanical oscillator with two masses is simulated next. Firstly, the physical constants for mass, damping and spring, and the matrices of the dynamical system have to be specified.
"

# ╔═╡ 6c32023a-88cd-11eb-0d6a-ff6dd06e2214
begin
	const d₁ = 0.5; # Damping
	const d₂ = 0.5;
	const c₁ = 0.1; # Spring constant
	const c₂ = 0.1;
	const m₁ = 1.0; # Mass
	const m₂ = 2.0;
end

# ╔═╡ 6bfb076c-88cd-11eb-1471-a943ecebebc6
# System matrix
A = [0 0 1 0; 0 0 0 1; (-1/m₁)*(c₁ + c₂) c₂/m₁ (-1/m₁)*(d₁ + d₂) d₂/m₁; c₂/m₂ -c₂/m₂ d₂/m₂ -d₂/m₂]

# ╔═╡ 6bce1c52-88cd-11eb-3c5b-55bd467ae59a
B = [0 0; 0 0; 1/m₁ 0; 0 1/m₂] # Input matrix

# ╔═╡ 6ba028f6-88cd-11eb-1e9c-95fb8366b120
C = [1 0 0 0; 0 1 0 0] # Output matrix

# ╔═╡ a9e5579e-88cd-11eb-1dad-013e29cc3e63
md"### Stability

Next, the stability is proved to guarantee a suitable behaviour."

# ╔═╡ a9afcae0-88cd-11eb-3e9c-631622f21fd9
ev = eigvals(A) # Eigenvalues

# ╔═╡ c2468bbe-88cd-11eb-23e5-399cec37223b
md"The real part of all Eigenvalues is smaller than zero, thus the uncontrolled system is stable. Furthermore, the imaginary part of some Eigenvalues is not zero and so the system dynamics tend to oscillating behaviour.

### Defining the Ordinary Differential Equation

The ordinary differential equation is defined as a Julia function, the initial value $x_{0} = \left(1, 2, 0, 0\right)^{\top}$ and the simulation time range $t \in \left[ 0, 100 \right]$ are set, and the ODE prolem is built. Here, the system has an input of two scaled step functions

$u_{1}(t) = 0.5 \quad \text{for} ~ t \geq 0$

and

$u_{2}(t) = 1.5 \quad \text{for} ~ t \geq 0 \text{.}$"

# ╔═╡ a97ab478-88cd-11eb-0010-f3e4d25003cd
# Definition of ODE
function mech_oscillator(dx,x,p,t)
  
  u1 = 0.5; # 1. Control input
  u2 = 1.5; # 2. Control input
  u = [u1; u2]

  dx .= A*x + B*u # Right-hand side of ODE
end

# ╔═╡ a9470d78-88cd-11eb-22ad-7991d1d625a9
x₀ = [1.0; 2.0; 0.0; 0.0]; # Initial values

# ╔═╡ a874973a-88cd-11eb-2e30-8913f2b2a47d
tspan = (0.0, 100.0); # Time span

# ╔═╡ f804df1c-88cd-11eb-2f3f-b703f5566953
#Build ODE Problem
prob = ODEProblem(mech_oscillator,x₀,tspan, A);

# ╔═╡ fe353364-88cd-11eb-14af-0d1459b2971d
md"### Results

Finally, the ODE problem is solved and the output is plotted."

# ╔═╡ 5f779748-88ce-11eb-1ba1-1f9e9f4a2963
sol = solve(prob); # Solve ODE Problem

# ╔═╡ 616c49f4-88ce-11eb-18b1-d7e36edd68d1
y = C * sol; # Calculate system output

# ╔═╡ 68c851ac-88ce-11eb-0b8d-bb9395f79fce
plot(sol.t, transpose(y), title="System response", xaxis="Time [s]", yaxis="Position [m]")

# ╔═╡ 6cf0c5b0-88cf-11eb-1483-2bd0f4bb69b6
md"## Eigenvalue decomposition"

# ╔═╡ ad8b2da0-88ce-11eb-055d-b92e0f7b0a12
eigvals(A)

# ╔═╡ b1020e40-88ce-11eb-3926-1759fdbb5b29
evecs = eigvecs(A)

# ╔═╡ d9c5f364-88ce-11eb-0632-936ce7c1a5e5
Ã = inv(evecs) * A * evecs

# ╔═╡ cfc4e5d6-88cf-11eb-1c4f-8db0d5a49cd8
round.(real(Ã), digits=4)

# ╔═╡ 06c7076a-88d0-11eb-2f68-c1f58512bfc8
round.(imag(Ã), digits=4)

# ╔═╡ Cell order:
# ╟─128e3a06-88cc-11eb-0925-a79101048c01
# ╠═8acd485a-88d0-11eb-3c08-630c88acf4c6
# ╠═8da2d48c-88d0-11eb-3d0c-d56d38a43dca
# ╟─5f55af18-88cc-11eb-2d8a-efe679337085
# ╠═69119f16-88cd-11eb-3c64-23d024729ffa
# ╠═6c32023a-88cd-11eb-0d6a-ff6dd06e2214
# ╠═6bfb076c-88cd-11eb-1471-a943ecebebc6
# ╠═6bce1c52-88cd-11eb-3c5b-55bd467ae59a
# ╠═6ba028f6-88cd-11eb-1e9c-95fb8366b120
# ╠═a9e5579e-88cd-11eb-1dad-013e29cc3e63
# ╠═a9afcae0-88cd-11eb-3e9c-631622f21fd9
# ╠═c2468bbe-88cd-11eb-23e5-399cec37223b
# ╠═a97ab478-88cd-11eb-0010-f3e4d25003cd
# ╠═a9470d78-88cd-11eb-22ad-7991d1d625a9
# ╠═a874973a-88cd-11eb-2e30-8913f2b2a47d
# ╠═f804df1c-88cd-11eb-2f3f-b703f5566953
# ╠═fe353364-88cd-11eb-14af-0d1459b2971d
# ╠═5f779748-88ce-11eb-1ba1-1f9e9f4a2963
# ╠═616c49f4-88ce-11eb-18b1-d7e36edd68d1
# ╠═68c851ac-88ce-11eb-0b8d-bb9395f79fce
# ╠═6cf0c5b0-88cf-11eb-1483-2bd0f4bb69b6
# ╠═ad8b2da0-88ce-11eb-055d-b92e0f7b0a12
# ╠═b1020e40-88ce-11eb-3926-1759fdbb5b29
# ╠═d9c5f364-88ce-11eb-0632-936ce7c1a5e5
# ╠═cfc4e5d6-88cf-11eb-1c4f-8db0d5a49cd8
# ╠═06c7076a-88d0-11eb-2f68-c1f58512bfc8
