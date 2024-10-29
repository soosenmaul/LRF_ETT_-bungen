#* Übung 0
# In dieser Übung werden die wichtigsten Befehle und Vorgehen geübt.

# Laden der notwendigen Packages, s.u.
using Parameters
using LinearAlgebra
using NLsolve
# using DifferentialEquations
using Plots
using NumericalIntegration
plotly()                    # wähle den plotly als plot backend

include("UE0_fun.jl") #--> dort werden alle Hilfsfunktionen gesammelt

# Dokumentation:
# https://docs.julialang.org/en/v1/

##! Aufgabe 0.1: Löse ein lineares Gleichungssystem
# Problem:
# Gl. (1)  a + 2 d = 1
# Gl. (2)  b + d   = 3
# In Matrix Vektor Schreibweise überführen:
# A * c = b
# Massen/ Koeffizienten/Problem Matrix A
A = [1.0 2.0 ; 1.0 1.0] 
# ACHTUNG kein Komma zwischen den Zahlen in den Spalten!
# Unbekannten Vektor
# c = [a, d]
# ACHTUNG bei Vektoren werden die Zeilen mit Komma abgetrennt
# Quelle Vektor
b = [1.0, 3.0]
# Lösen den linearen Gleichungssystems mit dem Paket LinearAlgebra

# using Pkg
# Pkg.add("LinearAlgebra")
x = A\b # \b bedeutet: Teilen von links also x = A^-1*b

##! Aufgabe 0.2: Löse ein nicht lineare Gleichung
# Löse numerisch die Gleichung
# p = 10^(A-B/(C+T)) 
# nach T auf
# Antoine Parameter
A = 8.07
B = 1730.63
C = 233.426
# gewünschter Druck in 
p = 0.8 * 0.00133322^(-1) # bar * 0.00133322^(-1) mmHg/bar = mmHg




# Weil nicht linear benötigen wir einen ersten Ratewert
T_guess= [100.0] # °C !
# ACHTUNG T_guess MUSS ein Vektor sein!!! Hier nur mit einem Eintrag
# Lösen der Gleichung f=0 mit Hilfe des Pakets NLsolve
# https://github.com/JuliaNLSolvers/NLsolve.jl

# Dies ist die von nlsolve empfohlene Variante, da schneller (y ist vordefiniert)
sol = nlsolve((y,T)->EQ_psat!(y,T,A,B,C,p),T_guess) # (y,T)->EQ_psat!(y,T,A,B,C,p) bedeutet: y,T sind die Variablen, A,B,C,p sind Parameter, die aus dem Workspace übernommen werden
Tsol= sol.zero                                      # lese die Nullstelle aus dem Lösungsobjekt aus

# Alternative, wenn die Funtion y über return zurückgibt. Das ist die besser verständliche Methode. Theoretisch aber etwas langsamer, da y bei jedem Aufruf erst in der Funktion ertsellt wird.
sol = nlsolve((T)->EQ_psat2!(T,A,B,C,p),T_guess)    # löse die Gleichung f=0 mit dem Solver nlsolve
Tsol= sol.zero                                      # lese die Nullstelle aus dem Lösungsobjekt aus



##! Aufgabe 0.3 Löse ein nichtlineares Gleichungssystem

# Setzte den Ratewert
X_guess=[0.5;0.5;0.5]               # Vektor mit den Ratewerten für x,y,z
sol = nlsolve(EQ_NonLin!,X_guess)   # Löse das Gleichungssystem
X_sol= sol.zero                     # Lese die Nullstelle aus dem Lösungsobjekt aus


##! Aufgabe 0.4: Lösen einer ODE für ein zeitliches Problem
# Problem
# Abkühlung eines Tanks

# Stoff- und Prozessdaten
p = (;
cp      = 2200.0  ,         # J/(kg K)
k       = 700   ,           # W/(m^2 K s)
A       = 15.0  ,           # m^2
M       = 1.e4,             # kg
TUmg    = 20.0 + 273.15,    # K
T0      = 85.0 + 273.15,    # K
tspan   = (0.0,10000)       # s
)

# Definition des Problems
prob = ODEProblem(RHS_InstTank,p.T0,p.tspan,p)          # https://diffeq.sciml.ai/stable/
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)    # Tsit5 ist ein Standardsolver wie in Matlab ode45
                                                        # Liste von Solvern: https://diffeq.sciml.ai/stable/solvers/ode_solve/

# Zeit, an denen Julia die Lösung berechnet (wird vom Solver bestimmt)
t_sol = sol.t


p1 = plot(sol, label="num. Lösung")

# Darstellung der Lösung zu gewünschten Zeiten
t_plot = 0:100:8000
T_plot = sol(t_plot) # interpoliert die Lösung an den gewünschten Zeiten

# Analytische Lösung
T_analyt = T_Tank(t_sol,p)
plot!(p1, t_sol,T_analyt,xlabel="t / s", ylabel="T / K", label="Analytische Lösung")
gui(p1)

##! Aufgabe 0.5: Lösung eine ODE 
#! Problem: Wärmeleitung in ebener Wand mit Wärmequelle

# Stoff- und Prozessdaten
p = (;
λ      = 50.0  , # W/(m K)
Φ_0     = 5000   ,           # W/(m^3)
Ngrid   = 5,             # Anzahl Gitterpunkte
Tlinks  = 293.0,    # K
T_ref   = 300.0,    # K
L       = 1.0        # m
)

#! Variante 1: Diskretisierung des Ortes, d..h. ODE --> Nicht lineares Gleichungssystem
# Berechnung des Abstandes der Gitterpunkte
Delta_z = p.L/(p.Ngrid - 1)
z = 0:Delta_z:p.L

# Fügt das zu den Parametern hinzu --> es muss nicht in der ODE berechnet werden
p=merge(p,(;Δz=Delta_z))

# Berechnung der Ratwerte: --> Achtung: da hier links die Temp. vorgegeben wird, wird T_1 nicht mehr berechnet, d.h. wir brauchen nur Ngrid-1 Punkte
T_rate = ones(p.Ngrid-1) .+ p.Tlinks .+ 10.0

sol    = nlsolve((T)->EQ_WL_ebeneWand(T,p),T_rate,ftol=1e-6) # Löse das Gleichungssystem        

if ~sol.f_converged # Gibt eine Fehlermeldung raus, wenn der Solver nicht KONVERGIERT
    error("Löser ist nicht konvergiert. Versuche anderen Ratwert")
end
# Erstellen des Lösungsvektors: jetzt muss T_1 = T_links angehängt werden
T_plot = [p.Tlinks;sol.zero]
# Erstellen des Plots
p2     = plot(z, T_plot, xlabel = "z / m", ylabel = "T / K")        
display(p2)         # Anzeigen des Plots

# Überprüfung des integralen 1. Hauptsatzes
# Alles was über die Wärmequelle frei wird, muss über den linken Rand abgegeben werden
qpkt_ab_links = -p.λ * (T_plot[2] - T_plot[1])./Delta_z             # W/m2

# Berechnung der mittleren Temp. zwischen den Gitterpzunkten
T_m         = (T_plot[1:end-1] + T_plot[2:end])./2
Φ_m         = Phi(T_m,p.Φ_0,p.T_ref)
qpkt_quelle = sum(Φ_m)*Delta_z                              # W/m2
rel_fehler  = (qpkt_ab_links + qpkt_quelle)./qpkt_quelle

#! Variante 2: Umschreiben in ein system an DGLs, dazu Einführung von qpkt!
#! Das Problem ist ein Randwertproblem, Werte am linken und rechten Rand bekannt sind, aber nicht zwei RBs links (z=0)

Y_guess = [p.Tlinks, 0.0] # Ratwerte für die Lösung von bc()
zspan = (0.0, p.L) # Definiert z, 
bvp2 = BVProblem(DGL_WL_ebeneWand!, bc!, Y_guess,zspan,p)
sol_b = solve(bvp2, MIRK4(), dt = 0.01) # we need to use the MIRK4 solver for TwoPointBVProblem
Y = convert(Array, sol_b)
T    = Y[1,:]
qdot = Y[2,:]
z = sol_b.t

# bvp2 = TwoPointBVProblem(DGL_WL_ebeneWand!, (bc_links!, bc_rechts!) , Y_guess,zspan,p;
# bcresid_prototype = (zeros(1), zeros(1)))
# sol_b2 = solve(bvp2, MIRK4(), dt = 0.01) # we need to use the MIRK4 solver for TwoPointBVProblem
# Y2 = convert(Array, sol_b2)
# T2 = Y2[1,:]
# z2 = sol_b2.t

# Erstellen des Plots
plot!(p2, z, T, xlabel = "z / m", ylabel = "T / K")        
        # Anzeigen des Plots

        
        # Erstellen des Plots
#plot!(p2, z2, T2)
                # Anzeigen des Plots

p2b     = plot(z, qdot, xlabel = "z / m", ylabel = "qdot / W/m2") 
display(p2b)

# Überprüfen des Fehlers

Φ = Phi(T,p.Φ_0,p.T_ref)
qpkt_quelle_b = integrate(z,Φ)
rel_fehler_b  = (qdot[1] + qpkt_quelle_b)./qpkt_quelle_b


#! Aufgabe 0.6: Lösung einer PDE: instationäre Wärmeleitung
#* Lösung mit Finite Differenzen Methode
p = (;
λ      = 50.0  ,        # W/(m K)
c      = 500.,          # J/(kg K)
ρ      = 2500.,         # kg/m3
Φ_0     = 5000   ,      # W/(m^3)
Ngrid   = 100,          # Anzahl Gitterpunkte
Tlinks  = 293.0,        # K
T_ref   = 300.0,        # K
Tini    = 273.0,        # K
L       = 1.0,          # m
tspan   = (0, 2*3600),  # s
)

# Berechnung des Abstandes der Gitterpunkte
Delta_z = p.L/(p.Ngrid - 1)
z       = 0:Delta_z:p.L

# Fügt das zu den Parametern hinzu --> es muss nicht in der ODE berechnet werden
p=merge(p,(;Δz=Delta_z))

# Berechnung der Ratwerte: --> Achtung: da hier links die Temp. vorgegeben wird, wird T_1 nicht mehr berechnet, d.h. wir brauchen nur Ngrid-1 Punkte
Tini_vect = ones(p.Ngrid-1) .* p.Tini

prob = ODEProblem(RHS_WL_ebeneWand,Tini_vect,p.tspan,p)  
sol = solve(prob, QNDF(), reltol=1e-6, abstol=1e-6) # Alternativer Solver bei Method of Line
# In den Zeilen: Ort, In den Spalten: Zeit

t_plot = (0:20:120).*60 # s

# hinzu fügen der Randbedingung links
T_plot       = zeros(p.Ngrid, length(t_plot))    # Erstelle eine Matrix mit Ngrid Zeilen und length(t_plot) Spalten
sol_plot     = sol(t_plot)                       # Interpoliere die Lösung an den gewünschten Zeiten
T_plot[1,:]  = ones(1,length(t_plot)).*p.Tlinks  # Setze die Randbedingung links
T_plot[2:end,:] = convert(Array,sol_plot)                 # sieht komisch aus. Ist aber notwendig, damit es in die Matrix T_plot passt


# Erstellen des Plots
p3              = plot(z, T_plot, xlabel = "z / m", ylabel = "T / K")

#* Lösung mit ModelingToolkit (symbolisch)
using ModelingToolkit, MethodOfLines, OrdinaryDiffEq, LinearSolve, DomainSets

# Symbole des "symbolic PDE-System" definieren
@parameters x t     # Raum: x, Zeit: t
@variables T(..)    # Temperatur als Funktion von x und t
Dt = Differential(t)
Dx = Differential(x)
Dxx = Differential(x)^2

# Bereiche definieren
x_min = t_min = 0.0
x_max = 1.0
t_max = 3600.0 * 2

domains = [x ∈ Interval(x_min, x_max), t ∈ Interval(t_min, t_max)]

α = p.λ / (p.ρ * p.c)                  # m2/s

phi(x,t) = p.Φ_0 * (T(x, t) / p.T_ref)^0.5    # W/(m^3)
# Anfangsbedingungen
T_ini(x, t) = 273.0   # K

# Wärmeleitgleichung in symbolischer Form
eq = Dt(T(x, t)) ~ α * Dxx(T(x, t)) + phi(x,t) /(p.ρ * p.c)

# Randbedingungen
bcs = [T(x, 0) ~ T_ini(x, 0), T(0, t) ~ p.Tlinks, Dx(T(1, t)) ~ 0]

# PDE System definieren
@named pdesys = PDESystem(eq, bcs, domains, [x, t], [T(x, t)])

# Automatisierte symbolische Diskretisierung mit MethodOfLines
N = 100;     # Anzahl der Gitterpunkte
dx = (x_max - x_min)/N 
order = 2   # Geneuigkeit der Approximation
discretization = MOLFiniteDifference([x => dx], t, approx_order = order, grid_align = center_align)

# System diskretisieren, aus PDE in ODE Problem
prob = discretize(pdesys, discretization)

# Lösen der ODE -> Standardverfahren über ODE solver
sol = solve(prob, TRBDF2(), saveat = 0.1)
#sol = solve(prob, QNDF(), reltol=1e-6, abstol=1e-6)

# Ergebnisse über das "Symbolic Solution Interface"
discrete_x = sol[x]
discrete_t = sol[t]

solT = sol[T(x, t)]

# Plotten der Ergebnisse
using Plots

num_plots = 6  # Anzahl der Plots über den Zeitraum t_max
time_step = ceil(Int, length(discrete_t)/num_plots)
time_indices = 1:time_step:length(discrete_t)
plot(discrete_x, solT[:, time_indices], xlabel="Wand", ylabel="Temperatur", legend=:bottomleft)