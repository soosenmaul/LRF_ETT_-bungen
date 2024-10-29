#! Funktionen für Übung 0

#? Aufgabe 0.2
# Aufstellen der Lösungsfunktion in der Form f=0
function EQ_psat!(y,T,A,B,C,p) # Allgemeine Form Fun!(y,x, param1, param2, param3,...)
    # y soll Null sein! --> Der Solver sucht Nullstellen!!!
    y.= p .- 10.0.^(A.-B./(C.+T)) # Hier ist .= wichtig => damit weißt Julia die rechte Seite y zu. Bei y = ... wird ein "neues" y kreiiert. Das fürht zu einem Fehler
end

# Alternative Form
function EQ_psat2!(T,A,B,C,p)
    # y soll Null sein!
    y = p .- 10.0.^(A.-B./(C.+T)) 
    return y
end

#? Aufgabe 0.3
# Problem
# 2 x + y + 2 z^2   = 5
#     y^2 + 4 z = 4
#     x*y + z = exp^(z) 
# Definiere die Problem funktion
function EQ_NonLin!(Y,X)
    x = X[1]
    y = X[2]
    z = X[3]
    Y[1] = 2.0 * x + y + 2* z^2 - 5.0
    Y[2] = y^2 + 4* z - 4.0
    Y[3] = x*y +z - exp(z)
end


#? Aufgabe 0.4 
# Rechte Seite der Differentialgleichung 
function RHS_InstTank(T,p,t)
    # unpack die Parameter Kontainer
    @unpack cp, k, A, M, TUmg  = p
    # berechne den Wärmestrom
    Qab  = k*A*(T - TUmg)  # W
    # stelle die ODE in der Form dy/dt = ... auf 
    return dTdt = -1.0/(cp*M)*Qab  #  K/s
end

function T_Tank(t,p)
    # Analytische Lösung
    @unpack  cp, k, A, M, TUmg, T0  = p
    C = k*A/(M*cp)
    T = TUmg .+ (T0 .- TUmg) .* exp.(-C.*t) # .-Operator wichtig, damit hier mit Vektoren gerechnet werden kann
end


#? Aufgabe 0.5
# Stationäre Wärmeleitung in der Wand
function EQ_WL_ebeneWand(T,p)    
    @unpack λ, Φ_0, Tlinks, T_ref, Ngrid, Δz  = p # Griechische Bustaben werden wie in LaTeX über \Delta erzeugt

    # Lösung der Punkte i=2...Ngrid, da links Temp. gegeben, d.h. length(Y) = Ngrid-1
    # d.h. T = [T_2; T_3; ...; T_N] => length(T) = Ngrid-1

    # Berechnung des Dummypunktes N+1
    qpkt_rechts = 0.0 # adiabatischer Randbedingung
    T_Nplus = T[end] - qpkt_rechts*Δz/λ # Hier: T_N+1 = T_N
    
    T_erw = [Tlinks;T;T_Nplus] # T an den Stellen i=1... N+1

    Φ = Phi(T,Φ_0,T_ref) # Berechnung von Phi an allen Gitterpunkten abgesehen von T_1

    Y = λ/Δz^2 .*(T_erw[3:end] .- 2.0 .* T_erw[2:end-1] .+ T_erw[1:end-2]) .+ Φ

    return Y
end
# Stationäre Wärmeleitung in der Wand
function DGL_WL_ebeneWand!(dY_dz, Y,p,t)  

    T = Y[1]
    qdot = Y[2]

    Φ = Phi(T,p.Φ_0,p.T_ref)

    # Das ist äquivalent zu 0 = λ d^2T/dz^2 + Φ
    dT_dz = - qdot ./ p.λ
    dq_dz = Φ 

    dY_dz .= [dT_dz, dq_dz]
end

# Beschreibung der RBs für die stationäre WL in der Wand
function bc!(residual, Y, p, t)
    # RB ist erfüllt, wenn die residuen = 0 sind!

    residual[1] = Y[1][1] - p.Tlinks # T=Tlinks at z= 0 Y[Ort][Variable]
    residual[2] = Y[end][2]  # adiabatic boundary at qdot(z= L)= 0
end

function bc_links!(resi0, Y0, p)
    # RB ist erfüllt, wenn die residuen = 0 sind!
    resi0[1] = Y0[1] - p.Tlinks # T=Tlinks at z= 0 Y[Ort][Variable]
    
end

function bc_rechts!(resi_end, Yend, p)
    # RB ist erfüllt, wenn die residuen = 0 sind!
    resi_end[1] = Yend[2]  # adiabatic boundary at qdot(z= L)= 0
end

#? Aufgabe 0.6
# Instationäre Wärmeleitung in der Wand
function RHS_WL_ebeneWand(dTdt,T,p,t)    
    @unpack λ, c, ρ, Φ_0, Tlinks, T_ref, Ngrid, Δz  = p # Griechische Bustaben werden wie in LaTeX über \Delta erzeugt

    # Lösung der Punkte i=2...Ngrid, da links Temp. gegeben, d.h. length(Y) = Ngrid-1
    # d.h. T = [T_2; T_3; ...; T_N] => length(T) = Ngrid-1

    # Berechnung des Dummypunktes N+1
    qpkt_rechts = 0.0                   # adiabatischer Randbedingung
    T_Nplus = T[end] - qpkt_rechts*Δz/λ # Hier: T_N+1 = T_N
    
    T_erw = [Tlinks;T;T_Nplus]          # T an den Stellen i=1... N+1

    Φ = Phi(T,Φ_0,T_ref)                # Berechnung von Phi an allen Gitterpunkten abgesehen von T_1

    dTdt .= (λ/Δz^2 .*(T_erw[3:end] .- 2.0 .* T_erw[2:end-1] .+ T_erw[1:end-2]) .+ Φ)./(ρ*c)

end

# Hilfsfunktion für Φ
function Phi(T,Φ_0,T_ref)
    Φ = Φ_0.*(T./T_ref).^(0.5)
    return Φ
end