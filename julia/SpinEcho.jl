"""
Nitrogen and carbon system. Simulation for the following
experiments:
    * Ramsey experiment
    * SpinEcho experiment with XYXYYX sequence
    * SpinEcho experiment with XYXYYX sequence and fixed τ
"""

using Plots
using LinearAlgebra
using Statistics
using QuantumInformation
using TensorOperations

#####################################################################################
# operators

# Electron operators
Sₓ = 1/2 * [0 1;1 0]
Sₓ = Sₓ ⊗ Matrix(1.0I,2,2)
S_y = im*0.5*[0 1;-1 0]
S_y = S_y ⊗ Matrix(1.0I,2,2)
S_z = [0 0;0 -1]
S_z = S_z ⊗ Matrix(1.0I,2,2)

# Nucleous operators
Iₓ = 1/2 * [0 1;1 0]
Iₓ = Matrix(1.0I,2,2) ⊗ Iₓ
I_y = 1/(2*im) * [0 1;-1 0]
I_y = Matrix(1.0I,2,2) ⊗ I_y
I_z = 1/2 * [1 0; 0 -1]
I_z = Matrix(1.0I,2,2) ⊗ I_z

# Density matrix
ρₙᵥ = [1 0; 0 0]
ρ₁₃ = [0.5 0; 0 0.5]
ρ = kron(ρₙᵥ, ρ₁₃)

#####################################################################################
# Parameters

# Energy between spin state 0 and +-1 (zfs in ground state)(MHz)
D = 2870
# Magnetic field (G)
B = 500
# gyromagnetic ratio of the NV center’s electron spin(MHz/G)
γₙᵥ = 2.8
# gyromagnetic ratio of the nucleus carbon spin(MHz/G)
γ₁₃ = 1.07e-3
Ω = 40

# Hyperfine
Axx = 50e-3
Ayy = 0
Azz = 10e-3

# Omega electron (MHz)
ω = D - γₙᵥ * B
#####################################################################################
# Hamiltonians
#To get a Hamiltonian with the energy levels in the diagonal (zfs = zero field splitting)
H_zfs = D * S_z^2
# Interactions between electron and nucleous. Zeeman= when +1, -1 and 0 are
# splittig, if not -1 and +1 are in the same position.
H_zeemanₙᵥ = γₙᵥ * B * S_z
H_zeeman₁₃ = γ₁₃ * B * I_z
# To rotate in X or Y axis, Omega bc is the "line of the wave"
H_driveₓ = Ω * Sₓ
H_drive_y = Ω * S_y
# ω means the different between state 1 and -1 when a magnetic flied is
# aplied, we multiply it by same Hamiltonian as H_zfs to represent the
# different spin states.
H_drive_freq = - ω * S_z^2
# Hyperfine Hamiltonian. Describes hyperfine interaction between the
# NV centre’s electron spin and the nitrogen’s nuclear spin
H_hyperfine = Axx * S_z * Iₓ + Ayy * S_z * I_y + Azz * S_z * I_z


# Hamiltonians of the system

# Full electron hamiltonian with both drives
H = H_zfs + H_zeemanₙᵥ + H_zeeman₁₃ + H_driveₓ + H_drive_y + H_hyperfine + H_drive_freq

# Electron hamiltonian with no drives
H_free = H_zfs + H_zeemanₙᵥ + H_zeeman₁₃ + H_hyperfine + H_drive_freq

# Electron hamiltonian with one drive
Hₓ = H_zfs + H_zeemanₙᵥ + H_zeeman₁₃ + H_driveₓ + H_hyperfine + H_drive_freq

H_y = H_zfs + H_zeemanₙᵥ + H_zeeman₁₃ + H_drive_y + H_hyperfine + H_drive_freq

#####################################################################################
# Functions

"""
Calculate time evolution
Arguments:
* `hm` - Hamiltonian of the system
* `ρ` - Density operator
* `t` - Time when time evolution happens
Returns:
* New hamiltonian after time evolution
"""
function timeEvolution(hm::Union{Matrix{Complex{Float64}},Matrix{Float64}},
                       ρ::Union{Matrix{Complex{Float64}},Matrix{Float64}},
                       t::Float64)::Matrix{Complex{Float64}}
    U = exp(-im*2*π*hm*t)
    return U*ρ*adjoint(U)
end

"""
Calculate proyection of a Hamiltonian and do its trace i.e. do a measurement
Arguments:
* `hm` - Hamiltonian of the system
Returns:
* Measurement result
"""
function traceWithProy(hm::Matrix{Complex{Float64}})::Complex{Float64}
    P₀ = [1 0;0 0]
    return tr(P₀ ⊗ Matrix(1.0I,2,2)*hm)
end 

"""
 Does a Ramsey experiment and displays it.
 Input: 
   * `t` - Time for the time evolution. Usually 1/(40)/4, which means we
   use a quater of the Omega, since a full period is 1/Omega we need
   1/Omega/4, where Omega = 40
   * `ρ` - Density matrix `ρ` (electron)
"""
function ramsey(t::Float64)
    tmp = []
    for τ in 0:1.0:1000
        ρ₁ = timeEvolution(Hₓ,ρ,t)
        ρ₂ = timeEvolution(H_free,ρ₁,τ)
        ρ₃ = timeEvolution(Hₓ,ρ₂,3*t)
        append!(tmp,traceWithProy(ρ₃))
    end
    display(plot(0:1.0:1000, real.(tmp), label ="Ramsey"))
end

"""
 Auxiliary function where to times evolutions are doing. First od tehm applies
 a rotation, and second one is a free evolution.
 Input:
   * `ρ₀` - Density matrix.
   * `t` - Time for the time evolution.
   * `hm` - Hamiltonian for time evolution, with x or y drive.
"""
function πxy(ρ₀::Matrix{Complex{Float64}}, t::Float64,
             hm::Matrix{Complex{Float64}},
             τ::Float64)::Matrix{Complex{Float64}}
ρ₁ = timeEvolution(hm,ρ₀,t)
return timeEvolution(H_free,ρ₁,τ)
end

"""
 Sequence XYXYYXYX is applied the desired number of times.
 Input: 
   * `ρ₁` - Density matrix `rho` (electron)
   * `t` - Time for the time evolution.
   * `τ` - Time given for the loop of last function
 Returns: density matrix after sequence. This sequence is usefull for
 decopling electro from enviroment and nucleus (less noise)
"""
function XYXYYX(ρ₁::Matrix{Complex{Float64}}, t::Float64,
                τ::Float64)::Matrix{Complex{Float64}}
ρ₂ = πxy(ρ₁,t,Hₓ,τ)
ρ₃ = πxy(ρ₂,t,H_y,τ)
ρ₄ = πxy(ρ₃,t,Hₓ,τ)
ρ₅ = πxy(ρ₄,t,H_y,τ)
ρ₆ = πxy(ρ₅,t,H_y,τ)
ρ₇ = πxy(ρ₆,t,Hₓ,τ)
ρ₈ = πxy(ρ₇,t,H_y,τ)
ρ₉ = timeEvolution(Hₓ,ρ₈,t)
return timeEvolution(H_free,ρ₉,τ/2)
end

"""
 Does an SpinEcho experiment and displays it. Spin echo is used to
 decouple the electro sping from the nuecleus spin and the enviroment 
 (XY8 secuence), also with the dephasing it increase the decay time so we
 have more time to apply gates.
 Input: 
   * `t` - Time for the time evolution. Usually 1/(2*Omega)
   * `times` - Times the sequence XYXYYXYX will be aplied
"""
function spinEcho(t::Float64, times::Int)
    tmp = []
    for τ in 0.7:1.e-3:1.2
        ρ₁ = timeEvolution(H_y,ρ,0.5*t)
        ρ₂ = timeEvolution(H_free,ρ₁,τ/2)
        map(x -> ρ₂ = XYXYYX(ρ₂, t, τ), 1:times)
        ρ₃ = timeEvolution(H_y,ρ₂,3*0.5*t)
        append!(tmp,traceWithProy(ρ₃))
    end
    display(plot(0.7:1e-3:1.2, real.(tmp), label ="spinEcho"))
end

"""
 Same as SpinEcho but  with a fixed τ and being repeited the experiments
 the desired times.
 Input: 
   * `t` - Time for the time evolution. Usually 1/(2*Omega)
   * `times` - Times the secuence XYXYYXYX will be aplied
"""
function spinEchoN(t::Float64, times::Int)
    ω = sqrt((γ₁₃*B-Azz)^2 + Axx^2)
    τ = 1/(2*ω)
    tmp = []
    for i in 1:times
        ρ₁ = timeEvolution(H_y,ρ,0.5*t)
        ρ₂ = timeEvolution(H_free,ρ₁,τ/2)
        map(x -> ρ₂ = XYXYYX(ρ₂, t, τ), 1:i)
        ρ₃ = timeEvolution(H_y,ρ₂,3*0.5*t)
    append!(tmp,traceWithProy(ρ₃))
    end
    display(plot(1:times, real.(tmp), label ="spinEchoN"))
end

#####################################################################################
# Run

function main()
    # ramsey(1/Ω/4)
    # spinEcho(1/(2*Ω),10)
    # spinEchoN(1/(2*Ω),300)
end

main()