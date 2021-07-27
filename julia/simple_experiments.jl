""" 
Three level spin system. Simulation for the following experiments in a simple way 
i.e. not taking into account carbon13 nucleus:
    * Observe Spins under magnetic field
    * ODMR experiment
    * Rabi experiment 
"""

using Plots: length
using Plots
using LinearAlgebra
using QuantumInformation
using TensorOperations

#####################################################################################
# operators

# Electron operators
Sₓ = 1/2 * [0 1 0; 1 0 1; 0 1 0]
S_z = [1 0 0; 0 0 0; 0 0 -1]

# Density matrix
ρ = [0.0 0.0 0.0; 0.0 1.0 0.0;0.0 0.0 0.0]

#####################################################################################
# Parameters

# Energy between spin state 0 and +-1 (zfs in ground state)(MHz)
const D = 2870
# Magnetic field (G)
B = 1500
# gyromagnetic ratio of the NV center’s electron spin(MHz/G)
const γₙᵥ = 2.8
const Ω = 40;

# Hyperfine
const A = 2.16


#####################################################################################
# Hamiltonians

#To get a Hamiltonian with the energy levels in the diagonal (zfs = zero field splitting)
H_zfs = D * S_z^2
# Interactions between electron and nucleous. Zeeman= when +1, -1 and 0 are
# splittig, if not -1 and +1 are in the same position.
H_zeemanₙᵥ(;B=B) = γₙᵥ * B * S_z;
# ω means the different between state 1 and -1 when a magnetic flied is
# aplied, we multiply it by same Hamiltonian as H_zfs to represent the
# different spin states.
H_drive_freq(ω) = - ω * S_z^2;
# To rotate in X axis, Ω bc is the "line of the wave"
H_driveₓ = Ω*Sₓ

# Hamiltonians of the system

# Drive Hamiltonian
H_drive(ω) =  H_driveₓ + H_drive_freq(ω)
# Electron hamiltonian with no drives
H_free(;B=B) = H_zfs + H_zeemanₙᵥ(;B)

#####################################################################################
# Functions

"""
Displays the energy levels in function of the magnetic field
"""
function energyB()
    dg::Array{Float64} = zeros(length(1:B),3)
    map(b -> dg[b,:] = deepcopy(eigvals(H_free(;B=b))), 1:B)
    display(plot(1:B,dg, 
        title="Energy levels",
            xlabel="Magnetic field", 
            ylabel="Energy",
            labels=[1 -1 0]))
end

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
                       ρ::Matrix{Float64},t::Float64)::Matrix{Complex{Float64}}
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
function traceWithProy(hm::Matrix{Complex{Float64}})::Float64
    P₀ = zeros(3,3)
    P₀[5] = 1
    return tr(P₀*hm)
end 

"""
Does a Rabi experiment and displays it
Returns
* Time when it's the first lowest point in the graphic.
"""
function rabi()::Float64
    H = H_drive(0.0)
    time = 0.001:0.001:0.05
    tmp = zeros(length(time))
    t = 0
    check = true
    for (j,i) in enumerate(time)
        tmEv = timeEvolution(H,ρ,i)
        tmp[j] = traceWithProy(tmEv)
        if j != 1 && check && tmp[j] >  tmp[j-1]
            t = i
            check = false
        end 
    end
    display(plot(0.001:0.001:0.05, tmp, 
            title="Rabi",
            xlabel="Time [μ]",
            ylabel="Population"))
    return t
end

"""
Does an ODMR experiment and displays  it
Input:
* `t` - Time where the first lower point in rabi
"""
function odmr(t::Float64)
    @show t
    energies = D-300:D+300
    tmp = zeros(length(energies))
    for (j,i) in enumerate(energies)
        ρ₁ = timeEvolution(H_free(;B=1)+H_drive(i),ρ,t)
        tmp[j] = traceWithProy(ρ₁)
    end
    display(plot(1:601, tmp, label ="ODRM"))
end
#####################################################################################
# Run

function main()
    # getms()
    # t = rabi()
    # odmr(t)
end

main()