using Plots: length, append!, height, hascolorbar
"""
Simulation of a 15 μm laser into diamond, swhowing a rabi experiment.
"""

using Plots
using LinearAlgebra
using Statistics
using QuantumInformation
using TensorOperations
using Interpolations

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
B = 250
Δ₀ = 1
Δᵣ(r::Float64) = Δ₀ * r
# gyromagnetic ratio of the NV center’s electron spin(MHz/G)
γₙᵥ = 2.8
# gyromagnetic ratio of the nucleus carbon spin(MHz/G)
γ₁₃ = 1.07e-3

# Omega in function of `r` (position)
# Ω₀ =  16
# 
Ω(r::Float64, Ω₀::Float64) = (Ω₀*r^2)/(r^2+depth^2)^(3/2)
#Ω(r::Float64, Ω₀::Float64) = (Ω₀)/(r^2+1)


# Hyperfine
Axx = 50e-3
Ayy = 0
Azz = 10e-3

# Phase
γ = 81.3

# Omega electron (MHz)
ω = D - γₙᵥ * B

# Laser dimensions in μm
width = 15.0
depth = 1.0000000

# Distance of laser to the wire
r = 0.01
#####################################################################################
# Hamiltonians
#To get a Hamiltonian with the energy levels in the diagonal (zfs = zero field splitting)
H_zfs = D * S_z^2
# Interactions between electron and nucleous. Zeeman= when +1, -1 and 0 are
# splittig, if not -1 and +1 are in the same position.
H_zeemanₙᵥ = γₙᵥ * B * S_z
H_zeeman₁₃ = γ₁₃ * B * I_z
# To rotate in X or Y axis, Omega bc is the "line of the wave"
H_driveₓ(r::Float64, Ω₀::Float64) = Ω(r,Ω₀) * Sₓ
H_drive_y(r::Float64, Ω₀::Float64) = Ω(r,Ω₀) * S_y
H_drive_z(r::Float64, Ω₀::Float64) = Δᵣ(r)* S_z
# ω means the different between state 1 and -1 when a magnetic flied is
# aplied, we multiply it by same Hamiltonian as H_zfs to represent the
# different spin states.
H_drive_freq = - ω * S_z^2
# Hyperfine Hamiltonian. Describes hyperfine interaction between the
# NV centre’s electron spin and the nitrogen’s nuclear spin
H_hyperfine = Axx * S_z * Iₓ + Ayy * S_z * I_y + Azz * S_z * I_z


# Hamiltonians of the system

# Full electron hamiltonian with both drives
H(r::Float64, Ω₀::Float64) = H_zfs + H_zeemanₙᵥ + H_zeeman₁₃ + H_drive_z(r, Ω₀) + H_driveₓ(r, Ω₀) +H_drive_y(r, Ω₀)+ H_hyperfine + H_drive_freq
# H(r::Float64, t::Float64) = H_zfs + H_zeemanₙᵥ + H_zeeman₁₃ + H_drive_z(r) + H_drive_timeₓ(r,t) + H_hyperfine + H_drive_freq
# Electron hamiltonian with no drives
H_free = H_zfs + H_zeemanₙᵥ + H_zeeman₁₃ + H_hyperfine + H_drive_freq

# Electron hamiltonian with one drive
Hₓ(r::Float64) = H_zfs + H_zeemanₙᵥ + H_zeeman₁₃ + H_driveₓ(r) + H_hyperfine + H_drive_freq

H_y(r::Float64) = H_zfs + H_zeemanₙᵥ + H_zeeman₁₃ + H_drive_y(r) + H_hyperfine + H_drive_freq

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
                       t::Union{Float64, Int})::Matrix{Complex{Float64}}
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


function mw_element(r::Float64, t::Float64, γ::Float64)::Array{Complex{Float64}}
    Ω₀ = Ω(r,t)
    S = [0 0.5*exp(γ*im) ; 0.5*exp(-γ*im) 0]
    return Ω₀*S ⊗ Matrix(1.0I,2,2) + H_zfs + H_zeemanₙᵥ + H_drive_freq
end

function rabi_t(time::Array{Float64}, amplitudes::Array{Float64}, Δt::Float64,ampltPoints::Array{Float64})
    ampl = LinearInterpolation(time,amplitudes, extrapolation_bc=Flat())(ampltPoints)
    positions = 0.0:0.01:width
    tmp::Array{Complex{Float64}} = zeros(length(time),length(positions))
    ρ₀ = ρ
    for (i,Ωᵢ) ∈ enumerate(ampl)
        for (j,r) ∈ enumerate(positions)
            hm = mw_element(r,Ωᵢ,0.0)
            ρ₀ = timeEvolution(hm,ρ₀,Δt)
            tmp[i,j] = traceWithProy(ρ₀)*r*depth
        end
    end
    display(plot(time, real(real(sum(tmp,dims=2))), label ="rabi"))
end

"""
Does a Rabi experiment and displays it
"""
function rabi()
    time = 0.001:0.005:1
    positions = 0.0:0.01:width
    tmp = zeros(length(time),length(positions))
    for (i,t) ∈ enumerate(time)
        for (j,r) ∈ enumerate(positions)
            tmEv = timeEvolution(H(r,16.0),ρ,t)
            tmp[i,j] = traceWithProy(tmEv)*r*depth
        end
    end
    plt2 = plot(time, real(sum(tmp,dims=2)), 
                title="Rabi",
                xlabel="Time [μ]",
                ylabel="Population")
    display(plot(plt2,layout=(2,1)))
end


function getAmplitudesOrPhaseAndTimes(path::String)::Tuple{Array{Float64},Array{Float64}}
    points::Array{Float64} = []
    times::Array{Float64} = []
    open(path) do f
        while ! eof(f) 

           tmp = parse.(Float64, split(readline(f), ' ') )
           append!(points,tmp[2])
           append!(times,tmp[1])
        end
    end
    
    return points,times
end


#####################################################################################
# Run

function main()
    # rabi()
    # @show t

    # amplitudes::Array{Float64}  = 16:1:4019
    # phase::Array{Float64}  = 0:0.08992:360
# 
    # ampltPoints,times = getAmplitudesOrPhaseAndTimes("GuessPulse_pi20ns_amplitude.txt")
    # phsPoints,_ = getAmplitudesOrPhaseAndTimes("GuessPulse_pi20ns_phase.txt")
    # tmp::Array{Float64} = times[end]+0.002e-8: 0.002e-8: 8.006e-8
    # append!(times,tmp)
    # @show length(times)
    # #@show times
# 
    # append!(phsPoints,phsPoints)
    # append!(ampltPoints,ampltPoints)
    # append!(phsPoints,phsPoints)
    # append!(ampltPoints,ampltPoints)



    
    amplitudes::Array{Float64}  = 16:1:1016
    phase::Array{Float64}  = 0:0.36:360

    times::Array{Float64} = 0.0:1e-6:1e-3

    ampltPoints::Array{Float64} = zeros(length(amplitudes)) * 16
    phsPoints::Array{Float64} = zeros(length(amplitudes))

    @show length(amplitudes)
    @show length(times)
    @show length(phase)
    @show length(phsPoints)
    @show length(ampltPoints)
    rabi_t(times, amplitudes,1e-6,ampltPoints)

end

main()
