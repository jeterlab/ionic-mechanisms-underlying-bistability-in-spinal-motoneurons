# neuron_ode.jl: ODE system for neuron simulations

include("neuron.jl")

"""
neuron_ode!(dx, x, neuron, t)

Function for integrating the neuron ODEs. 

The Julia integrators in the DifferentialEquations package expect a function that takes the state of the state variables,
a vector for their derivatives, the current value of the time variable t, and additional parameters used in your system.


# Arguments
- `dx::Vector{Float64}`    - Vector with the values of the derivatives.
- `x::Vector{Float64}`     - Vector with the current values of the state variables.
- `neuron::Neuron`         - Neuron model, modified in-place with `I1`, `I2`, `T`.
- `t::Float64`             - Value of t for computing the derivatives of the state variables.

# Returns
- Nothing

# Notes
- This function operates on the inputs in place, hence the ! notation.
"""
function neuron_ode!(dx::Vector{Float64}, x::Vector{Float64}, neuron::Neuron, t::Float64)
    
    v = x[1]      # Membrane potential
    Ca = x[2]     # Calcium concentration
    m_CaL = x[3]  # High voltage Calcium channel activation
    h_CaL = x[4]  # High voltage Calcium channel inactivation
    h_NaF = x[5]  # Fast Sodium channel inactivation
    m_Kdr = x[6]  # Potassium rectifier channel activation
    h_Kv1 = x[7]  # Kv1.2 channel inactivation

    # Calcium currents
    m_CaLInf = 1.0 / (1.0 + exp(-(v + 27.5) / 5.7))
    h_CaLInf = 1.0 / (1.0 + exp((v + 52.4) / 5.2))
    I_Ca = neuron.g_Ca * m_CaL * h_CaL
    I_KCa = neuron.g_KCa * Ca / (Ca + neuron.Kd)
    I_CAN = neuron.g_CAN * Ca / (Ca + neuron.K_CAN)

    # Sodium and potassium currents
    m_NaFInf = 1.0 / (1.0 + exp(-(v + 35.0) / 7.8))
    m_NaPInf = 1.0 / (1.0 + exp(-(v + 53.0) / 3.0))
    m_Kv1Inf = 1.0 / (1.0 + exp(-(v + 46.0) / 6.9))
    I_NaF = neuron.g_NaF * m_NaFInf^3 * h_NaF
    I_NaP = neuron.g_NaP * m_NaPInf
    I_Kdr = neuron.g_Kdr * m_Kdr^4
    I_Kv1 = neuron.g_Kv1 * m_Kv1Inf * h_Kv1

    # Differential equations
    dx[1] = (-I_KCa * (v - neuron.E_K) - I_CAN * (v - neuron.E_CAN) - I_Ca * (v - neuron.E_Ca) - 
             neuron.g_L * (v - neuron.E_leak) + I(neuron, t) - 
             (I_NaF + I_NaP) * (v - neuron.E_Na) - (I_Kdr + I_Kv1) * (v - neuron.E_K)) / neuron.C
    dx[2] = -neuron.alpha_Ca * I_Ca * (v - neuron.E_Ca) - Ca / neuron.tau_Ca
    dx[3] = (m_CaLInf - m_CaL) / 0.5
    dx[4] = (h_CaLInf - h_CaL) / 18.0

    # Fast sodium channel dynamics
    h_NaFInf = 1.0 / (1.0 + exp((v + 55.0) / 7.0))
    h_NaFTau = 30.0 / (exp((v + 50.0) / 15.0) + exp(-(v + 50.0) / 16.0))
    dx[5] = (h_NaFInf - h_NaF) / h_NaFTau

    # Potassium rectifier channel dynamics
    m_KdrInf = 1.0 / (1.0 + exp(-(v + 28.0) / 15.0))
    m_KdrTau = 7.0 / (exp((v + 40.0) / 40.0) + exp(-(v + 40.0) / 50.0))
    dx[6] = (m_KdrInf - m_Kdr) / m_KdrTau

    # Kv1.2 channel dynamics
    h_Kv1Inf = 1.0 / (1.0 + exp((v + 54.0) / 7.1))
    h_Kv1Tau = (3.737 / (0.00015 * exp(-(v + 13.0) / 15.0) + 0.06 / (1.0 + exp(-(v + 68.0) / 12.0)))) * 20.0
    dx[7] = (h_Kv1Inf - h_Kv1) / h_Kv1Tau
end