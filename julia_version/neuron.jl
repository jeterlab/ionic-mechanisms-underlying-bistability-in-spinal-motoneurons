# neuron.jl: Neuron model struct, initial state, and injected current function.

mutable struct Neuron
    
    # Conductances (instance-specific)
    g_NaP::Float64
    g_CAN::Float64
    g_Kv1::Float64
    
    # Conductances (shared, static in C++)
    g_KCa::Float64
    g_NaF::Float64
    g_Kdr::Float64
    g_L::Float64
    g_Ca::Float64
    tau_Ca::Float64
    alpha_Ca::Float64
    
    # Reversal potentials and constants
    K_in::Float64
    K_out::Float64
    E_K::Float64
    E_leak::Float64
    C::Float64
    E_Na::Float64
    E_CAN::Float64
    E_Ca::Float64
    K_CAN::Float64
    Kd::Float64

    # Simulation parameters
    I1::Float64
    I2::Float64
    T::Float64

    # Constructor for the neuron struct.
    function Neuron(g_CAN, g_NaP, g_Kv1; 
                    g_KCa=0.0, g_NaF=120.0, g_Kdr=100.0, g_L=0.1,
                    g_Ca=0.01, tau_Ca=250.0, alpha_Ca=5e-6,
                    K_in=140.0, K_out=4.0)
        if g_CAN < 0 || g_NaP < 0 || g_Kv1 < 0 || g_KCa < 0 || g_NaF < 0 || g_Kdr < 0 || g_L < 0 || g_Ca < 0 || alpha_Ca < 0
            throw(ArgumentError("Conductances and alpha_Ca must be non-negative"))
        end
        if tau_Ca <= 0
            throw(ArgumentError("tau_Ca must be positive"))
        end
        if K_in <= 0 || K_out <= 0
            throw(ArgumentError("K_in and K_out must be positive"))
        end
        E_K = 26.54 * log(K_out / K_in)
        new(g_NaP, g_CAN, g_Kv1, g_KCa, g_NaF, g_Kdr, g_L, g_Ca, tau_Ca, alpha_Ca, 
            K_in, K_out, E_K, -80.0, 1.0, 55.0, 0.0, 80.0, 0.74e-3, 0.2e-3, 
            0.0, 0.0, 0.0)
    end
end

# Setters for the variables of the neuron struct that are manipulatable.
function set_g_KCa!(neuron::Neuron, value::Float64)
    if value < 0
        throw(ArgumentError("g_KCa must be non-negative"))
    end
    neuron.g_KCa = value
end

function set_g_NaF!(neuron::Neuron, value::Float64)
    if value < 0
        throw(ArgumentError("g_NaF must be non-negative"))
    end
    neuron.g_NaF = value
end

function set_g_Kdr!(neuron::Neuron, value::Float64)
    if value < 0
        throw(ArgumentError("g_Kdr must be non-negative"))
    end
    neuron.g_Kdr = value
end

function set_g_L!(neuron::Neuron, value::Float64)
    if value < 0
        throw(ArgumentError("g_L must be non-negative"))
    end
    neuron.g_L = value
end

function set_g_Ca!(neuron::Neuron, value::Float64)
    if value < 0
        throw(ArgumentError("g_Ca must be non-negative"))
    end
    neuron.g_Ca = value
end

function set_tau_Ca!(neuron::Neuron, value::Float64)
    if value <= 0
        throw(ArgumentError("tau_Ca must be positive"))
    end
    neuron.tau_Ca = value
end

function set_alpha_Ca!(neuron::Neuron, value::Float64)
    if value < 0
        throw(ArgumentError("alpha_Ca must be non-negative"))
    end
    neuron.alpha_Ca = value
end

function set_K_in!(neuron::Neuron, value::Float64)
    if value <= 0
        throw(ArgumentError("K_in must be positive"))
    end
    neuron.K_in = value
    neuron.E_K = 26.54 * log(neuron.K_out / neuron.K_in)
end

function set_K_out!(neuron::Neuron, value::Float64)
    if value <= 0
        throw(ArgumentError("K_out must be positive"))
    end
    neuron.K_out = value
    neuron.E_K = 26.54 * log(neuron.K_out / neuron.K_in)
end

"""
I(neuron, t)

Current function that will return the current at a given time step given the two injected current values for the neuron.

# Arguments
- `neuron::Neuron`    - Neuron model with set parameters.

- `t::Float64`        - Time point to get the value of the injected current

# Returns
- `Float64`   - Value of the injected current for the given neuron at the given time.

# Notes
- This function is able to be used for all injected current protocols described in the manuscript
  by varying the values of I1 and I2.

"""
function I(neuron::Neuron, t::Float64)
    return neuron.I1 + (neuron.I2 - neuron.I1) * t / neuron.T
end

"""
default_initial_state(neuron)

Helper function creating the initial conditions array for the state variables.

# Arguments
- `neuron::Neuron`    - Neuron model with set parameters.

# Returns
- `Vector{Float64}`   - Initial conditions for state variables

"""
function default_initial_state(neuron::Neuron)
    x = zeros(7)
    x[1] = neuron.E_leak
    x[3] = 1e-6  # m_CaL
    x[7] = 1.0   # h_Kv1
    return x
end