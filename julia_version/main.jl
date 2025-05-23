# main.jl: Example simulations for integrating the neuron model for the various injected current protocols.

import Pkg
Pkg.activate("./environment")
using Tables
using CSV
using DelimitedFiles
using Plots
using LaTeXStrings
using Measures
using DifferentialEquations

include("neuron.jl")
include("neuron_ode.jl")

# Default parameter dictionary
const parameters = Dict(
    "T" => 10.0,              # Simulation time (s, converted to ms)
    "dt" => .1,               # Output time step (ms)
    "g_NaP" => 0.0,           # Persistent sodium conductance
    "g_CAN" => 1.0,           # CAN channel conductance
    "g_Kv1" => 1.0,           # Kv1 channel conductance
    "g_Kdr" => 100.0,         # Delayed rectifier potassium conductance
    "g_NaF" => 120.0,         # Fast sodium conductance
    "g_KCa" => 0.25,          # Calcium-activated potassium conductance
    "g_L" => 0.1,             # Leak conductance
    "alpha_Ca" => 5e-6,       # Calcium influx factor
    "tau_Ca" => 250.0,        # Calcium decay time constant
    "g_Ca" => 0.01,           # Calcium channel conductance
    "K_out" => 4.0,           # External potassium concentration
    "I" => [0.0, 5.0, 20.0]   # Injected current parameters (I_min, I_step, I_max)
)

"""
integrate_segments!(neuron, filename, segments, u0, dt)

Helper function for integrating over multiple segments of the time axis for different injected current values.
The outputs are stored in a tab delimited file.

This function integrates the ODE defined by `neuron_ode!` over multiple time segments,
updating the neuron's injected current parameters (`I1`, `I2`) and simulation duration (`T`)
for each segment. The time and state variables (`t`, `V`, `Ca`, `I(t)`) are saved to a tab-separated file without headers. 
Each segment continues from the final state of the previous segment.

# Arguments
- `neuron::Neuron`                  - Neuron model, modified in-place with `I1`, `I2`, `T`.
- `filename::String`                - Output file name (e.g., `"tent_data.txt"`).
- `segments::Vector{Tuple{Tuple{Float64,Float64},Float64,Float64}}`
                                    - List of segments, each with `(t_start, t_end)`, `I1`, `I2`.
- `u0::Vector{Float64}`             - Initial state vector (7 elements: `V`, `Ca`, etc.).
- `dt::Float64`                     - Output time step (ms).

# Returns
- `Vector{Float64}`                 - Final state vector after the last segment.

# Notes
- Uses `Tsit5()` solver with `adaptive=true`, `saveat=dt`, and internal step `dt=0.1` ms.
- Output file format: Tab-separated columns `[t, V, Ca, I(t)]`
"""
function integrate_segments!(neuron::Neuron, filename::String, segments::Vector{Tuple{Tuple{Float64,Float64},Float64,Float64}}, u0::Vector{Float64}, dt::Float64)
    data = []
    u = u0
    for (tspan, I1, I2) in segments
        neuron.I1 = I1
        neuron.I2 = I2
        prob = ODEProblem(neuron_ode!, u, tspan, neuron)
        sol = solve(prob, Tsit5(), dt=0.1, adaptive=true, saveat=dt)
        for j in 1:length(sol.t)
            push!(data, [sol.t[j], sol.u[j][1], sol.u[j][2], I(neuron, sol.t[j])])
        end
        u = sol.u[end]
    end

    data_matrix = reduce(vcat, data')
    CSV.write(filename, Tables.table(data_matrix), delim='\t', writeheader=false, floatformat="%.6f")
    return u
end

"""
plot_time_series(filename, save_to_file)

Helper function for plotting the time series saved in the files generated by integrate_segments!.

# Arguments
- `filename::String`               - Filename for the data to be plotted.
- `save_to_file::Bool = false`     - Bool indicating whether you want to save the plot to file or display it.

# Returns
- Nothing

# Notes
- Uses gr() backend for plotting.
- Deletes data file after plotting.
"""
function plot_time_series(filename::String; save_to_file::Bool = false)

    data = readdlm(filename, '\t', header=false)

    t  = data[:, 1]
    v  = data[:, 2]
    ca = data[:,3]
    I  = data[:, 4]

    plot1 = plot(t, v, linewidth = 2, label = "v", legend = true, margin = 10mm)

    xlabel!(L"$t$")
    ylabel!(L"$v(t)$")
    plot!(size=(1200,800))
    plot!(titlefontsize=18, guidefontsize=18, tickfontsize=16, legendfontsize=12)

    plot2 = plot(t, ca, linewidth = 2, label = "I", legend = true, margin = 10mm)

    xlabel!(L"$t$")
    ylabel!(L"$Ca(t)$")
    plot!(size=(1200,800))
    plot!(titlefontsize=18, guidefontsize=18, tickfontsize=16, legendfontsize=12)


    plot3 = plot(t, I, linewidth = 2, label = "I", legend = true, margin = 10mm)

    xlabel!(L"$t$")
    ylabel!(L"$I(t)$")
    plot!(size=(1200,800))
    plot!(titlefontsize=18, guidefontsize=18, tickfontsize=16, legendfontsize=12)

    plot4 = plot(plot1, plot2, plot3, layout = (3,1))

    if save_to_file
        savefig(filename * ".png")
    else
        display(plot4)
        readline()
    end
    rm(filename)
end

function main()
    
    T = parameters["T"] * 1000.0
    dt = parameters["dt"]
    g_NaP = parameters["g_NaP"]
    g_CAN = parameters["g_CAN"]
    g_Kv1 = parameters["g_Kv1"]
    I_min, I_step, I_max = parameters["I"]

    if T <= 0
        error("T must be positive")
    end
    if dt <= 0
        error("dt must be positive")
    end
    if g_NaP < 0 || g_CAN < 0 || g_Kv1 < 0
        error("g_NaP, g_CAN, and g_Kv1 must be non-negative")
    end

    # Create neuron with given parameters
    neuron = Neuron(g_CAN, g_NaP, g_Kv1)
    neuron.T = T
    set_g_Kdr!(neuron, parameters["g_Kdr"])
    set_g_NaF!(neuron, parameters["g_NaF"])
    set_g_KCa!(neuron, parameters["g_KCa"])
    set_g_L!(neuron, parameters["g_L"])
    set_alpha_Ca!(neuron, parameters["alpha_Ca"])
    set_tau_Ca!(neuron, parameters["tau_Ca"])
    set_g_Ca!(neuron, parameters["g_Ca"])
    set_K_out!(neuron, parameters["K_out"])

    ramp_protocol = "step"

    if ramp_protocol == "up"
        # Constant, linear injected current ramp up protocol
        u0 = default_initial_state(neuron)
        segments = [((0.0, T), I_min, I_max)]
        integrate_segments!(neuron, "up_data.txt", segments, u0, dt)

    elseif ramp_protocol == "tent"
        # Tent injected current ramp up protocol
        u0 = default_initial_state(neuron)
        segments = [
            ((0.0, T/2.0), I_min, I_min + (I_max - I_min) * 2.0)
            ((T/2.0, T), I_min + (I_max - I_min) * 2.0, I_min)
        ]
        integrate_segments!(neuron, "tent_data.txt", segments, u0, dt)
        plot_time_series("tent_data.txt"; save_to_file = true)

    elseif  ramp_protocol == "step"
        # Step injected current ramp up protocol
        u0 = default_initial_state(neuron)
        segments = [
            ((0.0, T/6.0), I_min, I_min),
            ((T/6.0, 2.0*T/6.0), I_step, I_step),
            ((2.0*T/6.0, 3.0*T/6.0), I_max, I_max),
            ((3.0*T/6.0, 5.0*T/6.0), I_step, I_step),
            ((5.0*T/6.0, T), I_min, I_min)
        ]
        integrate_segments!(neuron, "step_data.txt", segments, u0, dt)
        plot_time_series("step_data.txt", save_to_file = true)
    end
end

main()