import Pkg
Pkg.activate("../julia_version/environment")
Pkg.add("Peaks")
using DelimitedFiles
using Plots
using LaTeXStrings
using Measures
using Peaks
using DataFrames
using CSV

data = readdlm("./current_step_data.txt", '\t', header=false)

t = data[:, 1]
v = data[:, 2]
ca = data[:,3]
I = data[:, 4]

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

savefig("plot.png");