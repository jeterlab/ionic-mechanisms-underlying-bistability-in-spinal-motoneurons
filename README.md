# Ionic Mechanisms Underlying Bistability in Spinal Motoneurons

This repository is the companion code for the computational model described in the paper *"Ionic Mechanisms Underlying Bistability in Spinal Motoneurons"*. It provides implementations of the model in both C++ and Julia.

## Repository Structure

The repository is organized into two main directories:

- **`cpp_version/`**: Houses the C++ implementation of the computational model, including a script to compile, run, and plot results.
- **`julia_version/`**: Contains the Julia implementation of the model, including a main script and the environment setup.

### Directory Structure Diagram

```
.
├── cpp_version/
│   ├── main.cpp
│   ├── neuron.h
│   ├── run_code.sh
│   └── plot.jl
├── julia_version/
│   ├── environment/
│   │   ├── Project.toml
│   │   └── Manifest.toml
│   ├── main.jl
│   ├── neuron_ode.jl
│   └── neuron.jl
```

- **`cpp_version/`**: Contains C++ source files (`main.cpp`, `neuron.h`), a Bash script (`run_code.sh`) to compile and run the model, and a Julia script (`plot.jl`) for plotting.
- **`julia_version/`**: Includes the Julia main script (`main.jl`), model implementation (`neuron.jl` and `neuron_ode.jl`), and the `environment/` subdirectory with Julia project files (`Project.toml`, `Manifest.toml`).


## Prerequisites

- **C++ Implementation**:
  - A C++ compiler (e.g., `g++`) 
    - OpenMP is recommended for running simulations in parallel
  - Visualization software (currently the example uses Julia to plot the results, but gnuplot and Python also work well.)
- **Julia Implementation**:
  - Julia (version 1.6 or later recommended)
  - Julia dependencies listed in `environment/Project.toml`

## Running the C++ Implementation

1. Navigate to the `cpp_version/` directory:
   ```bash
   cd cpp_version
   ```
2. Run the provided Bash script to compile the C++ code, execute the model with example command-line flags, and generate plots using Julia:
   ```bash
   bash run_code.sh
   ```
3. The script will output results and plots as specified in the example configuration.

## Running the Julia Implementation

1. Ensure the Julia environment is set up by activating the project environment in the `environment/` directory:
   ```
   cd environment
   julia --project=. -e "using Pkg; Pkg.instantiate()"
   ```
2. Navigate to the `julia_version/` directory:
   ```
   cd ../julia_version
   ```
3. Run the Julia model:
   ```
   julia main.jl
   ```
4. The script will integrate the neuron model with the parameter set and injected current ramping protocol outlined in `main.jl` and produce a plot.

## Notes

- Ensure all dependencies in the `environment/` directory are installed before running the Julia code.
- The C++ implementation requires a compatible compiler and Julia for plotting. Modify `run_code.sh` if you need to adjust command-line flags, compiler settings, or the visualization software.
- For detailed model descriptions, refer to the associated paper.

## Contact

For questions or issues, please open an issue on this repository or contact [Russell Jeter](math.gsu.edu/rjeter) at [rjeter@gsu.edu](mailto:rjeter@gsu.edu).