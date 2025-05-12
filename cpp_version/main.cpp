#include <iostream>
#include <fstream>
#include <vector>
#include <stdexcept>
#ifdef _OPENMP
    #include <omp.h>
#endif
#include <boost/numeric/odeint.hpp>
#include "neuron.h"

using namespace boost::numeric::odeint;

struct OutputObserver {
    Neuron& neuron;
    std::ostream& out;
    double dt, t_last = 0.0, v_prev = 0.0, v_prev_prev = 0.0;

    OutputObserver(Neuron& n, std::ostream& os, double step = 1.0) : neuron(n), out(os), dt(step) {}

    void operator()(const state_type& x, double t) {
        bool is_extrema = (v_prev > v_prev_prev && v_prev > x[0]) || (v_prev < v_prev_prev && v_prev < x[0]);
        if (t >= t_last + dt || is_extrema) {
            out << t << '\t' << x[0] << '\t' << x[1] << '\t' << neuron.I(t) << '\n';
            t_last = t;
        }
        v_prev_prev = v_prev;
        v_prev = x[0];
    }
};

// Parse command-line arguments and update parameters
void parse_arguments(int argc, char** argv, double& gCAN, double& gNaP, double& gKv1,
                    double& Imin, double& Imax, double& Ib, double& T, double& dt) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-T") {
            if (++i >= argc) throw std::invalid_argument("Missing value for -T");
            T = std::atof(argv[i]) * 1000.0;
            if (T <= 0) throw std::invalid_argument("T must be positive");
        }
        else if (arg == "-dt") {
            if (++i >= argc) throw std::invalid_argument("Missing value for -dt");
            dt = std::atof(argv[i]);
            if (dt <= 0) throw std::invalid_argument("dt must be positive");
        }
        else if (arg == "-nap") {
            if (++i >= argc) throw std::invalid_argument("Missing value for -nap");
            gNaP = std::atof(argv[i]);
            if (gNaP < 0) throw std::invalid_argument("gNaP must be non-negative");
        }
        else if (arg == "-can") {
            if (++i >= argc) throw std::invalid_argument("Missing value for -can");
            gCAN = std::atof(argv[i]);
            if (gCAN < 0) throw std::invalid_argument("gCAN must be non-negative");
        }
        else if (arg == "-ka") {
            if (++i >= argc) throw std::invalid_argument("Missing value for -ka");
            gKv1 = std::atof(argv[i]);
            if (gKv1 < 0) throw std::invalid_argument("gKv1 must be non-negative");
        }
        else if (arg == "-kdr") {
            if (++i >= argc) throw std::invalid_argument("Missing value for -kdr");
            Neuron::set_gKdr(std::atof(argv[i]));
        }
        else if (arg == "-naf") {
            if (++i >= argc) throw std::invalid_argument("Missing value for -naf");
            Neuron::set_gNaF(std::atof(argv[i]));
        }
        else if (arg == "-kca") {
            if (++i >= argc) throw std::invalid_argument("Missing value for -kca");
            Neuron::set_gKCa(std::atof(argv[i]));
        }
        else if (arg == "-gl") {
            if (++i >= argc) throw std::invalid_argument("Missing value for -gl");
            Neuron::set_g_L(std::atof(argv[i]));
        }
        else if (arg == "-a") {
            if (++i >= argc) throw std::invalid_argument("Missing value for -a");
            Neuron::set_alphaCa(std::atof(argv[i]));
        }
        else if (arg == "-t") {
            if (++i >= argc) throw std::invalid_argument("Missing value for -t");
            Neuron::set_tauCa(std::atof(argv[i]));
        }
        else if (arg == "-ca") {
            if (++i >= argc) throw std::invalid_argument("Missing value for -ca");
            Neuron::set_gCa(std::atof(argv[i]));
        }
        else if (arg == "-kout") {
            if (++i >= argc) throw std::invalid_argument("Missing value for -kout");
            Neuron::set_Kout(std::atof(argv[i]));
        }
        else if (arg == "-I") {
            if (i + 3 >= argc) throw std::invalid_argument("Missing values for -I (requires Imin, Ib, Imax)");
            Imin = std::atof(argv[++i]);
            Ib = std::atof(argv[++i]);
            Imax = std::atof(argv[++i]);
        }
        else {
            throw std::invalid_argument("Unknown option: " + arg);
        }
    }
}

int main(int argc, char** argv) {
    // Default parameters
    double gCAN = 0.0, gNaP = 0.0, gKv1 = 0.0;
    double Imin = 0.0, Imax = 0.0, Ib = 0.0, T = 10000.0, dt = 1.0;

    // Parse command-line arguments
    try {
        parse_arguments(argc, argv, gCAN, gNaP, gKv1, Imin, Imax, Ib, T, dt);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    const int dim = 7;

    #pragma omp parallel sections
    {
        //Constant, linear injected current ramp up protocol
        #pragma omp section
        {
            Neuron neuron(gCAN, gNaP, gKv1);
            neuron.T = T;
            neuron.I1 = Imin;
            neuron.I2 = Imax;
            state_type x = Neuron::default_initial_state();
            std::ofstream os("dat0");
            integrate(neuron, x, 0.0, T, 0.1, OutputObserver(neuron, os, dt));
        }
        //Tent injected current ramp up protocol
        #pragma omp section
        {
            Neuron neuron(gCAN, gNaP, gKv1);
            neuron.T = T;
            state_type x = Neuron::default_initial_state();
            std::ofstream os("dat1");
            neuron.I1 = Imin;
            neuron.I2 = Imax * 2.0;
            integrate(neuron, x, 0.0, T / 2.0, 0.1, OutputObserver(neuron, os, dt));
            neuron.I1 = Imax * 2.0;
            neuron.I2 = Imin;
            integrate(neuron, x, T / 2.0, T, 0.1, OutputObserver(neuron, os, dt));
        }
        #pragma omp section
        {
            Neuron neuron(gCAN, gNaP, gKv1);
            neuron.T = T;
            state_type x = Neuron::default_initial_state();
            x[0] = -35.0;
            x[2] = 1e-5;
            std::ofstream ns("/dev/null");
            neuron.I1 = Imax;
            neuron.I2 = Imax;
            integrate(neuron, x, 0.0, 1000.0, 0.1, OutputObserver(neuron, ns, dt));
            std::ofstream os("dat2");
            neuron.I1 = Imax;
            neuron.I2 = Imin;
            integrate(neuron, x, 0.0, T, 0.1, OutputObserver(neuron, os, dt));
        }
        //Step injected current ramp up protocol
        #pragma omp section
        {
            Neuron neuron(gCAN, gNaP, gKv1);
            neuron.T = T;
            state_type x = Neuron::default_initial_state();
            std::ofstream os("dat3");
            neuron.I1 = neuron.I2 = Imin;
            integrate(neuron, x, 0.0, T / 6.0, 0.1, OutputObserver(neuron, os, dt));
            neuron.I1 = neuron.I2 = Ib;
            integrate(neuron, x, T / 6.0, 2.0 * T / 6.0, 0.1, OutputObserver(neuron, os, dt));
            neuron.I1 = neuron.I2 = Imax;
            integrate(neuron, x, 2.0 * T / 6.0, 3.0 * T / 6.0, 0.1, OutputObserver(neuron, os, dt));
            neuron.I1 = neuron.I2 = Ib;
            integrate(neuron, x, 3.0 * T / 6.0, 5.0 * T / 6.0, 0.1, OutputObserver(neuron, os, dt));
            neuron.I1 = neuron.I2 = Imin;
            integrate(neuron, x, 5.0 * T / 6.0, T, 0.1, OutputObserver(neuron, os, dt));
        }
    }

    return 0;
}