#include <vector>
#include <cmath>
#include <stdexcept>

typedef std::vector<double> state_type;

class Neuron {
private:
    // Conductances
    double gNaP, gCAN, gKv1;
    static inline double gKCa = 0.0, gNaF = 120.0, gKdr = 100.0, g_L = 0.1;
    static inline double gCa = 0.01, tauCa = 250.0, alphaCa = 5e-6;

    // Reversal potentials and constants
    static inline double Kin = 140.0, Kout = 4.0;
    static inline double EK = 26.54 * std::log(Kout / Kin);
    static constexpr double E_leak = -80.0;
    const double C = 1.0, ENa = 55.0, ECAN = 0.0, ECa = 80.0;
    const double K_CAN = 0.74e-3, Kd = 0.2e-3;

    // Update EK when Kin or Kout changes
    static void update_EK() {
        if (Kout <= 0 || Kin <= 0) {
            throw std::invalid_argument("Kout and Kin must be positive");
        }
        EK = 26.54 * std::log(Kout / Kin);
    }

public:
    double I1, I2, T;

    // Getters
    static double get_gKCa() { return gKCa; }
    static double get_gNaF() { return gNaF; }
    static double get_gKdr() { return gKdr; }
    static double get_g_L() { return g_L; }
    static double get_gCa() { return gCa; }
    static double get_tauCa() { return tauCa; }
    static double get_alphaCa() { return alphaCa; }
    static double get_Kin() { return Kin; }
    static double get_Kout() { return Kout; }
    static double get_EK() { return EK; }
    static const double get_E_leak() { return E_leak; }

    // Setters with validation
    static void set_gKCa(double value) {
        if (value < 0) throw std::invalid_argument("gKCa must be non-negative");
        gKCa = value;
    }
    static void set_gNaF(double value) {
        if (value < 0) throw std::invalid_argument("gNaF must be non-negative");
        gNaF = value;
    }
    static void set_gKdr(double value) {
        if (value < 0) throw std::invalid_argument("gKdr must be non-negative");
        gKdr = value;
    }
    static void set_g_L(double value) {
        if (value < 0) throw std::invalid_argument("g_L must be non-negative");
        g_L = value;
    }
    static void set_gCa(double value) {
        if (value < 0) throw std::invalid_argument("gCa must be non-negative");
        gCa = value;
    }
    static void set_tauCa(double value) {
        if (value <= 0) throw std::invalid_argument("tauCa must be positive");
        tauCa = value;
    }
    static void set_alphaCa(double value) {
        if (value < 0) throw std::invalid_argument("alphaCa must be non-negative");
        alphaCa = value;
    }
    static void set_Kin(double value) {
        if (value <= 0) throw std::invalid_argument("Kin must be positive");
        Kin = value;
        update_EK();
    }
    static void set_Kout(double value) {
        if (value <= 0) throw std::invalid_argument("Kout must be positive");
        Kout = value;
        update_EK();
    }

    // Default initial state
    static state_type default_initial_state() {
        state_type x(7, 0.0);
        x[0] = get_E_leak();
        x[2] = 1e-6;
        x[6] = 1.0;
        return x;
    }

    Neuron(double gCAN, double gNaP, double gKv1) : gCAN(gCAN), gNaP(gNaP), gKv1(gKv1) {
        if (gCAN < 0 || gNaP < 0 || gKv1 < 0) {
            throw std::invalid_argument("Conductances must be non-negative");
        }
    }

    inline double I(double t) const { return I1 + (I2 - I1) * t / T; }

    void operator()(const state_type& x, state_type& dxdt, double t) const {
        // State variables
        const double& v = x[0];    // Membrane potential
        const double& Cain = x[1]; // Calcium concentration
        const double& mc = x[2];   // Calcium channel activation
        const double& hc = x[3];   // Calcium channel inactivation
        const double& h = x[4];    // Sodium channel inactivation
        const double& n = x[5];    // Potassium channel activation
        const double& z = x[6];    // Kv1 channel activation

        // Calcium currents
        double mcinf = 1.0 / (1.0 + std::exp(-(v + 27.5) / 5.7));
        double hcinf = 1.0 / (1.0 + std::exp((v + 52.4) / 5.2));
        double ICa = gCa * mc * hc;
        double IKCa = gKCa * Cain / (Cain + Kd);
        double ICAN = gCAN * Cain / (Cain + K_CAN);

        // Sodium and potassium currents
        double minf = 1.0 / (1.0 + std::exp(-(v + 35.0) / 7.8));
        double mpinf = 1.0 / (1.0 + std::exp(-(v + 53.0) / 3.0)); // Brocard et al. 2013
        double winf = 1.0 / (1.0 + std::exp(-(v + 46.0) / 6.9));  // Bos et al. 2021
        double INaF = gNaF * minf * minf * minf * h;
        double INaP = gNaP * mpinf;
        double IKdr = gKdr * n * n * n * n;
        double IKv1 = gKv1 * winf * z;

        // Differential equations
        dxdt[0] = (-IKCa * (v - EK) - ICAN * (v - ECAN) - ICa * (v - ECa) - g_L * (v - E_leak) + I(t)
                  - (INaF + INaP) * (v - ENa) - (IKdr + IKv1) * (v - EK)) / C;
        dxdt[1] = -alphaCa * ICa * (v - ECa) - Cain / tauCa;
        dxdt[2] = (mcinf - mc) / 0.5;
        dxdt[3] = (hcinf - hc) / 18.0;

        // Sodium and potassium channel dynamics
        double hinf = 1.0 / (1.0 + std::exp((v + 55.0) / 7.0));
        double htau = 30.0 / (std::exp((v + 50.0) / 15.0) + std::exp(-(v + 50.0) / 16.0));
        dxdt[4] = (hinf - h) / htau;

        double ninf = 1.0 / (1.0 + std::exp(-(v + 28.0) / 15.0));
        double ntau = 7.0 / (std::exp((v + 40.0) / 40.0) + std::exp(-(v + 40.0) / 50.0));
        dxdt[5] = (ninf - n) / ntau;

        double zinf = 1.0 / (1.0 + std::exp((v + 54.0) / 7.1)); // Bos et al. 2018 + 2021
        double hktau = (3.737 / (0.00015 * std::exp(-(v + 13.0) / 15.0) + 0.06 / (1.0 + std::exp(-(v + 68.0) / 12.0)))) * 20.0;
        dxdt[6] = (zinf - z) / hktau;
    }
};