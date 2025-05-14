#include <vector>
#include <cmath>
#include <stdexcept>

typedef std::vector<double> state_type;

class Neuron {
private:
    // Conductances
    double g_NaP, g_CAN, g_Kv1;
    static inline double g_KCa = 0.0, g_NaF = 120.0, g_Kdr = 100.0, g_L = 0.1;
    static inline double g_Ca = 0.01, tau_Ca = 250.0, alpha_Ca = 5e-6;

    // Reversal potentials and constants
    static inline double K_in = 140.0, K_out = 4.0;
    static inline double E_K = 26.54 * std::log(K_out / K_in);
    static constexpr double E_leak = -80.0;
    const double C = 1.0, E_Na = 55.0, E_CAN = 0.0, E_Ca = 80.0;
    const double K_CAN = 0.74e-3, Kd = 0.2e-3;

    // Update E_K when K_in or K_out changes
    static void update_EK() {
        if (K_out <= 0 || K_in <= 0) {
            throw std::invalid_argument("K_out and K_in must be positive");
        }
        E_K = 26.54 * std::log(K_out / K_in);
    }

public:
    double I1, I2, T;

    // Getters
    static double get_g_KCa() { return g_KCa; }
    static double get_g_NaF() { return g_NaF; }
    static double get_g_Kdr() { return g_Kdr; }
    static double get_g_L() { return g_L; }
    static double get_g_Ca() { return g_Ca; }
    static double get_tau_Ca() { return tau_Ca; }
    static double get_alpha_Ca() { return alpha_Ca; }
    static double get_K_in() { return K_in; }
    static double get_K_out() { return K_out; }
    static double get_E_K() { return E_K; }
    static const double get_E_leak() { return E_leak; }

    // Setters for private variables that are updated through command line arguments
    static void set_g_KCa(double value) {
        if (value < 0) throw std::invalid_argument("g_KCa must be non-negative");
        g_KCa = value;
    }
    static void set_g_NaF(double value) {
        if (value < 0) throw std::invalid_argument("g_NaF must be non-negative");
        g_NaF = value;
    }
    static void set_g_Kdr(double value) {
        if (value < 0) throw std::invalid_argument("g_Kdr must be non-negative");
        g_Kdr = value;
    }
    static void set_g_L(double value) {
        if (value < 0) throw std::invalid_argument("g_L must be non-negative");
        g_L = value;
    }
    static void set_g_Ca(double value) {
        if (value < 0) throw std::invalid_argument("g_Ca must be non-negative");
        g_Ca = value;
    }
    static void set_tau_Ca(double value) {
        if (value <= 0) throw std::invalid_argument("tau_Ca must be positive");
        tau_Ca = value;
    }
    static void set_alpha_Ca(double value) {
        if (value < 0) throw std::invalid_argument("alpha_Ca must be non-negative");
        alpha_Ca = value;
    }
    static void set_K_in(double value) {
        if (value <= 0) throw std::invalid_argument("K_in must be positive");
        K_in = value;
        update_EK();
    }
    static void set_K_out(double value) {
        if (value <= 0) throw std::invalid_argument("K_out must be positive");
        K_out = value;
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

    Neuron(double g_CAN, double g_NaP, double g_Kv1) : g_CAN(g_CAN), g_NaP(g_NaP), g_Kv1(g_Kv1) {
        if (g_CAN < 0 || g_NaP < 0 || g_Kv1 < 0) {
            throw std::invalid_argument("Conductances must be non-negative");
        }
    }

    inline double I(double t) const { return I1 + (I2 - I1) * t / T; }

    void operator()(const state_type& x, state_type& dxdt, double t) const {
        // State variables
        const double& v = x[0];    // Membrane potential
        const double& Ca = x[1]; // Calcium concentration
        const double& m_CaL = x[2];   // High voltage Calcium channel activation
        const double& h_CaL = x[3];   // High voltage Calcium channel inactivation
        const double& h_NaF = x[4];    // Fast Sodium channel inactivation
        const double& m_Kdr = x[5];    // Potassium rectifier channel activation
        const double& h_Kv1 = x[6];    // Kv1.2 channel inactivation

        // Calcium currents
        double m_CaLInf = 1.0 / (1.0 + std::exp(-(v + 27.5) / 5.7));
        double h_CaLInf = 1.0 / (1.0 + std::exp((v + 52.4) / 5.2));
        double ICa = g_Ca * m_CaL * h_CaL;
        double IKCa = g_KCa * Ca / (Ca + Kd);
        double ICAN = g_CAN * Ca / (Ca + K_CAN);

        // Sodium and potassium currents
        double m_NaFInf = 1.0 / (1.0 + std::exp(-(v + 35.0) / 7.8));
        double m_NaPInf = 1.0 / (1.0 + std::exp(-(v + 53.0) / 3.0)); // Brocard et al. 2013
        double m_Kv1Inf = 1.0 / (1.0 + std::exp(-(v + 46.0) / 6.9));  // Bos et al. 2021
        double I_NaF = g_NaF * m_NaFInf * m_NaFInf * m_NaFInf * h_NaF;
        double I_NaP = g_NaP * m_NaPInf;
        double I_Kdr = g_Kdr * m_Kdr * m_Kdr * m_Kdr * m_Kdr;
        double I_Kv1 = g_Kv1 * m_Kv1Inf * h_Kv1;

        // Membrane potential dynamics
        dxdt[0] = (-IKCa * (v - E_K) - ICAN * (v - E_CAN) - ICa * (v - E_Ca) - g_L * (v - E_leak) + I(t)
                  - (I_NaF + I_NaP) * (v - E_Na) - (I_Kdr + I_Kv1) * (v - E_K)) / C;
        //Calcium concentration dynamics
        dxdt[1] = -alpha_Ca * ICa * (v - E_Ca) - Ca / tau_Ca;
        
        // High voltage calcium channel dynamics
        dxdt[2] = (m_CaLInf - m_CaL) / 0.5;
        dxdt[3] = (h_CaLInf - h_CaL) / 18.0;

        // Fast sodium channel dynamics
        double h_NaFInf = 1.0 / (1.0 + std::exp((v + 55.0) / 7.0));
        double h_NaFTau = 30.0 / (std::exp((v + 50.0) / 15.0) + std::exp(-(v + 50.0) / 16.0));
        dxdt[4] = (h_NaFInf - h_NaF) / h_NaFTau;

        // Potassium rectifier channel dynamics
        double m_KdrInf = 1.0 / (1.0 + std::exp(-(v + 28.0) / 15.0));
        double m_KdrTau = 7.0 / (std::exp((v + 40.0) / 40.0) + std::exp(-(v + 40.0) / 50.0));
        dxdt[5] = (m_KdrInf - m_Kdr) / m_KdrTau;
        
        // Kv1.2 channel dynamics
        double h_Kv1Inf = 1.0 / (1.0 + std::exp((v + 54.0) / 7.1)); // Bos et al. 2018 + 2021
        double h_Kv1Tau = (3.737 / (0.00015 * std::exp(-(v + 13.0) / 15.0) + 0.06 / (1.0 + std::exp(-(v + 68.0) / 12.0)))) * 20.0;
        dxdt[6] = (h_Kv1Inf - h_Kv1) / h_Kv1Tau;
    }
};