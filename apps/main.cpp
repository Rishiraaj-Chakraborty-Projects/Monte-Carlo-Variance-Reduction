#include "monte_carlo/simulator.hpp"
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

void print_usage(const char* prog_name) {
    std::cout << "Usage: " << prog_name << " [-S spot] [-K strike] [-T years] [-r rate] [-v vol] [-n paths]\n";
    std::cout << "  -S spot    Underlying asset price (default: 100.0)\n";
    std::cout << "  -K strike  Option strike price (default: 105.0)\n";
    std::cout << "  -T years   Time to expiration in years (default: 1.0)\n";
    std::cout << "  -r rate    Risk-free interest rate (default: 0.03)\n";
    std::cout << "  -v vol     Annualized volatility (default: 0.25)\n";
    std::cout << "  -n paths   Number of simulation paths (default: 1000000)\n";
}

int main(int argc, char** argv) {
    monte_carlo::OptionSpec spec{100.0, 105.0, 1.0, 0.03, 0.25};
    std::uint64_t paths = 1'000'000;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help") {
            print_usage(argv[0]);
            return 0;
        } else if (arg == "-S" && i + 1 < argc) {
            spec.spot = std::stod(argv[++i]);
        } else if (arg == "-K" && i + 1 < argc) {
            spec.strike = std::stod(argv[++i]);
        } else if (arg == "-T" && i + 1 < argc) {
            spec.time_to_expiry = std::stod(argv[++i]);
        } else if (arg == "-r" && i + 1 < argc) {
            spec.risk_free_rate = std::stod(argv[++i]);
        } else if (arg == "-v" && i + 1 < argc) {
            spec.volatility = std::stod(argv[++i]);
        } else if (arg == "-n" && i + 1 < argc) {
            paths = static_cast<std::uint64_t>(std::stoull(argv[++i]));
        }
    }

    const double exact_bs = monte_carlo::MonteCarloEngine::black_scholes_call(spec);

    std::cout << "========================================================\n";
    std::cout << "   MONTE CARLO OPTION PRICER WITH VARIANCE REDUCTION   \n";
    std::cout << "========================================================\n\n";

    std::cout << "Contract Specification:\n";
    std::cout << "  Underlying Spot (S):  $" << std::fixed << std::setprecision(2) << spec.spot << "\n";
    std::cout << "  Strike Price (K):     $" << std::fixed << std::setprecision(2) << spec.strike << "\n";
    std::cout << "  Time to Expiry (T):   " << std::fixed << std::setprecision(2) << spec.time_to_expiry << " years\n";
    std::cout << "  Risk-Free Rate (r):   " << std::fixed << std::setprecision(2) << (spec.risk_free_rate * 100.0) << "%\n";
    std::cout << "  Volatility (sigma):   " << std::fixed << std::setprecision(2) << (spec.volatility * 100.0) << "%\n";
    std::cout << "  Closed-Form Price:    $" << std::fixed << std::setprecision(6) << exact_bs << "\n\n";

    const auto plain = monte_carlo::MonteCarloEngine::run_simulation(paths, 42, false, false, spec);
    const auto anti  = monte_carlo::MonteCarloEngine::run_simulation(paths, 42, true,  false, spec);
    const auto ctrl  = monte_carlo::MonteCarloEngine::run_simulation(paths, 42, false, true,  spec);
    const auto both  = monte_carlo::MonteCarloEngine::run_simulation(paths, 42, true,  true,  spec);

    std::cout << std::left << std::setw(24) << "Scheme"
              << std::right << std::setw(16) << "Price Est ($)"
              << std::setw(16) << "Std Error"
              << std::setw(16) << "Abs Error"
              << std::setw(18) << "Speedup Factor" << "\n";
    std::cout << std::string(90, '-') << "\n";

    auto print_row = [&](const char* name, const monte_carlo::SimulationResult& res) {
        const double variance_ratio = (plain.standard_error / res.standard_error) * (plain.standard_error / res.standard_error);
        std::cout << std::left << std::setw(24) << name
                  << std::right << std::setw(16) << std::fixed << std::setprecision(6) << res.price_estimate
                  << std::setw(16) << std::fixed << std::setprecision(6) << res.standard_error
                  << std::setw(16) << std::scientific << std::setprecision(2) << std::abs(res.price_estimate - exact_bs)
                  << std::setw(18) << std::fixed << std::setprecision(1) << variance_ratio << "x\n";
    };

    print_row("Plain Monte Carlo", plain);
    print_row("Antithetic Variates", anti);
    print_row("Control Variate (S_T)", ctrl);
    print_row("Antithetic + Control", both);

    return 0;
}
