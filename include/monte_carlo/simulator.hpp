#ifndef MONTE_CARLO_SIMULATOR_HPP
#define MONTE_CARLO_SIMULATOR_HPP

#include <cmath>
#include <cstdint>
#include <cstddef>
#include <vector>

namespace monte_carlo {

struct OptionSpec {
    double spot{100.0};            // Asset spot price S
    double strike{105.0};          // Strike price K
    double time_to_expiry{1.0};    // Time to expiration in years T
    double risk_free_rate{0.03};   // Risk-free interest rate r
    double volatility{0.25};       // Annualized volatility sigma
};

struct SimulationResult {
    double price_estimate{0.0};
    double standard_error{0.0};
    std::uint64_t normal_draws{0};
};

/**
 * Monte Carlo Option Pricing Engine with Variance Reduction Techniques:
 * 1. Antithetic Variates
 * 2. Control Variates (using underlying terminal asset S_T)
 */
class MonteCarloEngine {
public:
    /**
     * Runs Monte Carlo simulation pricing European Call options.
     */
    [[nodiscard]] static SimulationResult run_simulation(std::uint64_t num_paths,
                                                          unsigned int seed,
                                                          bool use_antithetic,
                                                          bool use_control_variate,
                                                          const OptionSpec& spec) noexcept;

    /**
     * Estimates control variate optimal beta coefficient via a pilot simulation.
     */
    [[nodiscard]] static double compute_pilot_beta(const OptionSpec& spec,
                                                   unsigned int seed,
                                                   std::size_t pilot_paths = 20000) noexcept;

    /**
     * Analytical Black-Scholes European Call closed-form reference.
     */
    [[nodiscard]] static double black_scholes_call(const OptionSpec& spec) noexcept;
};

} // namespace monte_carlo

#endif // MONTE_CARLO_SIMULATOR_HPP
