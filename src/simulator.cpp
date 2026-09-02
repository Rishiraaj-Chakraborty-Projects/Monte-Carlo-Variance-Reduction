#include "monte_carlo/simulator.hpp"
#include <algorithm>
#include <cmath>
#include <random>

namespace monte_carlo {

static double norm_cdf(double x) noexcept {
    return 0.5 * std::erfc(-x * 0.70710678118654752440);
}

double MonteCarloEngine::black_scholes_call(const OptionSpec& s) noexcept {
    const double vol_sqrt_t = s.volatility * std::sqrt(s.time_to_expiry);
    if (vol_sqrt_t < 1e-14) {
        const double df = std::exp(-s.risk_free_rate * s.time_to_expiry);
        const double forward = s.spot * std::exp(s.risk_free_rate * s.time_to_expiry);
        return df * std::max(forward - s.strike, 0.0);
    }
    const double d1 = (std::log(s.spot / s.strike) + (s.risk_free_rate + 0.5 * s.volatility * s.volatility) * s.time_to_expiry) / vol_sqrt_t;
    const double d2 = d1 - vol_sqrt_t;
    return s.spot * norm_cdf(d1) - s.strike * std::exp(-s.risk_free_rate * s.time_to_expiry) * norm_cdf(d2);
}

double MonteCarloEngine::compute_pilot_beta(const OptionSpec& s, unsigned int seed, std::size_t pilot_paths) noexcept {
    std::mt19937_64 rng(seed + 99);
    std::normal_distribution<double> dist(0.0, 1.0);

    const double drift = (s.risk_free_rate - 0.5 * s.volatility * s.volatility) * s.time_to_expiry;
    const double vol = s.volatility * std::sqrt(s.time_to_expiry);
    const double df = std::exp(-s.risk_free_rate * s.time_to_expiry);

    double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_xx = 0.0;
    const double n = static_cast<double>(pilot_paths);

    for (std::size_t i = 0; i < pilot_paths; ++i) {
        const double w = dist(rng);
        const double ST = s.spot * std::exp(drift + vol * w);
        const double y = df * std::max(ST - s.strike, 0.0);

        sum_x += ST;
        sum_y += y;
        sum_xy += ST * y;
        sum_xx += ST * ST;
    }

    const double mean_x = sum_x / n;
    const double mean_y = sum_y / n;
    const double cov_xy = (sum_xy / n) - (mean_x * mean_y);
    const double var_x  = (sum_xx / n) - (mean_x * mean_x);

    return var_x > 0.0 ? cov_xy / var_x : 0.0;
}

SimulationResult MonteCarloEngine::run_simulation(std::uint64_t num_paths,
                                                   unsigned int seed,
                                                   bool use_antithetic,
                                                   bool use_control_variate,
                                                   const OptionSpec& s) noexcept {
    if (num_paths == 0) return {0.0, 0.0, 0};

    std::mt19937_64 rng(seed);
    std::normal_distribution<double> dist(0.0, 1.0);

    const double drift = (s.risk_free_rate - 0.5 * s.volatility * s.volatility) * s.time_to_expiry;
    const double vol = s.volatility * std::sqrt(s.time_to_expiry);
    const double df = std::exp(-s.risk_free_rate * s.time_to_expiry);

    // Known expectation E[S_T] for the control variate
    const double expected_ST = s.spot * std::exp(s.risk_free_rate * s.time_to_expiry);
    const double beta = use_control_variate ? compute_pilot_beta(s, seed, 20000) : 0.0;

    const std::uint64_t draws = use_antithetic ? num_paths / 2 : num_paths;
    double sum = 0.0;
    double sum_sq = 0.0;

    for (std::uint64_t i = 0; i < draws; ++i) {
        const double w = dist(rng);
        double accum_payoff = 0.0;
        const int k = use_antithetic ? 2 : 1;

        for (int j = 0; j < k; ++j) {
            const double z_val = (j == 1) ? -w : w;
            const double ST = s.spot * std::exp(drift + vol * z_val);
            double payoff = df * std::max(ST - s.strike, 0.0);
            if (use_control_variate) {
                payoff -= beta * (ST - expected_ST);
            }
            accum_payoff += payoff;
        }

        const double path_val = accum_payoff / static_cast<double>(k);
        sum += path_val;
        sum_sq += path_val * path_val;
    }

    const double N = static_cast<double>(draws);
    const double mean = sum / N;
    const double variance = std::max((sum_sq / N) - (mean * mean), 0.0);
    const double standard_error = std::sqrt(variance / N);

    return {mean, standard_error, draws};
}

} // namespace monte_carlo
