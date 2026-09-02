#include "monte_carlo/simulator.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>

static int check(const char* name, bool ok) {
    std::cout << "  [ " << (ok ? "PASS" : "FAIL") << " ] " << name << std::endl;
    return ok ? 0 : 1;
}

int main() {
    int failures = 0;
    std::cout << "Running Monte Carlo Variance Reduction Engine Unit Tests...\n" << std::endl;

    const monte_carlo::OptionSpec spec{100.0, 105.0, 1.0, 0.03, 0.25};
    const double exact_bs = monte_carlo::MonteCarloEngine::black_scholes_call(spec);

    // Test 1: Plain MC Estimate Bounds
    {
        const auto res = monte_carlo::MonteCarloEngine::run_simulation(2000000, 1, false, false, spec);
        const double err = std::abs(res.price_estimate - exact_bs);
        failures += check("Plain MC estimate within 3 standard errors of Black-Scholes", err < 3.0 * res.standard_error);
    }

    // Test 2: Antithetic Standard Error Reduction
    {
        const auto plain = monte_carlo::MonteCarloEngine::run_simulation(1000000, 2, false, false, spec);
        const auto anti  = monte_carlo::MonteCarloEngine::run_simulation(1000000, 2, true,  false, spec);
        failures += check("Antithetic variates reduces standard error vs plain MC", anti.standard_error < plain.standard_error);
    }

    // Test 3: Control Variate Standard Error Reduction
    {
        const auto plain = monte_carlo::MonteCarloEngine::run_simulation(1000000, 3, false, false, spec);
        const auto ctrl  = monte_carlo::MonteCarloEngine::run_simulation(1000000, 3, false, true,  spec);
        failures += check("Control variate reduces standard error vs plain MC", ctrl.standard_error < plain.standard_error);
        failures += check("Control variate cuts standard error by > 50%", ctrl.standard_error < 0.5 * plain.standard_error);
    }

    // Test 4: Pricing Consistency Across All 4 Schemes
    {
        double worst_diff = 0.0;
        for (int mode = 0; mode < 4; ++mode) {
            const bool anti = (mode & 1);
            const bool ctrl = (mode & 2);
            const auto res = monte_carlo::MonteCarloEngine::run_simulation(1000000, 11, anti, ctrl, spec);
            worst_diff = std::max(worst_diff, std::abs(res.price_estimate - exact_bs));
        }
        failures += check("All 4 variance reduction schemes converge to exact price within 0.02", worst_diff < 0.02);
    }

    // Test 5: 1/sqrt(N) Convergence Law Verification
    {
        const auto res_small = monte_carlo::MonteCarloEngine::run_simulation(100000, 5, false, false, spec);
        const auto res_large = monte_carlo::MonteCarloEngine::run_simulation(1600000, 5, false, false, spec);
        const double se_ratio = res_small.standard_error / res_large.standard_error;
        failures += check("16x path count reduces standard error by ~4x (1/sqrt(N) law)", se_ratio > 3.4 && se_ratio < 4.6);
    }

    // Test 6: Zero Volatility Forward Payoff
    {
        const monte_carlo::OptionSpec zero_vol{100.0, 90.0, 1.0, 0.03, 1e-9};
        const auto res = monte_carlo::MonteCarloEngine::run_simulation(20000, 5, false, false, zero_vol);
        const double expected_forward_payoff = std::exp(-zero_vol.risk_free_rate * zero_vol.time_to_expiry) *
            std::max(zero_vol.spot * std::exp(zero_vol.risk_free_rate * zero_vol.time_to_expiry) - zero_vol.strike, 0.0);
        failures += check("Zero volatility yields exact discounted forward payoff", std::abs(res.price_estimate - expected_forward_payoff) < 1e-6);
    }

    // Test 7: Deep In-The-Money Linear Collapse
    {
        const monte_carlo::OptionSpec deep_itm{100.0, 20.0, 1.0, 0.03, 0.25};
        const double expected_val = deep_itm.spot - deep_itm.strike * std::exp(-deep_itm.risk_free_rate * deep_itm.time_to_expiry);

        const auto plain = monte_carlo::MonteCarloEngine::run_simulation(200000, 5, false, false, deep_itm);
        const auto ctrl  = monte_carlo::MonteCarloEngine::run_simulation(200000, 5, false, true,  deep_itm);

        failures += check("Deep ITM plain MC matches forward less strike within 3 SE", std::abs(plain.price_estimate - expected_val) < 3.0 * plain.standard_error);
        failures += check("Control variate collapses standard error on linear payoff (< 1% of plain SE)", ctrl.standard_error < 0.01 * plain.standard_error);
        failures += check("Controlled deep ITM matches closed form to 1e-6", std::abs(ctrl.price_estimate - expected_val) < 1e-6);
    }

    // Test 8: Seed Reproducibility
    {
        const auto res1 = monte_carlo::MonteCarloEngine::run_simulation(500000, 21, false, false, spec);
        const auto res2 = monte_carlo::MonteCarloEngine::run_simulation(500000, 21, false, false, spec);
        failures += check("Identical PRNG seed produces bitwise identical price estimate", res1.price_estimate == res2.price_estimate);
    }

    std::cout << "\nSummary: " << (failures == 0 ? "ALL TESTS PASSED" : "TEST FAILURES DETECTED") << std::endl;
    return failures == 0 ? 0 : 1;
}
