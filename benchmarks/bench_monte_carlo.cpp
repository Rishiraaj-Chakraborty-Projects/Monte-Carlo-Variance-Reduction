#include "monte_carlo/simulator.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>

int main() {
    constexpr std::uint64_t PATHS = 5'000'000;
    std::cout << "Starting Monte Carlo Variance Reduction Benchmark (" << PATHS << " paths)...\n";

    const monte_carlo::OptionSpec spec{100.0, 105.0, 1.0, 0.03, 0.25};

    auto run_bench = [&](const char* name, bool anti, bool ctrl) {
        const auto t0 = std::chrono::high_resolution_clock::now();
        const auto res = monte_carlo::MonteCarloEngine::run_simulation(PATHS, 42, anti, ctrl, spec);
        const auto t1 = std::chrono::high_resolution_clock::now();

        const std::chrono::duration<double, std::milli> duration_ms = t1 - t0;
        const double m_paths_per_sec = (static_cast<double>(PATHS) / duration_ms.count()) * 1e3 / 1e6;

        std::cout << std::left << std::setw(24) << name
                  << std::right << std::setw(16) << std::fixed << std::setprecision(6) << res.price_estimate
                  << std::setw(16) << std::fixed << std::setprecision(6) << res.standard_error
                  << std::setw(16) << std::fixed << std::setprecision(2) << duration_ms.count() << " ms"
                  << std::setw(18) << std::fixed << std::setprecision(1) << m_paths_per_sec << " M paths/s\n";
    };

    std::cout << std::left << std::setw(24) << "Scheme"
              << std::right << std::setw(16) << "Price ($)"
              << std::setw(16) << "Std Error"
              << std::setw(19) << "Latency"
              << std::setw(28) << "Throughput" << "\n";
    std::cout << std::string(103, '-') << "\n";

    run_bench("Plain Monte Carlo", false, false);
    run_bench("Antithetic Variates", true, false);
    run_bench("Control Variate (S_T)", false, true);
    run_bench("Antithetic + Control", true, true);

    return 0;
}
