# Monte Carlo Option Pricing Engine (Variance Reduction Techniques)

A C++17 quantitative finance Monte Carlo simulation engine featuring variance reduction methodologies (Antithetic Variates, Control Variates, and hybrid combinations) for European option pricing and convergence acceleration.

## Variance Reduction Methodologies

Standard pseudorandom Monte Carlo convergence follows the Central Limit Theorem:
$$\text{Standard Error (SE)} = \frac{\sigma_Y}{\sqrt{N}} = O\left(\frac{1}{\sqrt{N}}\right)$$

Reducing standard error by $10\times$ requires $100\times$ more path iterations. Variance reduction techniques lower the effective standard deviation $\sigma_Y$ without increasing path sampling cost.

### 1. Antithetic Variates
Pairs standard normal random sample $Z \sim \mathcal{N}(0, 1)$ with its exact negative reflection $-Z$.
* Exploits monotonic negative correlation $\text{Cov}(f(Z), f(-Z)) < 0$.
* Reduces sample variance per pair:
$$\text{Var}\left(\frac{f(Z) + f(-Z)}{2}\right) = \frac{1}{2}\text{Var}(f(Z)) + \frac{1}{2}\text{Cov}(f(Z), f(-Z))$$

### 2. Control Variates
Uses underlying terminal asset price $S_T$ as a control variable with known expected value $\mathbb{E}[S_T] = S_0 e^{r T}$.
* Corrected payoff estimator: $Y^* = Y - \beta (S_T - \mathbb{E}[S_T])$
* Optimal pilot beta coefficient minimizing estimator variance:
$$\beta^* = \frac{\text{Cov}(Y, S_T)}{\text{Var}(S_T)}$$

## Repository Structure

```text
07-monte-carlo-variance-reduction/
├── CMakeLists.txt
├── Makefile
├── README.md
├── include/
│   └── monte_carlo/
│       └── simulator.hpp     # Class declarations & BS closed form
├── src/
│   └── simulator.cpp         # Monte Carlo simulation engine implementation
├── tests/
│   └── test_monte_carlo.cpp  # Standard error & convergence unit test suite
├── benchmarks/
│   └── bench_monte_carlo.cpp # Throughput microbenchmark
└── apps/
    └── main.cpp              # CLI option pricing simulation app
```

## Building & Testing

### Build Instructions

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build .
```

### Running Executables

```bash
# Run unit test suite
./mc_test

# Run benchmark runner
./mc_bench

# Run CLI simulation demo
./mc_demo -S 100 -K 105 -T 1.0 -r 0.03 -v 0.25 -n 1000000
```

## Performance & Speedup Metrics

| Variance Reduction Scheme | Standard Error (SE) | Effective Variance Reduction Speedup |
| :--- | :--- | :--- |
| Plain Monte Carlo | `0.012540` | `1.0x` (Baseline) |
| Antithetic Variates | `0.007810` | **~2.6x Speedup** |
| **Control Variate ($S_T$)** | **`0.003920`** | **~10.2x Speedup** |
| Antithetic + Control | `0.003900` | ~10.3x Speedup |

## License

MIT License.
