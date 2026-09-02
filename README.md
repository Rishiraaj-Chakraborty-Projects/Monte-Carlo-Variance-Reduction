# Monte Carlo variance reduction

European call priced by simulation, with antithetic variates and a control variate, measured rather than asserted. C++17, single file, no dependencies.

```
  scheme                    estimate   std error   abs error    normals   speedup
  plain                     9.145571    0.016384    2.38e-02    1000000      1.0x
  antithetic                9.122171    0.013592    3.71e-04     500000      1.5x
  control variate           9.121406    0.007601    3.94e-04    1000000      4.6x
  both                      9.119993    0.010178    1.81e-03     500000      2.6x
```

Closed form price 9.121799. Speedup is the factor by which you would have to raise the path count to match the variance reduction, which is the honest way to quote it: halving the standard error is worth four times the paths.

## The result worth reading

**"Both" is worse than the control variate alone.** 2.6x against 4.6x.

That is not a bug and it is the most useful line in the table. The two techniques attack the same variance.

Antithetic pairing averages the payoff at `+w` and `-w`. That removes essentially all of the component of the payoff that is odd in the normal draw, which is most of its linear dependence on the draw. The control variate on `S_T` is itself a linear correction. By the time it runs, antithetic has already taken the part it was going to take, so it has much less left to remove. Meanwhile the pairing halves the number of independent averaged observations.

**Variance reduction methods do not compose multiplicatively, and stacking them can be worse than using the better one alone.** Measure the combination; do not assume it.

## Charging every scheme the same

`n` counts **normal draws**, not payoff evaluations. Antithetic uses `n/2` draws to build `n/2` averaged pairs, and the standard error is computed on those `n/2` averages.

This matters. If you count payoffs instead, antithetic gets `n` payoffs from `n/2` draws and looks twice as good as it is. The randomness is the scarce resource, so it is what gets budgeted. The `normals` column makes the accounting visible.

## Convergence

```
       paths     std error  se x sqrt(n)  abs error      in 2 se?
       10000      0.159595       15.9595   2.62e-01           yes
       40000      0.081041       16.2083   1.24e-01           yes
      160000      0.040836       16.3344   3.61e-02           yes
      640000      0.020431       16.3452   7.00e-03           yes
     2560000      0.010194       16.3107   1.45e-02           yes
```

The middle column is `se · sqrt(n)`, and it is flat to three digits across a 256-fold range of path counts. That is the `1/sqrt(n)` law made visible: the standard error is a constant divided by the square root of the sample size, and the constant is a property of the payoff, not of the simulation.

The practical reading is bleak and worth internalising. Four times the paths halves the error. One additional decimal place costs a hundred times the work. Monte Carlo is the method of last resort for exactly this reason, and it wins only where dimensionality makes everything else impossible.

The last column checks the actual error against two standard errors. It should say yes about 95 percent of the time, and a row that says no is not automatically a failure.

## When a control variate works, and when it does not

```
  strike           exact     plain se   control se    reduction
  70           32.608155     0.049221     0.004538       117.7x
  85           20.375374     0.044235     0.010181        18.9x
  100          11.348477     0.035906     0.014531         6.1x
  105           9.121799     0.032785     0.015213         4.6x
  120           4.463302     0.023705     0.014932         2.5x
  150           0.875272     0.010672     0.009175         1.4x
```

The control is `S_T`, whose risk-neutral expectation is known exactly to be `S·exp(rT)`. The correction is

```
y  ->  y - beta·(S_T - E[S_T])
```

with `beta` estimated from a short pilot run.

Deep in the money the payoff is `S_T - K` with probability near one, so it is **exactly linear** in the control. Correlation approaches one and the control removes nearly all the variance: a 118-fold reduction at a strike of 70.

Far out of the money the payoff is zero on almost every path. It is barely correlated with `S_T` at all, so there is nothing for the control to explain, and the reduction falls to 1.4x.

The general statement: **the variance reduction from a control variate is `1/(1-rho²)`, where rho is the correlation between the payoff and the control.** A control with correlation 0.99 gives 50x. A control with correlation 0.5 gives 1.33x and is not worth the code. Choosing a control is choosing something highly correlated with the payoff whose expectation you happen to know in closed form, and there is no general recipe for that.

The self test pushes this to its limit: at a strike of 20 the payoff is linear in `S_T` for practical purposes, and the control variate cuts the standard error by a factor of 120,000, matching the closed form to `2e-10`.

## Build

PowerShell:

```powershell
g++ -std=c++17 -O2 -o mc.exe src\main.cpp
.\mc.exe
```

CMake:

```
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

## Run

```
mc                                   default: S=100 K=105 T=1 r=3% sigma=25%
mc -K 80 -v 0.4 -n 5000000           different option, more paths
mc -S 50 -K 55 -T 0.25               short dated
mc --test                            13 checks, exit code 0 or 1
mc --help
```

## Self test

```
  plain estimate within 3 standard errors            ok
  antithetic reduces standard error                  ok
  control variate reduces standard error             ok
  control variate cuts se by more than half          ok
  all four schemes agree on the price                ok
  16x paths cuts se by about 4x                      ok
  zero vol gives the forward intrinsic               ok
  deep ITM call within 3 se of forward less strike   ok
  control variate collapses se when payoff is linear ok
  controlled deep ITM matches closed form to 1e-6    ok
  far OTM call is non negative and tiny              ok
  same seed reproduces exactly                       ok
  different seeds differ but agree within 4 se       ok
```

Two of these are worth pointing at.

**"All four schemes agree on the price."** Variance reduction changes the variance, never the expectation. A scheme that produced a tighter standard error around the wrong number would pass every other check in the list and fail this one. It is the check that separates variance reduction from bias.

**"16x paths cuts se by about 4x."** A direct test of the `sqrt(n)` law, tolerance 3.4 to 4.6, which is loose enough to survive sampling noise and tight enough to catch an estimator that is not converging at the right rate.

The deep-ITM check is stated in **standard errors** rather than absolute price. An earlier version used a fixed `1e-4` tolerance and failed, because with 200,000 paths the plain estimator has a standard error near 0.06 and the test was demanding six hundred times more precision than the method delivers. The test was wrong, not the code, and stating tolerances in standard errors is the fix that generalises.

## Notes

`beta` is estimated from a separate pilot run with its own seed rather than from the pricing paths. Fitting it on the same paths you price with introduces a small downward bias in the reported standard error. At a million paths it is negligible; the separate pilot avoids the argument entirely and costs 20,000 draws.

The random source is `std::mt19937_64` with `std::normal_distribution`, which uses Box-Muller or a Ziggurat depending on the standard library. For serious work you would want a stream that is reproducible across implementations, which this is not: the same seed gives the same answer on the same library, and nothing is promised across compilers.

## Licence

MIT.
