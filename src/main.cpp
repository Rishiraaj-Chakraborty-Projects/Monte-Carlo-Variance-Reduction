#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <vector>

double ncdf(double x) { return 0.5 * std::erfc(-x * 0.70710678118654752440); }

double bs_call(double S, double K, double T, double r, double sig) {
    const double v = sig * std::sqrt(T);
    const double d1 = (std::log(S / K) + (r + 0.5 * sig * sig) * T) / v;
    return S * ncdf(d1) - K * std::exp(-r * T) * ncdf(d1 - v);
}

struct Res { double est, se; long normals; };

struct Spec { double S, K, T, r, sig; };

// beta from a short pilot run. estimating it on the same paths you price with
// introduces a bias, small at a million paths and not worth the argument.
static double pilot_beta(const Spec& s, unsigned seed, int m) {
    std::mt19937_64 p(seed + 99);
    std::normal_distribution<double> z(0, 1);
    const double drift = (s.r - 0.5 * s.sig * s.sig) * s.T;
    const double vol = s.sig * std::sqrt(s.T);
    const double df = std::exp(-s.r * s.T);
    double sx = 0, sy = 0, sxy = 0, sxx = 0;
    for (int i = 0; i < m; ++i) {
        const double w = z(p);
        const double ST = s.S * std::exp(drift + vol * w);
        const double y = df * std::max(ST - s.K, 0.0);
        sx += ST; sy += y; sxy += ST * y; sxx += ST * ST;
    }
    const double cov = sxy / m - (sx / m) * (sy / m);
    const double var = sxx / m - (sx / m) * (sx / m);
    return var > 0 ? cov / var : 0.0;
}

// n counts normal draws, so every scheme is charged the same for its
// randomness. antithetic uses n/2 draws to build n/2 averaged pairs.
Res mc(long n, unsigned seed, bool antithetic, bool control, const Spec& s) {
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> z(0, 1);
    const double drift = (s.r - 0.5 * s.sig * s.sig) * s.T;
    const double vol = s.sig * std::sqrt(s.T);
    const double df = std::exp(-s.r * s.T);
    const double EST = s.S * std::exp(s.r * s.T);   // known mean of S_T, the control

    const double beta = control ? pilot_beta(s, seed, 20000) : 0.0;
    const long draws = antithetic ? n / 2 : n;

    double sum = 0, sumsq = 0;
    for (long i = 0; i < draws; ++i) {
        const double w = z(rng);
        double acc = 0;
        const int k = antithetic ? 2 : 1;
        for (int j = 0; j < k; ++j) {
            const double ww = j ? -w : w;
            const double ST = s.S * std::exp(drift + vol * ww);
            double y = df * std::max(ST - s.K, 0.0);
            if (control) y -= beta * (ST - EST);
            acc += y;
        }
        const double v = acc / k;
        sum += v; sumsq += v * v;
    }
    const double N = double(draws);
    const double mean = sum / N;
    const double var = std::max(sumsq / N - mean * mean, 0.0);
    return {mean, std::sqrt(var / N), draws};
}

static void table(const Spec& s, long n, double exact) {
    const Res plain = mc(n, 42, false, false, s);
    const Res anti  = mc(n, 42, true,  false, s);
    const Res ctrl  = mc(n, 42, false, true,  s);
    const Res both  = mc(n, 42, true,  true,  s);

    std::printf("  %-22s %11s %11s %11s %10s %9s\n",
                "scheme", "estimate", "std error", "abs error", "normals", "speedup");
    auto row = [&](const char* nm, const Res& x) {
        const double sp = (plain.se / x.se) * (plain.se / x.se);
        std::printf("  %-22s %11.6f %11.6f %11.2e %10ld %8.1fx\n",
                    nm, x.est, x.se, std::abs(x.est - exact), x.normals, sp);
    };
    row("plain", plain);
    row("antithetic", anti);
    row("control variate", ctrl);
    row("both", both);
}

static void convergence(const Spec& s, double exact) {
    std::printf("  %10s %13s %13s %10s %13s\n",
                "paths", "std error", "se x sqrt(n)", "abs error", "in 2 se?");
    for (long n : {10000L, 40000L, 160000L, 640000L, 2560000L}) {
        const Res x = mc(n, 7, false, false, s);
        const double e = std::abs(x.est - exact);
        std::printf("  %10ld %13.6f %13.4f %10.2e %13s\n",
                    n, x.se, x.se * std::sqrt(double(n)), e,
                    e < 2 * x.se ? "yes" : "no");
    }
}

static void moneyness(const Spec& base, long n) {
    std::printf("  %-10s %11s %12s %12s %12s\n",
                "strike", "exact", "plain se", "control se", "reduction");
    for (double K : {70.0, 85.0, 100.0, 105.0, 120.0, 150.0}) {
        Spec s = base; s.K = K;
        const double exact = bs_call(s.S, s.K, s.T, s.r, s.sig);
        const Res p = mc(n, 42, false, false, s);
        const Res c = mc(n, 42, false, true,  s);
        std::printf("  %-10.0f %11.6f %12.6f %12.6f %11.1fx\n",
                    K, exact, p.se, c.se, (p.se / c.se) * (p.se / c.se));
    }
}

static int demo(Spec s, long n) {
    const double exact = bs_call(s.S, s.K, s.T, s.r, s.sig);
    std::printf("European call  S %.1f  K %.1f  T %.2f  r %.3f  sigma %.3f\n",
                s.S, s.K, s.T, s.r, s.sig);
    std::printf("closed form price %.6f, %ld normal draws per scheme\n\n", exact, n);

    table(s, n, exact);
    std::printf("\n  speedup is the factor by which you would have to increase the\n");
    std::printf("  path count to match the variance reduction.\n");

    std::printf("\nwhy \"both\" is not the best row\n\n");
    std::printf("  the two techniques attack the same variance. antithetic pairing\n");
    std::printf("  already removes most of the linear dependence of the payoff on\n");
    std::printf("  the normal draw, and the control variate on S_T is a linear\n");
    std::printf("  correction, so by the time it is applied there is much less\n");
    std::printf("  left for it to remove. stacking variance reduction methods is\n");
    std::printf("  not multiplicative and can be worse than the better one alone.\n");
    std::printf("  measure, do not assume.\n");

    std::printf("\nconvergence of the plain estimator\n\n");
    convergence(s, exact);
    std::printf("\n  se x sqrt(n) is flat, which is the 1/sqrt(n) law. four times the\n");
    std::printf("  paths halves the error. one more decimal place costs a hundred\n");
    std::printf("  times the work.\n");

    std::printf("\ncontrol variate effectiveness by moneyness\n\n");
    moneyness(s, n / 4);
    std::printf("\n  the control is S_T, whose expectation is known exactly. it works\n");
    std::printf("  when the payoff is close to linear in S_T, which is deep in the\n");
    std::printf("  money. far out of the money the payoff is mostly zero and almost\n");
    std::printf("  uncorrelated with S_T, so there is nothing for the control to\n");
    std::printf("  explain. a control variate is only as good as its correlation\n");
    std::printf("  with the thing you are pricing.\n");
    return 0;
}

static int check(const char* name, bool ok, double detail = 0.0) {
    std::printf("  %-50s %-4s %.3e\n", name, ok ? "ok" : "FAIL", detail);
    return ok ? 0 : 1;
}

static int selftest() {
    int bad = 0;
    std::printf("self test\n\n");
    const Spec s{100, 105, 1.0, 0.03, 0.25};
    const double exact = bs_call(s.S, s.K, s.T, s.r, s.sig);

    {
        const Res x = mc(2000000, 1, false, false, s);
        bad += check("plain estimate within 3 standard errors",
                     std::abs(x.est - exact) < 3 * x.se, std::abs(x.est - exact) / x.se);
    }
    {
        const Res a = mc(1000000, 2, false, false, s);
        const Res b = mc(1000000, 2, true, false, s);
        bad += check("antithetic reduces standard error", b.se < a.se, a.se - b.se);
    }
    {
        const Res a = mc(1000000, 3, false, false, s);
        const Res b = mc(1000000, 3, false, true, s);
        bad += check("control variate reduces standard error", b.se < a.se, a.se - b.se);
        bad += check("control variate cuts se by more than half", b.se < 0.5 * a.se,
                     b.se / a.se);
    }
    {
        // every scheme must agree on the price whatever it does to the variance
        double worst = 0;
        for (int mode = 0; mode < 4; ++mode) {
            const Res x = mc(1000000, 11, mode & 1, mode & 2, s);
            worst = std::max(worst, std::abs(x.est - exact));
        }
        bad += check("all four schemes agree on the price", worst < 0.02, worst);
    }
    {
        const Res a = mc(100000, 5, false, false, s);
        const Res b = mc(1600000, 5, false, false, s);
        const double ratio = a.se / b.se;
        bad += check("16x paths cuts se by about 4x", ratio > 3.4 && ratio < 4.6, ratio);
    }
    {
        // zero volatility gives the discounted intrinsic on the forward
        const Spec z{100, 90, 1.0, 0.03, 1e-9};
        const Res x = mc(20000, 5, false, false, z);
        const double want = std::exp(-z.r * z.T) *
                            std::max(z.S * std::exp(z.r * z.T) - z.K, 0.0);
        bad += check("zero vol gives the forward intrinsic",
                     std::abs(x.est - want) < 1e-6, std::abs(x.est - want));
    }
    {
        // deep in the money the payoff is S_T - K with probability ~1, so the
        // call is the forward less the discounted strike. plain MC still has
        // the full variance of S_T, so this is checked in standard errors.
        const Spec d{100, 20, 1.0, 0.03, 0.25};
        const double want = d.S - d.K * std::exp(-d.r * d.T);
        const Res x = mc(200000, 5, false, false, d);
        bad += check("deep ITM call within 3 se of forward less strike",
                     std::abs(x.est - want) < 3 * x.se, std::abs(x.est - want) / x.se);
        // with the payoff exactly linear in S_T the control variate should
        // remove essentially all of the variance
        const Res c = mc(200000, 5, false, true, d);
        bad += check("control variate collapses se when payoff is linear",
                     c.se < 0.01 * x.se, c.se / x.se);
        bad += check("controlled deep ITM matches closed form to 1e-6",
                     std::abs(c.est - want) < 1e-6, std::abs(c.est - want));
    }
    {
        const Spec f{100, 400, 1.0, 0.03, 0.25};
        const Res x = mc(200000, 5, false, false, f);
        bad += check("far OTM call is non negative and tiny",
                     x.est >= 0.0 && x.est < 0.01, x.est);
    }
    {
        const Res a = mc(500000, 21, false, false, s);
        const Res b = mc(500000, 21, false, false, s);
        bad += check("same seed reproduces exactly", a.est == b.est,
                     std::abs(a.est - b.est));
    }
    {
        const Res a = mc(500000, 31, false, false, s);
        const Res b = mc(500000, 32, false, false, s);
        bad += check("different seeds differ but agree within 4 se",
                     a.est != b.est && std::abs(a.est - b.est) < 4 * (a.se + b.se),
                     std::abs(a.est - b.est));
    }
    std::printf("\n%s\n", bad ? "FAILURES PRESENT" : "all checks passed");
    return bad ? 1 : 0;
}

int main(int argc, char** argv) {
    Spec s{100, 105, 1.0, 0.03, 0.25};
    long n = 1000000;
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--test")) return selftest();
        if (!std::strcmp(argv[i], "--help")) {
            std::printf("usage: mc [-S spot] [-K strike] [-T years] [-r rate] "
                        "[-v vol] [-n paths] [--test]\n");
            return 0;
        }
        if (i + 1 >= argc) continue;
        if (!std::strcmp(argv[i], "-S")) s.S = std::atof(argv[++i]);
        else if (!std::strcmp(argv[i], "-K")) s.K = std::atof(argv[++i]);
        else if (!std::strcmp(argv[i], "-T")) s.T = std::atof(argv[++i]);
        else if (!std::strcmp(argv[i], "-r")) s.r = std::atof(argv[++i]);
        else if (!std::strcmp(argv[i], "-v")) s.sig = std::atof(argv[++i]);
        else if (!std::strcmp(argv[i], "-n")) n = std::atol(argv[++i]);
    }
    return demo(s, n);
}
