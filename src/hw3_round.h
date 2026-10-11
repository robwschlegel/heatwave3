#ifndef HW3_ROUND_H
#define HW3_ROUND_H

// R-compatible rounding and summaries, so heatwave3 computes exactly as
// heatwaveR (which uses R's round(), mean(), sum() and sd()). Ports of R 4.x
// src/nmath/fround.c (fround with R_pow_di), of R's mean and sum for doubles
// (summary.c: long double accumulation, mean with one refinement pass) and of
// stats::var (cov.c), on which sd() is built.
//
// R's round(x, d) picks whichever of floor(x * 10^d) / 10^d and
// ceil(x * 10^d) / 10^d is closer to x, the even one on an exact tie. It
// differs from std::round(x * 10^d) / 10^d for numbers such as 0.15
// (R: round(0.15, 1) == 0.1).

#include <cfloat>
#include <cmath>
#include <limits>

namespace hw3 {

// R's R_pow_di(): x^n by repeated squaring
inline double r_pow_di(double x, int n) {
    double pow = 1.0;
    if (std::isnan(x)) return x;
    if (n != 0) {
        if (!std::isfinite(x)) return std::pow(x, static_cast<double>(n));
        if (n < 0) { n = -n; x = 1 / x; }
        for (;;) {
            if (n & 01) pow *= x;
            if (n >>= 1) x *= x; else break;
        }
    }
    return pow;
}

// R's round(x, digits)
inline double r_round(double x, int digits) {
    const int max10e = DBL_MAX_10_EXP;
    if (std::isnan(x)) return x;
    if (!std::isfinite(x)) return x;
    if (digits > max10e + DBL_DIG || x == 0.) return x;
    if (digits < -max10e) return 0.;
    if (digits == 0) return std::nearbyint(x);

    int dig = digits;
    double sgn = +1.;
    if (x < 0.) { sgn = -1.; x = -x; }
    double l10x = 0.301029995663981195213738894724 * (0.5 + std::logb(x));  // M_LOG10_2
    if (l10x + dig > DBL_DIG) return sgn * x;

    double pow10, x10, i10, xd, xu;
    if (dig <= max10e) {
        pow10 = r_pow_di(10., dig);
        x10 = x * pow10;
        i10 = std::floor(x10);
        xd = i10 / pow10;
        xu = std::ceil(x10) / pow10;
    } else {
        int e10 = dig - max10e;
        double p10 = r_pow_di(10., e10);
        pow10 = r_pow_di(10., max10e);
        x10 = (x * pow10) * p10;
        i10 = std::floor(x10);
        xd = i10 / pow10 / p10;
        xu = std::ceil(x10) / pow10 / p10;
    }
    double du = xu - x, dd = x - xd;
    return sgn * ((du < dd || (std::fmod(i10, 2.) == 1 && du == dd)) ? xu : xd);
}

// R's mean() of n doubles (no missing values): long double sum, then one
// refinement pass, as R's summary.c
inline double r_mean(const double* x, int n) {
    long double s = 0.0;
    for (int i = 0; i < n; ++i) s += x[i];
    s /= n;
    if (std::isfinite(static_cast<double>(s))) {
        long double t = 0.0;
        for (int i = 0; i < n; ++i) t += (x[i] - s);
        s += t / n;
    }
    return static_cast<double>(s);
}

// R's sum() of n doubles: long double accumulation
inline double r_sum(const double* x, int n) {
    long double s = 0.0;
    for (int i = 0; i < n; ++i) s += x[i];
    return static_cast<double>(s);
}

// R's sd(): sqrt of var(), which takes the refined mean (as double), then
// sums (x - mean)^2 in long double (each square in double) over n - 1;
// NaN (R's NA) for fewer than two values
inline double r_sd(const double* x, int n) {
    if (n < 2) return std::numeric_limits<double>::quiet_NaN();
    double xm = r_mean(x, n);
    long double s = 0.0;
    for (int i = 0; i < n; ++i) {
        double d = x[i] - xm;
        s += d * d;
    }
    return std::sqrt(static_cast<double>(s / (n - 1)));
}

} // namespace hw3

#endif
