#ifndef LMOMENTS_UTILS_H
#define LMOMENTS_UTILS_H

#include <cmath>
#include <stdexcept>
#include <vector>

namespace lmoments {

template<typename T>
inline T dlgama(T x) {  // natural logarithm of the gamma function
    static const T small = 1e-7;
    static const T c5 = 8.41750841750841751e-4;
    static const T c6 = -.00191752691752691753;
    static const T c7 = .00641025641025641026;
    static const T s1 = -.577215664901532861;  // -Euler's Constant
    static const T s2 = .822467033424113218;   // pi^2/12
    static const T crit = 13.;
    static const T big = 1e9;
    static const T toobig = 2e36;
    static const T c0 = .918938533204672742;
    static const T c1 = .0833333333333333333;
    static const T c2 = -.00277777777777777778;
    static const T c3 = 7.93650793650793651e-4;
    static const T c4 = -5.95238095238095238e-4;

    if (x <= 0 || x > toobig) {
        throw std::invalid_argument("argument out of range");
    }
    T xx;
    if (std::abs(x - 2) <= small) {
        xx = x - 2;
        return std::log(x - 1) + xx * (s1 + xx * s2);
    }
    if (std::abs(x - 1) <= small) {
        xx = x - 1;
        return xx * (s1 + xx * s2);
    }
    if (x <= small) {
        return -std::log(x) + s1 * x;
    }
    T sum1 = 0;
    T y = x;
    if (y < crit) {
        T z__ = 1;
        while (y < crit) {
            z__ *= y;
            y += 1;
        }
        sum1 -= std::log(z__);
    }
    sum1 = sum1 + (y - 0.5) * std::log(y) - y + c0;
    if (y >= big) {
        return sum1;
    }
    T z__ = 1 / (y * y);
    T sum2 = ((((((c7 * z__ + c6) * z__ + c5) * z__ + c4) * z__ + c3) * z__ + c2) * z__ + c1) / y;
    return sum1 + sum2;
}

template<typename T, typename DataVector = std::vector<T>>
inline void lmoment_ratios_unbiased(const DataVector& x, std::vector<T>& xmom) {  // samlmu
    const std::size_t nmom = xmom.size();
    const int n = x.size();  // type must be signed
    std::vector<T> coef1(nmom);
    std::vector<T> coef2(nmom);

    if (nmom > 100) {
        throw std::invalid_argument("too many lmoment ratios requested");
    }
    for (std::size_t j = 0; j < nmom; ++j) {
        xmom[j] = 0.;
    }
    if (nmom <= 2) {
        T sum1 = 0.;
        T sum2 = 0.;
        T temp = -n + 1;
        for (int i = 0; i < n; ++i) {
            sum1 += x[i];
            sum2 += x[i] * temp;
            temp += 2.;
        }
        xmom[0] = sum1 / n;
        if (nmom == 2) {
            xmom[1] = sum2 / (n * (n - 1));
        }
    } else {
        for (std::size_t j = 2; j < nmom; ++j) {
            T temp = 1. / static_cast<T>(j * (n - j));
            coef1[j] = (2 * j - 1) * temp;
            coef2[j] = ((j - 1) * (n + j - 1)) * temp;
        }
        T temp = -n - 1;
        std::size_t nhalf = n / 2;
        for (std::size_t i = 0; i < nhalf; ++i) {
            temp += 2.;
            T xi = x[i];
            T xii = x[n - 1 - i];
            T termp = xi + xii;
            T termn = xi - xii;
            xmom[0] += termp;
            T s1 = 1.;
            T s = temp / (n - 1);
            xmom[1] += s * termn;
            for (std::size_t j = 2; j < nmom; j += 2) {
                T s2 = s1;
                s1 = s;
                s = coef1[j] * temp * s1 - coef2[j] * s2;
                xmom[j] += s * termp;
                if (j < nmom - 1) {
                    s2 = s1;
                    s1 = s;
                    s = coef1[j + 1] * temp * s1 - coef2[j + 1] * s2;
                    xmom[j + 1] += s * termn;
                }
            }
        }
        if (n % 2 == 1) {
            T term = x[nhalf];
            T s = 1.;
            xmom[0] += term;
            for (std::size_t j = 2; j < nmom; j += 2) {
                s = -coef2[j] * s;
                xmom[j] += s * term;
            }
        }
        xmom[0] /= n;
        if (xmom[1] == 0.) {  // all data values are equal
            for (std::size_t j = 2; j < nmom; ++j) {
                xmom[j] = 0.;
            }
        } else {
            for (std::size_t j = 2; j < nmom; ++j) {
                xmom[j] /= xmom[1];
            }
            xmom[1] /= n;
        }
    }
}

template<typename T>
class distribution {
  public:
    virtual void set_lmoment_ratios(const std::vector<T>& xmom) = 0;  // pel*
    virtual void get_lmoment_ratios(std::vector<T>& xmom) const = 0;  // lmr*
    virtual T cdf(T x) const = 0;                                     // cdf*
    virtual T quantile(T f) const = 0;                                // qua*
    virtual std::vector<T> get_parameters() const = 0;
    virtual ~distribution() = default;
};

}  // namespace lmoments

#endif
