#ifndef LMOMENTS_GUM_H
#define LMOMENTS_GUM_H

#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <vector>

#include "utils.h"

namespace lmoments {

template<typename T>
class GUM : public distribution<T> {
  protected:
    T u = 0;
    T a = 0;

  public:
    GUM(){};
    template<typename DataVector = std::vector<T>>
    explicit GUM(const DataVector& x) { fit_data(x); }

    template<typename DataVector = std::vector<T>>
    inline void fit_data(const DataVector& x) {
        std::vector<T> xmom(2);
        lmoment_ratios_unbiased(x, xmom);
        set_lmoment_ratios(xmom);
    }

    inline void set_lmoment_ratios(const std::vector<T>& xmom) override {
        static const T eu = .577215664901532861;   // Euler's Constant
        static const T dl2 = .693147180559945309;  // log(2)
        if (xmom[1] <= 0) {
            throw std::invalid_argument("lmoments invalid");
        }
        a = xmom[1] / dl2;
        u = xmom[0] - eu * a;
    }

    inline void get_lmoment_ratios(std::vector<T>& xmom) const override {
        const std::size_t nmom = xmom.size();
        static const T zmom[20] = {.577215664901532861,   .693147180559945309,   .169925001442312363,   .150374992788438185,   .0558683500577583138,
                                   .0581100239999710876,  .0276242584297309125,  .0305563766579053126,  .0164650282258328802,  .0187846624298170912,
                                   .0109328215063027148,  .012697312667632953,   .00778982818057231804, .00914836179621999726, .00583332389328363588,
                                   .00690104287590348154, .00453267970180679549, .00538916811326595459, .0036240776777236879,  .00432387608605538096};

        if (nmom > 20) {
            throw std::invalid_argument("too many lmoment ratios requested");
        }
        xmom[0] = u + a * zmom[0];
        if (nmom > 1) {
            xmom[1] = a * zmom[1];
            if (nmom > 2) {
                for (std::size_t j = 2; j < nmom; ++j) {
                    xmom[j] = zmom[j];
                }
            }
        }
    }

    inline T cdf(T x) const override {
        const T y = (x - u) / a;
        return std::exp(-std::exp(-y));
    }

    inline T quantile(T f) const override {
        if (f <= 0 || f >= 1) {
            throw std::invalid_argument("argument out of range");
        }
        return u - a * std::log(-std::log(f));
    }

    inline std::vector<T> get_parameters() const override { return std::vector<T>({u, a}); }
};
}  // namespace lmoments

#endif
