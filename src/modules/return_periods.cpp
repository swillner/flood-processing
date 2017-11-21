/*
  Copyright (C) 2017 Sven Willner <sven.willner@pik-potsdam.de>

  This file is part of flood-processing.

  flood-processing is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as
  published by the Free Software Foundation, either version 3 of
  the License, or (at your option) any later version.

  flood-processing is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with flood-processing.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "modules/return_periods.h"
#include "lmoments.h"
#include "nvector.h"
#include "settingsnode.h"
#ifdef FLOOD_PROCESSING_WITH_TQDM
#include "tqdm/tqdm.h"
#endif

namespace flood_processing {
namespace modules {

// http://www.itl.nist.gov/div898/handbook/eda/section3/eda35e.htm
template<typename T>
template<typename Distribution>
inline T ReturnPeriods<T>::anderson_darling_test(const std::vector<T>& data, const Distribution& d) {
    // TODO sort data
    const auto N = static_cast<T>(data.size());
    T res = N;
    for (std::size_t i = 0; i < data.size(); ++i) {
        res += (2 * i + 1) * (std::log(d.cdf(data[i])) + std::log(1 - d.cdf(data[N - 1 - i]))) / N;
    }
    return std::sqrt(-res);
}

// http://www.itl.nist.gov/div898/handbook/eda/section3/eda35g.htm
template<typename T>
template<typename Distribution>
inline T ReturnPeriods<T>::kolmogorov_smirnov_test(const std::vector<T>& data, const Distribution& d) {
    T res = 0;
    const auto N = static_cast<T>(data.size());
    for (std::size_t i = 0; i < data.size(); ++i) {
        res = std::max(res, std::max(d.cdf(data[i]) - i / N, (i + 1) / N - d.cdf(data[i])));
    }
    return res;
}

template<typename T>
nvector::Vector<T, 3> ReturnPeriods<T>::return_periods(nvector::Vector<T, 3>& history_discharge,
                                                       nvector::Vector<T, 3>& projection_discharge,
                                                       std::size_t from,
                                                       std::size_t to) {
    const auto size = projection_discharge.template size<0>();
    const auto lat_count = projection_discharge.template size<1>();
    const auto lon_count = projection_discharge.template size<2>();
    if (from == 0 && to == 0) {
        to = size - 1;
    }
    if (to < from || to >= history_discharge.template size<0>()) {
        throw std::runtime_error("return period parameters invalid");
    }
    if (lat_count != history_discharge.template size<1>() || lon_count != history_discharge.template size<2>()) {
        throw std::runtime_error("return period grids differ");
    }
#ifdef FLOOD_PROCESSING_WITH_TQDM
    const auto total = lat_count * lon_count;
    tqdm::Params p;
    p.desc = "Return periods";
    p.ascii = "";
    p.f = stdout;
    tqdm::RangeTqdm<int> it{tqdm::RangeIterator<int>(total), tqdm::RangeIterator<int>(total, total), p};
#endif
    auto result_grid = nvector::Vector<T, 3>(std::numeric_limits<T>::quiet_NaN(), size, lat_count, lon_count);
    nvector::foreach_split_parallel<nvector::Split<true, false, false>>(
        std::make_tuple(history_discharge, projection_discharge, result_grid),
        [&](std::size_t lat, std::size_t lon, nvector::View<T, 1>& history_series, nvector::View<T, 1>& projection_series, nvector::View<T, 1>& result_series) {
            (void)lat;
            (void)lon;
            if (history_series(0) >= 0 && history_series(0) < 1e10 && !std::isnan(history_series(0))) {
                try {
                    std::vector<T> view(to - from + 1);
                    std::partial_sort_copy(std::begin(history_series) + from, std::begin(history_series) + (to + 1), std::begin(view), std::end(view));
                    lmoments::GEV<T> d(view);
                    for (std::size_t i = 0; i < size; ++i) {
                        T return_period = 1. / (1. - d.cdf(projection_series(i)));
                        if (return_period > 1000) {
                            return_period = 1000;
                        }
                        result_series(i) = return_period;
                    }
                } catch (std::invalid_argument&) {
                    for (std::size_t i = 0; i < size; ++i) {
                        result_series(i) = 0;
                    }
                }
            }
#ifdef FLOOD_PROCESSING_WITH_TQDM
#pragma omp critical(output)
            { ++it; }
#endif
        });
#ifdef FLOOD_PROCESSING_WITH_TQDM
    it.close();
#endif
    return result_grid;
}

template<typename T>
ReturnPeriods<T>::ReturnPeriods(const settings::SettingsNode& settings) {
    from = settings["from"].as<std::size_t>();
    to = settings["to"].as<std::size_t>();
}

template<typename T>
void ReturnPeriods<T>::run(pipeline::Pipeline* p) {
    auto history_discharge = p->consume<nvector::Vector<T, 3>>("history_discharge");
    auto projection_discharge = p->consume<nvector::Vector<T, 3>>("projection_discharge");
    p->provide<nvector::Vector<T, 3>>("return_periods",
                                      std::make_shared<nvector::Vector<T, 3>>(return_periods(*history_discharge, *projection_discharge, from, to)));
}

template class ReturnPeriods<double>;
template class ReturnPeriods<float>;

}  // namespace modules
}  // namespace flood_processing