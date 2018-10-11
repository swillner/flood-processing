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

#include "modules/return_period_threshold.h"
#include "progressbar.h"
#include "nvector.h"

namespace flood_processing {
namespace modules {

template<typename T>
void ReturnPeriodThreshold<T>::run(pipeline::Pipeline* p) {
    const auto return_periods = p->consume<nvector::Vector<T, 3>>("return_periods");
    const auto return_levels = p->consume<nvector::Vector<T, 3>>("return_levels");
    const auto size = return_periods->template size<0>();
    const auto lat_count = return_periods->template size<1>();
    const auto lon_count = return_periods->template size<2>();
    auto return_levels_thresholded = std::make_shared<nvector::Vector<T, 3>>(std::numeric_limits<T>::quiet_NaN(), size, lat_count, lon_count);
    std::shared_ptr<nvector::Vector<T, 2>> raster;
    if (use_raster) {
        raster = p->consume<nvector::Vector<T, 2>>("raster");
    }

    progressbar::ProgressBar progress(lat_count * lon_count, "Return period threshold");
    nvector::foreach_split_parallel<nvector::Split<true, false, false>>(
        std::make_tuple(*return_periods, *return_levels, *return_levels_thresholded),
        [&](std::size_t lat, std::size_t lon, nvector::View<T, 1>& return_periods_l, nvector::View<T, 1>& return_levels_l,
            nvector::View<T, 1>& return_levels_thresholded_l) {
            if (!std::isnan(return_periods_l(0)) && !std::isnan(return_levels_l(0))) {
                T this_threshold = threshold;
                if (use_raster) {
                    this_threshold = (*raster)(lat, lon);
                }
                nvector::foreach_view(std::make_tuple(return_periods_l, return_levels_l, return_levels_thresholded_l),
                                      [&](std::size_t t, T& return_period, T& return_level, T& return_level_thresholded) {
                                          (void)t;
                                          if (return_period >= this_threshold) {
                                              return_level_thresholded = return_level;
                                          } else {
                                              return_level_thresholded = 0;
                                          }
                                          return true;
                                      });
            }
            ++progress;
        });
    p->provide<nvector::Vector<T, 3>>("return_levels_thresholded", return_levels_thresholded);
}

template class ReturnPeriodThreshold<float>;
template class ReturnPeriodThreshold<double>;

}  // namespace modules
}  // namespace flood_processing
