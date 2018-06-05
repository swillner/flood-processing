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

#include "modules/return_level_lookup.h"
#include "ProgressBar.h"
#include "nvector.h"

namespace flood_processing {
namespace modules {

template<typename T>
void ReturnLevelLookup<T>::run(pipeline::Pipeline* p) {
    const auto return_periods = p->consume<nvector::Vector<T, 3>>("return_periods");
    const auto return_levels_mapping = p->consume<nvector::Vector<T, 3>>("return_levels_mapping");
    const auto return_periods_mapping = p->consume<std::vector<double>>("return_periods_mapping");
    const auto lat_count = return_periods->template size<1>();
    const auto lon_count = return_periods->template size<2>();
    auto return_levels = std::make_shared<nvector::Vector<T, 3>>(std::numeric_limits<T>::quiet_NaN(), return_periods->template size<0>(), lat_count, lon_count);
    ProgressBar progress("Return level lookup", lat_count * lon_count);
    nvector::foreach_split_parallel<nvector::Split<true, false, false>>(
        std::make_tuple(*return_periods, *return_levels_mapping, *return_levels),
        [&](std::size_t lat, std::size_t lon, nvector::View<T, 1>& return_periods_l, nvector::View<T, 1>& return_levels_mapping_l,
            nvector::View<T, 1>& return_levels_l) {
            (void)lat;
            (void)lon;
            if (!std::isnan(return_periods_l(0)) && return_levels_mapping_l(0) >= 0) {
                nvector::foreach_view(std::make_tuple(return_periods_l, return_levels_l), [&](std::size_t t, T& return_period, T& return_level) {
                    (void)t;
                    std::size_t i;
                    for (i = 0; i < return_periods_mapping->size(); ++i) {
                        if (return_period <= (*return_periods_mapping)[i]) {
                            break;
                        }
                    }
                    if (i == 0) {
                        return_level = 0;
                    } else {
                        return_level = return_levels_mapping_l(i - 1);
                    }
                    return true;
                });
            }
            progress.tick();
        });
    p->provide<nvector::Vector<T, 3>>("return_levels", return_levels);
}

template class ReturnLevelLookup<float>;
template class ReturnLevelLookup<double>;

}  // namespace modules
}  // namespace flood_processing
