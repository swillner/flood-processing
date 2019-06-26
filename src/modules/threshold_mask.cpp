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

#include "modules/threshold_mask.h"
#include "nvector.h"
#include "progressbar.h"

namespace flood_processing {
namespace modules {

template<typename T>
void ThresholdMask<T>::run(pipeline::Pipeline* p) {
    const auto compare_to = p->consume<nvector::Vector<T, 2>>(compare_to_name);
    const auto apply_to = p->consume<nvector::Vector<T, 3>>(apply_to_name);
    const auto size = apply_to->template size<0>();
    const auto lat_count = apply_to->template size<1>();
    const auto lon_count = apply_to->template size<2>();
    auto output = std::make_shared<nvector::Vector<T, 3>>(std::numeric_limits<T>::quiet_NaN(), size, lat_count, lon_count);
    std::shared_ptr<nvector::Vector<T, 3>> raster;
    if (use_raster) {
        raster = p->consume<nvector::Vector<T, 3>>(raster_name);
    }

    progressbar::ProgressBar progress(lat_count * lon_count, "Threshold mask");
    nvector::foreach_split_parallel<nvector::Split<true, false, false>>(
        std::make_tuple(*apply_to, *output), [&](std::size_t lat, std::size_t lon, nvector::View<T, 1>& apply_to_l, nvector::View<T, 1>& output_l) {
            const auto c = (*compare_to)(lat, lon);
            if (!std::isnan(c) && !std::isnan(apply_to_l(0))) {
                nvector::foreach_view(std::make_tuple(apply_to_l, output_l), [&](std::size_t t, T& a, T& o) {
                    T this_threshold = threshold;
                    if (use_raster) {
                        this_threshold = (*raster)(t, lat, lon);
                    }
                    if ((c < this_threshold && less_than) || (c > this_threshold && !less_than) || (c == this_threshold && or_equal)) {
                        o = a;
                    } else {
                        o = mask_value;
                    }
                    return true;
                });
            }
            ++progress;
        });
    p->provide<nvector::Vector<T, 3>>(output_name, output);
}

template class ThresholdMask<float>;
template class ThresholdMask<double>;

}  // namespace modules
}  // namespace flood_processing
