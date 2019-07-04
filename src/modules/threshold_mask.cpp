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
#include <cmath>
#include "nvector.h"

namespace flood_processing {
namespace modules {

template<typename T>
ThresholdMask<T>::ThresholdMask(const settings::SettingsNode& settings) {
    const auto threshold_string = settings["threshold"].as<std::string>();
    std::cout << "threshold_string: " << threshold_string << std::endl;
    std::istringstream ss(threshold_string);
    ss >> std::noskipws >> threshold;
    std::cout << "threshold: " << threshold << std::endl;
    use_threshold_raster = !ss.eof() || ss.fail();
    std::cout << "use_threshold_raster: " << use_threshold_raster << std::endl;
    if (use_threshold_raster) {
        threshold_raster_name = threshold_string;
    }
    less_than = settings.has("le_than") || settings.has("less_than");
    or_equal = settings.has("le_than") || settings.has("ge_than");
    if (less_than) {
        if (or_equal) {
            compare_to_name = settings["le_than"].as<std::string>();
        } else {
            compare_to_name = settings["less_than"].as<std::string>();
        }
    } else {
        if (or_equal) {
            compare_to_name = settings["ge_than"].as<std::string>();
        } else {
            compare_to_name = settings["greater_than"].as<std::string>();
        }
    }
    apply_to_name = settings["apply_to"].as<std::string>();
    output_name = settings["output"].as<std::string>();
    mask_value = settings["mask_value"].as<T>();
}

template<typename T>
nvector::View<T, 3> duplicate_time_if_necessary(std::size_t size, nvector::Vector<T, 3>& v) {
    nvector::Slice slice = v.template slice<0>();
    if (slice.size == 1 && slice.begin == 0) {
        slice.stride = 0;
        slice.size = size;
    } else if (v.template size<0>() != size) {
        throw std::runtime_error("Cannot duplicate time dimensions");
    }
    return nvector::View<T, 3>(std::begin(v.data()), {slice, v.template slice<1>(), v.template slice<2>()});
}

template<typename T>
void ThresholdMask<T>::run(pipeline::Pipeline* p) {
    const auto apply_to = p->consume<nvector::Vector<T, 3>>(apply_to_name);
    const auto size = apply_to->template size<0>();
    const auto lat_count = apply_to->template size<1>();
    const auto lon_count = apply_to->template size<2>();
    auto output = std::make_shared<nvector::Vector<T, 3>>(std::numeric_limits<T>::quiet_NaN(), size, lat_count, lon_count);

    const auto compare_to_raw = p->consume<nvector::Vector<T, 3>>(compare_to_name);  // need to keep this if it's the last shared_ptr
    const auto compare_to = duplicate_time_if_necessary(size, *compare_to_raw);

    if (use_threshold_raster) {
        const auto threshold_raster_raw = p->consume<nvector::Vector<T, 3>>(threshold_raster_name);  // need to keep this if it's the last shared_ptr
        const auto threshold_raster = duplicate_time_if_necessary(size, *threshold_raster_raw);
        nvector::foreach_view_parallel(nvector::collect(*apply_to, compare_to, threshold_raster, *output), [&](std::size_t t, std::size_t lat, std::size_t lon,
                                                                                                               T a, T c, T this_threshold, T& out) {
            if (!std::isnan(a) && !std::isnan(c) && !std::isnan(this_threshold)) {
                if ((c < this_threshold && less_than) || (c > this_threshold && !less_than) || (c == this_threshold && or_equal)) {
                    out = a;
                } else {
                    out = mask_value;
                }
            }
            return true;
        });
    } else {
        nvector::foreach_view_parallel(nvector::collect(*apply_to, compare_to, *output),
                                       [&](std::size_t t, std::size_t lat, std::size_t lon, T a, T c, T& out) {
                                           if (!std::isnan(a) && !std::isnan(c)) {
                                               if ((c < threshold && less_than) || (c > threshold && !less_than) || (c == threshold && or_equal)) {
                                                   out = a;
                                               } else {
                                                   out = mask_value;
                                               }
                                           }
                                           return true;
                                       });
    }
    p->provide<nvector::Vector<T, 3>>(output_name, output);
}

template class ThresholdMask<float>;
template class ThresholdMask<double>;

}  // namespace modules
}  // namespace flood_processing
