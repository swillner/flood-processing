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

#include "modules/sum_per_region.h"
#include "nvector.h"

namespace flood_processing {
namespace modules {

template<typename T>
void SumPerRegion<T>::run(pipeline::Pipeline* p) {
    const auto region_index_raster_raw = p->consume<nvector::Vector<T, 3>>("region_index_raster");
    if (region_index_raster_raw->template size<0>() != 1) {
        throw std::runtime_error("Only 2d input grids supported");
    }
    const auto region_index_raster = region_index_raster_raw->template split<nvector::Split<false, true, true>>().at(0);
    const auto regions = p->consume<std::vector<std::string>>("regions");
    const auto regions_count = regions->size();
    if (filename.empty()) {
        const auto data = p->consume<nvector::Vector<T, 3>>(inputname);
        if (data->template size<1>() != region_index_raster.template size<0>() || data->template size<2>() != region_index_raster.template size<1>()) {
            throw std::runtime_error("grid sizes differ");
        }
        const auto time_count = data->template size<0>();
        auto output = std::make_shared<nvector::Vector<T, 2>>(0, time_count, regions_count);
        data->template split<false, true, true>().foreach_parallel([&](std::size_t index, nvector::View<T, 2>& grid) {
            nvector::foreach_view(nvector::collect(grid, region_index_raster), [&](std::size_t lat, std::size_t lon, T d, T region_index_l) {
                (void)lat;
                (void)lon;
                if (d > 0 && region_index_l >= 0) {
                    (*output)(index, region_index_l) += d;
                }
                return true;
            });
        });
        p->provide<nvector::Vector<T, 2>>(outputname, output);
    } else {
        netCDF::File file(filename, 'r');
        netCDF::NcVar var = file.var(varname);
        const auto time_count = file.size<0>(var);
        auto output = std::make_shared<nvector::Vector<T, 2>>(0, time_count, regions_count);
        const std::size_t lat_count = file.size<1>(var);
        const std::size_t lon_count = file.size<2>(var);
        if (lat_count != region_index_raster.template size<0>() || lon_count != region_index_raster.template size<1>()) {
            throw std::runtime_error("grid sizes differ");
        }
        nvector::Vector<T, 2> grid(0, lat_count, lon_count);
        for (std::size_t index = 0; index < time_count; ++index) {
            var.getVar({index, 0, 0}, {1, lat_count, lon_count}, &grid.data()[0]);
            nvector::foreach_view(nvector::collect(grid, region_index_raster), [&](std::size_t lat, std::size_t lon, T d, T region_index_l) {
                (void)lat;
                (void)lon;
                if (d > 0 && region_index_l >= 0) {
                    (*output)(index, region_index_l) += d;
                }
                return true;
            });
        }
        p->provide<nvector::Vector<T, 2>>(outputname, output);
    }
}

template class SumPerRegion<float>;
template class SumPerRegion<double>;

}  // namespace modules
}  // namespace flood_processing
