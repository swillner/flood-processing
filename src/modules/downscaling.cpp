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

#include "modules/downscaling.h"
#include <cmath>
#include <fstream>
#include "cudatools.h"
#include "netcdf/File.h"
#include "progressbar.h"

namespace flood_processing {
namespace modules {

template<typename T>
nvector::Vector<T, 2> Downscaling<T>::coarse_to_fine(const Area& area, const nvector::View<T, 2>& coarse_flddph) const {
    nvector::Vector<T, 2> result(std::numeric_limits<T>::quiet_NaN(), area.size.y, area.size.x);
    nvector::foreach_aligned_parallel(nvector::collect(area.gridx, area.gridy, area.flddif, result),
                                      [&](std::size_t i, std::int16_t x, std::int16_t y, float dif, T& out) {
                                          (void)i;
                                          if (x > 0) {
                                              const auto d = coarse_flddph(y - 1, x - 1);
                                              if (std::isnan(d)) {
                                                  out = 0.0;
                                              } else {
                                                  out = std::max<T>(0.0, d - dif);
                                              }
                                          }
                                      });
    return result;
}

template<typename T>
void Downscaling<T>::coarse_to_fine_gpu(const Area& area, const nvector::View<T, 2, T*>& coarse_flddph, nvector::View<T, 2, T*>& result) const {
    const auto* coarse_flddph_ptr = &coarse_flddph[0];
    const auto size_x = coarse_flddph.template size<1>();
    nvector::foreach_aligned_gpu(nvector::collect(area.gridx, area.gridy, area.flddif, result),
                                 [=] CUDA_DEVICE(std::size_t i, std::int16_t x, std::int16_t y, float dif, T& out) {
                                     (void)i;
                                     if (x > 0) {
                                         const auto d = coarse_flddph_ptr[(y - 1) * size_x + x - 1];
                                         if (std::isnan(d)) {
                                             out = 0.0;
                                         } else {
                                             if (d > dif) {
                                                 out = d - dif;
                                             } else {
                                                 out = 0.0;
                                             }
                                         }
                                     } else {
                                         out = -1;
                                     }
                                 });
}

template<typename T>
template<typename Function>
constexpr void Downscaling<T>::fine_to_med_dx_dy(
    std::size_t area_size_x, T const* p_fine_flddph, double lat, double lon, int dx, int dy, Function&& func) const {
    const int offset = dy * area_size_x + dx;
    const T dlon = fine_cell_size * dx;
    const T dlat = -fine_cell_size * dy;
    const auto lat_ = std::min(static_cast<std::size_t>(std::max((90. - lat) * inverse_target_cell_size + dlat * inverse_target_cell_size / 3., 0.0)),
                               180 * inverse_target_cell_size - 1);
    const auto lon_ = std::min(static_cast<std::size_t>(std::max((180. + lon) * inverse_target_cell_size + dlon * inverse_target_cell_size / 3., 0.0)),
                               360 * inverse_target_cell_size - 1);
    func(lat_, lon_, *p_fine_flddph, *(p_fine_flddph + offset));
}

template<typename T>
template<typename Function>
void Downscaling<T>::fine_to_med(const Area& area, const nvector::Vector<T, 2>& fine_flddph, Function&& func) {
    T const* p_fine_flddph = &fine_flddph.data()[0];
    for (std::size_t y = 0; y < area.size.y; ++y) {
        for (std::size_t x = 0; x < area.size.x; ++x) {
            if (!std::isnan(*p_fine_flddph)) {
                const double lon = area.origin.lon + fine_cell_size * (x + 0.5);
                const double lat = area.origin.lat - fine_cell_size * (y + 0.5);
                if (inverse_target_cell_size < 200) {
                    if (x > 0) {
                        if (y > 0) {
                            fine_to_med_dx_dy(area.size.x, p_fine_flddph, lat, lon, -1, -1, std::forward<Function>(func));
                        }
                        fine_to_med_dx_dy(area.size.x, p_fine_flddph, lat, lon, -1, 0, std::forward<Function>(func));
                        if (y < area.size.y - 1) {
                            fine_to_med_dx_dy(area.size.x, p_fine_flddph, lat, lon, -1, 1, std::forward<Function>(func));
                        }
                    }
                    if (y > 0) {
                        fine_to_med_dx_dy(area.size.x, p_fine_flddph, lat, lon, 0, -1, std::forward<Function>(func));
                    }
                    fine_to_med_dx_dy(area.size.x, p_fine_flddph, lat, lon, 0, 0, std::forward<Function>(func));
                    if (y < area.size.y - 1) {
                        fine_to_med_dx_dy(area.size.x, p_fine_flddph, lat, lon, 0, 1, std::forward<Function>(func));
                    }
                    if (x < area.size.x - 1) {
                        if (y > 0) {
                            fine_to_med_dx_dy(area.size.x, p_fine_flddph, lat, lon, 1, -1, std::forward<Function>(func));
                        }
                        fine_to_med_dx_dy(area.size.x, p_fine_flddph, lat, lon, 1, 0, std::forward<Function>(func));
                        if (y < area.size.y - 1) {
                            fine_to_med_dx_dy(area.size.x, p_fine_flddph, lat, lon, 1, 1, std::forward<Function>(func));
                        }
                    }
                } else {
                    const auto lat_ = std::min(static_cast<std::size_t>((90. - lat) / fine_cell_size), 180 * inverse_target_cell_size - 1);
                    const auto lon_ = std::min(static_cast<std::size_t>((180. + lon) / fine_cell_size), 360 * inverse_target_cell_size - 1);
                    func(lat_, lon_, *p_fine_flddph, *p_fine_flddph);
                }
            }
            ++p_fine_flddph;
        }
    }
}

template<typename T>
Downscaling<T>::Downscaling(const settings::SettingsNode& settings) {
    map_path = settings["map_path"].as<std::string>();
    flddph_filename = settings["downscaled_flood_depth"]["filename"].as<std::string>();
    flddph_varname = settings["downscaled_flood_depth"]["variable"].as<std::string>();
    fldfrc_filename = settings["downscaled_flood_fraction"]["filename"].as<std::string>();
    fldfrc_varname = settings["downscaled_flood_fraction"]["variable"].as<std::string>();

    if (settings.has("input")) {
        projection_times_name = settings["input"]["times"].as<std::string>("projection_times");
        return_levels_name = settings["input"]["return_levels"].as<std::string>("return_levels");
    } else {
        projection_times_name = "projection_times";
        return_levels_name = "return_levels";
    }

    inverse_target_cell_size = settings["inverse_target_cell_size"].as<std::size_t>();
    if (inverse_target_cell_size > 200 || inverse_target_cell_size == 0) {
        throw std::runtime_error("inverse_target_cell_size must be less than or equal 200");
    }
    if (inverse_target_cell_size < 200) {
        if (inverse_target_cell_size % 3 == 0) {
            if (fine_lat_count % (inverse_target_cell_size / 3) != 0) {
                throw std::runtime_error("invalid inverse_target_cell_size");
            }
        } else {
            if (fine_lat_count % inverse_target_cell_size != 0) {
                throw std::runtime_error("invalid inverse_target_cell_size");
            }
        }
    }

    from_lat = settings["from_lat"].as<double>(-90);
    to_lat = settings["to_lat"].as<double>(90);
    from_lon = settings["from_lon"].as<double>(-180);
    to_lon = settings["to_lon"].as<double>(180);
    if (from_lat < -90) {
        throw std::runtime_error("invalid from_lat");
    }
    if (to_lat > 90 || from_lat >= to_lat) {
        throw std::runtime_error("invalid to_lat");
    }
    if (from_lon < -180) {
        throw std::runtime_error("invalid from_lon");
    }
    if (to_lon > 180 || from_lon >= to_lon) {
        throw std::runtime_error("invalid to_lon");
    }

    // calculate these before adjusting lat/lon boundaries to avoid rounding issues
    target_lon_count =
        360 * inverse_target_cell_size - std::floor((180 - to_lon) * inverse_target_cell_size) - std::floor((from_lon - -180) * inverse_target_cell_size);
    target_lat_count =
        180 * inverse_target_cell_size - std::floor((90 - to_lat) * inverse_target_cell_size) - std::floor((from_lat - -90) * inverse_target_cell_size);

    from_lat = -90 + std::floor((from_lat - -90) * inverse_target_cell_size) / inverse_target_cell_size;
    from_lon = -180 + std::floor((from_lon - -180) * inverse_target_cell_size) / inverse_target_cell_size;
    to_lat = 90 - std::floor((90 - to_lat) * inverse_target_cell_size) / inverse_target_cell_size;
    to_lon = 180 - std::floor((180 - to_lon) * inverse_target_cell_size) / inverse_target_cell_size;
}

template<typename T>
void Downscaling<T>::run(pipeline::Pipeline* p) {
    for (auto& area : areas) {
        const auto name = map_path + "/" + area.name;
        {
            area.gridx.resize(0, area.size.y, area.size.x);
            area.gridy.resize(0, area.size.y, area.size.x);
            std::ifstream grid(name + ".catmxy", std::ios::in | std::ios::binary);
            area.gridx.data().read(grid);
            area.gridy.data().read(grid);
        }
        {
            area.flddif.resize(0, area.size.y, area.size.x);
            std::ifstream flddif(name + ".flddif", std::ios::in | std::ios::binary);
            area.flddif.data().read(flddif);
        }
    }
    auto chunk_lat_count = target_lat_count;
    auto chunk_lon_count = target_lon_count;
    while (chunk_lat_count * chunk_lon_count >= 1000000000) {
        chunk_lat_count /= 2;
        chunk_lon_count /= 2;
    }
    std::vector<std::size_t> chunks = {1, chunk_lat_count, chunk_lon_count};
    auto coarse_flddph = p->consume<nvector::View<T, 3>>(return_levels_name);
    const auto projection_times = p->consume<netCDF::DimVar<double>>(projection_times_name);

    netCDF::File flddph_file(flddph_filename, 'w');
    netCDF::NcVar flddph_var = flddph_file.var<T>(flddph_varname, {flddph_file.dimvar(*projection_times), flddph_file.lat(target_lat_count, from_lat, to_lat),
                                                                   flddph_file.lon(target_lon_count, from_lon, to_lon)});
    flddph_var.putAtt("grid_mapping", "crs");
    flddph_var.putAtt("units", "m");
    flddph_var.putAtt("long_name", "flood depth");
    flddph_var.setChunking(netCDF::NcVar::nc_CHUNKED, chunks);
    nc_del_att(flddph_file.getId(), flddph_file.var("time").getId(), "bounds");

    netCDF::File fldfrc_file(fldfrc_filename, 'w');
    netCDF::NcVar fldfrc_var = fldfrc_file.var<T>(fldfrc_varname, {fldfrc_file.dimvar(*projection_times), fldfrc_file.lat(target_lat_count, from_lat, to_lat),
                                                                   fldfrc_file.lon(target_lon_count, from_lon, to_lon)});
    fldfrc_var.putAtt("grid_mapping", "crs");
    fldfrc_var.putAtt("units", "1");
    fldfrc_var.putAtt("long_name", "flood fraction");
    fldfrc_var.setChunking(netCDF::NcVar::nc_CHUNKED, chunks);
    nc_del_att(fldfrc_file.getId(), fldfrc_file.var("time").getId(), "bounds");

    downscale(*coarse_flddph, flddph_file, flddph_var, fldfrc_file, fldfrc_var);

    for (auto& area : areas) {
        area.gridx.resize(0, 0, 0);
        area.gridy.resize(0, 0, 0);
        area.flddif.resize(0, 0, 0);
    }
}

template<typename T>
void Downscaling<T>::downscale_remapping(const Area& area,
                                         const nvector::Vector<T, 2>& fine_flddph,
                                         nvector::Vector<T, 2, cudatools::vector<T>>& flddph,
                                         nvector::Vector<T, 2, cudatools::vector<T>>& fldfrc,
                                         nvector::Vector<T, 2, cudatools::vector<T>>& fldnum) {
    fine_to_med(area, fine_flddph, [&](std::size_t med_lat, std::size_t med_lon, T cell_dph, T dcell_dph) {
        if (med_lat >= inverse_target_cell_size * (90 - to_lat) && med_lat < inverse_target_cell_size * (90 - from_lat)
            && med_lon >= inverse_target_cell_size * (180 + from_lon) && med_lon < inverse_target_cell_size * (180 + to_lon)) {
            const auto target_lat = static_cast<int>(med_lat) - inverse_target_cell_size * (90 - to_lat);
            const auto target_lon = static_cast<int>(med_lon) - inverse_target_cell_size * (180 + from_lon);
            if (target_lat >= flddph.template size<0>() || target_lon >= flddph.template size<1>()) {
                return;
            }
            if (dcell_dph > 0) {
                T& tmp = flddph(target_lat, target_lon);
                tmp = std::max(tmp, cell_dph);
                ++fldfrc(target_lat, target_lon);
            }
            ++fldnum(target_lat, target_lon);
        }
    });
}

template<typename T>
void Downscaling<T>::downscale(const nvector::View<T, 3>& timed_flddph,
                               netCDF::File& result_flddph,
                               netCDF::NcVar result_flddph_var,
                               netCDF::File& result_fldfrc,
                               netCDF::NcVar result_fldfrc_var) {
    // nvector::Vector<T, 3, cudatools::vector<T, true>> coarse_flddph_gpu(0, timed_flddph.slices());
    // coarse_flddph_gpu.data().set(&timed_flddph[0]);
    progressbar::ProgressBar progress(timed_flddph.template size<0>(), "Downscaling");
    nvector::Vector<T, 2, cudatools::vector<T>> flddph(0, target_lat_count, target_lon_count);
    nvector::Vector<T, 2, cudatools::vector<T>> fldfrc(0, target_lat_count, target_lon_count);
    nvector::Vector<T, 2, cudatools::vector<T>> fldnum(0, target_lat_count, target_lon_count);
    nvector::foreach_split<nvector::Split<false, true, true>>(nvector::collect(timed_flddph), [&](std::size_t index, const nvector::View<T, 2>& coarse_flddph) {
        flddph.reset(0);
        fldfrc.reset(0);
        fldnum.reset(0);
        if (inverse_target_cell_size % 3 == 0) {
#pragma omp parallel for default(shared) schedule(dynamic)
            for (std::size_t area_i = 0; area_i < areas.size(); ++area_i) {
                const auto& area = areas[area_i];
                const auto fine_flddph = coarse_to_fine(area, coarse_flddph);
                downscale_remapping(area, fine_flddph, flddph, fldfrc, fldnum);
            }
        } else {
            for (auto& area : areas) {
                // nvector::Vector<T, 2, cudatools::vector<T, true>> fine_flddph(std::numeric_limits<T>::quiet_NaN(), area.size.y, area.size.x);
                // coarse_to_fine_gpu(area, coarse_flddph_gpu.template split<nvector::Split<false, true, true>>()[index], fine_flddph);
                // coarse_to_fine_gpu(area, coarse_flddph, fine_flddph);
                const auto fine_flddph = coarse_to_fine(area, coarse_flddph);
                T* fine_flddph_ptr = &fine_flddph[0];
                const auto size_x = area.size.x;
                const auto size_y = area.size.y;
                const auto origin_lat = area.origin.lat;
                const auto origin_lon = area.origin.lon;
                const auto fine_inverse_cell_size_l = fine_inverse_cell_size;
                const auto inverse_target_cell_size_l = inverse_target_cell_size;
                const auto target_lon_count_l = target_lon_count;
                const auto from_lat_l = from_lat;
                const auto from_lon_l = from_lon;
                const auto to_lat_l = to_lat;
                const auto to_lon_l = to_lon;
                nvector::foreach_aligned_gpu(nvector::collect(flddph, fldfrc, fldnum), [=] CUDA_DEVICE(std::size_t i, T & dph, T & frc, T & num) {
                    const auto lat_index = i / target_lon_count_l;
                    const auto lon_index = i % target_lon_count_l;
                    const auto lat = to_lat_l - static_cast<double>(lat_index) / static_cast<double>(inverse_target_cell_size_l);
                    if (lat > origin_lat || lat < from_lat_l) {
                        return;
                    }
                    const auto lon = from_lon_l + static_cast<double>(lon_index) / static_cast<double>(inverse_target_cell_size_l);
                    if (lon < origin_lon || lon > to_lon_l) {
                        return;
                    }
                    const auto start_y = static_cast<std::size_t>((origin_lat - lat) * fine_inverse_cell_size_l);
                    if (start_y >= size_y) {
                        return;
                    }
                    const auto start_x = static_cast<std::size_t>((lon - origin_lon) * fine_inverse_cell_size_l);
                    if (start_x >= size_x) {
                        return;
                    }
                    for (std::size_t y = start_y; y < start_y + (fine_inverse_cell_size_l / inverse_target_cell_size_l); ++y) {
                        for (std::size_t x = start_x; x < start_x + (fine_inverse_cell_size_l / inverse_target_cell_size_l); ++x) {
                            const auto fine_dph = fine_flddph_ptr[y * size_x + x];
                            if (fine_dph >= 0) {
                                if (fine_dph > 0) {
                                    if (fine_dph > dph) {
                                        dph = fine_dph;
                                    }
                                    ++frc;
                                }
                                ++num;
                            }
                        }
                    }
                });
            }
        }
        const auto nan = std::numeric_limits<T>::quiet_NaN();
        nvector::foreach_aligned_gpu(nvector::collect(fldnum, flddph, fldfrc), [=] CUDA_DEVICE(std::size_t i, T num, T & dph, T & frc) {
            (void)i;
            if (num > 0) {
                frc /= num;
            } else {
                dph = nan;
                frc = nan;
            }
        });
        result_flddph.set<T, 2>(result_flddph_var, flddph, index);
        result_fldfrc.set<T, 2>(result_fldfrc_var, fldfrc, index);
        ++progress;
        return true;
    });
}

template class Downscaling<float>;
template class Downscaling<double>;

}  // namespace modules
}  // namespace flood_processing
