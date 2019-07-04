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
#include "FortranGrid.h"
#include "netcdf/File.h"
#include "progressbar.h"

namespace flood_processing {
namespace modules {

template<typename T>
inline void Downscaling<T>::coarse_to_fine(const Area& area, const nvector::View<T, 2>& coarse_flddph, nvector::Vector<T, 2>* fine_flddph) {
    const std::size_t fine_gridsize = area.size.x * area.size.y;
    const int* x = *area.grid;
    const int* y = &x[fine_gridsize];
    const float* dif = *area.flddif;
    T* p_fine_flddph = &fine_flddph->data()[0];
    for (std::size_t i = 0; i < fine_gridsize; ++i) {
        if (*x > 0) {
            const T& d = coarse_flddph(*y - 1, *x - 1);
            if (!std::isnan(d) && d > *dif) {
                *p_fine_flddph = d - *dif;
            } else {
                *p_fine_flddph = 0.0;
            }
        }
        ++p_fine_flddph;
        ++x;
        ++y;
        ++dif;
    }
}

template<typename T>
template<typename Function>
inline void Downscaling<T>::fine_to_med_dx_dy(
    std::size_t area_size_x, bool area_has_lonlat, T* p_fine_flddph, float const* lat, float const* lon, int dx, int dy, Function&& func) {
    const int offset = dy * area_size_x + dx;
    if (!area_has_lonlat || *(lon + offset) > -9999.0) {
        T dlon;
        T dlat;
        if (area_has_lonlat) {
            dlon = *(lon + offset) - *lon;
            dlat = *(lat + offset) - *lat;
            if (dlon > 180) {
                dlon -= 360.;
            } else if (dlon < -180) {
                dlon += 360.;
            }
        } else {
            dlon = fine_cell_size * dx;
            dlat = -fine_cell_size * dy;
        }
        const auto lat_ = std::min(static_cast<std::size_t>(std::max((90. - *lat) * inverse_target_cell_size + dlat * inverse_target_cell_size / 3., 0.0)),
                                   180 * inverse_target_cell_size - 1);
        const auto lon_ = std::min(static_cast<std::size_t>(std::max((180. + *lon) * inverse_target_cell_size + dlon * inverse_target_cell_size / 3., 0.0)),
                                   360 * inverse_target_cell_size - 1);
        func(lat_, lon_, *p_fine_flddph, *(p_fine_flddph + offset));
    }
}

template<typename T>
template<typename Function>
inline void Downscaling<T>::fine_to_med(Area& area, nvector::Vector<T, 2>* fine_flddph, Function&& func) {
    T* p_fine_flddph = &fine_flddph->data()[0];
    float* lon;
    float* lat;
    float lon_mem;
    float lat_mem;
    const bool area_has_lonlat = std::isnan(area.origin.lon);
    if (area_has_lonlat) {
        lon = *area.lonlat;
        lat = &lon[area.size.x * area.size.y];
    } else {
        lon = &lon_mem;
        lat = &lat_mem;
    }
    for (std::size_t y = 0; y < area.size.y; ++y) {
        for (std::size_t x = 0; x < area.size.x; ++x) {
            if (*p_fine_flddph >= 0.0) {
                if (area_has_lonlat) {
                    if (*lon > -9999.0) {
                        if (*lon < -180.0) {
                            *lon = *lon + 360.;
                        } else if (*lon > 180.0) {
                            *lon = *lon - 360.;
                        }
                    }
                } else {
                    *lon = area.origin.lon + fine_cell_size * (x + 0.5);
                    *lat = area.origin.lat - fine_cell_size * (y + 0.5);
                }
                if (inverse_target_cell_size < 200) {
                    if (x > 0) {
                        if (y > 0) {
                            fine_to_med_dx_dy(area.size.x, area_has_lonlat, p_fine_flddph, lat, lon, -1, -1, std::forward<Function>(func));
                        }
                        fine_to_med_dx_dy(area.size.x, area_has_lonlat, p_fine_flddph, lat, lon, -1, 0, std::forward<Function>(func));
                        if (y < area.size.y - 1) {
                            fine_to_med_dx_dy(area.size.x, area_has_lonlat, p_fine_flddph, lat, lon, -1, 1, std::forward<Function>(func));
                        }
                    }
                    if (y > 0) {
                        fine_to_med_dx_dy(area.size.x, area_has_lonlat, p_fine_flddph, lat, lon, 0, -1, std::forward<Function>(func));
                    }
                    fine_to_med_dx_dy(area.size.x, area_has_lonlat, p_fine_flddph, lat, lon, 0, 0, std::forward<Function>(func));
                    if (y < area.size.y - 1) {
                        fine_to_med_dx_dy(area.size.x, area_has_lonlat, p_fine_flddph, lat, lon, 0, 1, std::forward<Function>(func));
                    }
                    if (x < area.size.x - 1) {
                        if (y > 0) {
                            fine_to_med_dx_dy(area.size.x, area_has_lonlat, p_fine_flddph, lat, lon, 1, -1, std::forward<Function>(func));
                        }
                        fine_to_med_dx_dy(area.size.x, area_has_lonlat, p_fine_flddph, lat, lon, 1, 0, std::forward<Function>(func));
                        if (y < area.size.y - 1) {
                            fine_to_med_dx_dy(area.size.x, area_has_lonlat, p_fine_flddph, lat, lon, 1, 1, std::forward<Function>(func));
                        }
                    }
                } else {
                    const auto lat_ = std::min(static_cast<std::size_t>((90. - *lat) / fine_cell_size), 180 * inverse_target_cell_size - 1);
                    const auto lon_ = std::min(static_cast<std::size_t>((180. + *lon) / fine_cell_size), 360 * inverse_target_cell_size - 1);
                    func(lat_, lon_, *p_fine_flddph, *p_fine_flddph);
                }
            }
            if (area_has_lonlat) {
                ++lon;
                ++lat;
            }
            ++p_fine_flddph;
        }
    }
}

template<typename T>
Downscaling<T>::Downscaling(const settings::SettingsNode& settings) {
    map_path = settings["map_path"].as<std::string>();
    flddph_filename = settings["downscaled_flood_depth"]["filename"].as<std::string>();
    flddph_varname = settings["downscaled_flood_depth"]["varname"].as<std::string>();
    fldfrc_filename = settings["downscaled_flood_fraction"]["filename"].as<std::string>();
    fldfrc_varname = settings["downscaled_flood_fraction"]["varname"].as<std::string>();

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
        if (inverse_target_cell_size % 3 != 0) {
            throw std::runtime_error("inverse_target_cell_size must be 200 or a multiple of 3");
        }
        if (fine_lat_count % (inverse_target_cell_size / 3) != 0) {
            throw std::runtime_error("invalid inverse_target_cell_size");
        }
    }

    from_lat = settings["from_lat"].as<int>(-90);
    to_lat = settings["to_lat"].as<int>(90);
    from_lon = settings["from_lon"].as<int>(-180);
    to_lon = settings["to_lon"].as<int>(180);
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
    target_lon_count = (to_lon - from_lon) * inverse_target_cell_size;
    target_lat_count = (to_lat - from_lat) * inverse_target_cell_size;
}

template<typename T>
void Downscaling<T>::run(pipeline::Pipeline* p) {
    for (auto& area : areas) {
        const std::string name = map_path + "/" + area.name;
        area.grid.reset(new FortranGrid<int>(name + ".catmxy", 2 * area.size.x, area.size.y, 'r'));
        area.flddif.reset(new FortranGrid<float>(name + ".flddif", area.size.x, area.size.y, 'r'));
        if (std::isnan(area.origin.lon)) {
            area.lonlat.reset(new FortranGrid<float>(name + ".lonlat", 2 * area.size.x, area.size.y, 'r'));
        }
    }
    auto coarse_flddph = p->consume<nvector::View<T, 3>>(return_levels_name);
    const auto projection_times = p->consume<netCDF::DimVar<double>>(projection_times_name);
    netCDF::File flddph_file(flddph_filename, 'w');
    netCDF::File fldfrc_file(fldfrc_filename, 'w');
    netCDF::NcVar flddph_var = flddph_file.var<T>(flddph_varname, {flddph_file.dimvar(*projection_times), flddph_file.lat(target_lat_count, from_lat, to_lat),
                                                                   flddph_file.lon(target_lon_count, from_lon, to_lon)});
    netCDF::NcVar fldfrc_var = fldfrc_file.var<T>(fldfrc_varname, {fldfrc_file.dimvar(*projection_times), fldfrc_file.lat(target_lat_count, from_lat, to_lat),
                                                                   fldfrc_file.lon(target_lon_count, from_lon, to_lon)});
    downscale(*coarse_flddph, flddph_file, flddph_var, fldfrc_file, fldfrc_var);
    for (auto& area : areas) {
        area.grid.reset();
        area.flddif.reset();
        if (std::isnan(area.origin.lon)) {
            area.lonlat.reset();
        }
    }
}

template<typename T>
void Downscaling<T>::downscale(
    nvector::View<T, 3>& flddph, netCDF::File& result_flddph, netCDF::NcVar result_flddph_var, netCDF::File& result_fldfrc, netCDF::NcVar result_fldfrc_var) {
    progressbar::ProgressBar progress(flddph.template size<0>(), "Downscaling");
    flddph.template split<false, true, true>().foreach_element([&](std::size_t index, const nvector::View<T, 2>& coarse_flddph) {
        nvector::Vector<T, 2> flddph(0, target_lat_count, target_lon_count);
        nvector::Vector<T, 2> fldfrc(0, target_lat_count, target_lon_count);
        nvector::Vector<T, 2> fldnum(0, target_lat_count, target_lon_count);
#pragma omp parallel for default(shared) schedule(dynamic)
        for (std::size_t area_i = 0; area_i < areas.size(); ++area_i) {
            auto& area = areas[area_i];
            nvector::Vector<T, 2> fine_flddph(-9999.0, area.size.x, area.size.y);
            coarse_to_fine(area, coarse_flddph, &fine_flddph);
            fine_to_med(area, &fine_flddph, [&](std::size_t med_lat, std::size_t med_lon, T cell_dph, T dcell_dph) {
                if (med_lat >= inverse_target_cell_size * (90 - to_lat) && med_lat < inverse_target_cell_size * (90 - from_lat)
                    && med_lon >= inverse_target_cell_size * (180 + from_lon) && med_lon < inverse_target_cell_size * (180 + to_lon)) {
                    const auto target_lat = static_cast<int>(med_lat) - inverse_target_cell_size * (90 - to_lat);
                    const auto target_lon = static_cast<int>(med_lon) - inverse_target_cell_size * (180 + from_lon);
                    if (dcell_dph > 0) {
                        T& tmp = flddph(target_lat, target_lon);
                        tmp = std::max(tmp, cell_dph);
                        ++fldfrc(target_lat, target_lon);
                    }
                    ++fldnum(target_lat, target_lon);
                }
            });
        }
        T* num = &fldnum.data()[0];
        T* dph = &flddph.data()[0];
        T* frc = &fldfrc.data()[0];
        for (std::size_t i = 0; i < target_lon_count * target_lat_count; ++i) {
            if (*num > 0) {
                *frc /= *num;
            } else {
                *dph = std::numeric_limits<T>::quiet_NaN();
                *frc = std::numeric_limits<T>::quiet_NaN();
            }
            ++num;
            ++frc;
            ++dph;
        }
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
