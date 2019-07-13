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

#ifndef FLOOD_PROCESSING_DOWNSCALING_H
#define FLOOD_PROCESSING_DOWNSCALING_H

#include <array>
#include <cstdint>
#include <memory>
#include <tuple>
#include "cudatools.h"
#include "nvector.h"
#include "pipeline.h"
#include "settingsnode.h"

namespace netCDF {
class File;
class NcVar;
}  // namespace netCDF

namespace flood_processing {
namespace modules {

template<typename T>
class Downscaling : public pipeline::Module {
  public:
    struct Area {
        const char* name = "\0";
        struct {
            float lon;
            float lat;
        } origin;
        struct {
            std::size_t x;
            std::size_t y;
        } size;
        nvector::Vector<std::int16_t, 2, cudatools::vector<std::int16_t, true>> gridx;
        nvector::Vector<std::int16_t, 2, cudatools::vector<std::int16_t, true>> gridy;
        nvector::Vector<float, 2, cudatools::vector<float, true>> flddif;
    };

  protected:
    std::array<Area, 14> areas{
        Area{"sa1", {-85, 15}, {11000, 15000}, {}, {}}, Area{"ca1", {-120, 40}, {12000, 7000}, {}, {}}, Area{"na1", {-130, 60}, {16000, 7000}, {}, {}},
        Area{"af1", {5, 35}, {11000, 14000}, {}, {}},   Area{"eu1", {-20, 60}, {8000, 12000}, {}, {}},  Area{"eu2", {5, 60}, {13000, 8000}, {}, {}},
        Area{"as1", {55, 60}, {9000, 11000}, {}, {}},   Area{"as2", {90, 60}, {12000, 8000}, {}, {}},   Area{"as3", {90, 35}, {13000, 10000}, {}, {}},
        Area{"oc1", {110, -10}, {14000, 8000}, {}, {}}, Area{"na2", {-170, 75}, {23000, 5000}, {}, {}}, Area{"eu3", {0, 80}, {14000, 7000}, {}, {}},
        Area{"si1", {55, 80}, {12000, 7000}, {}, {}},   Area{"si2", {100, 75}, {19000, 5000}, {}, {}},
    };
    std::size_t target_lon_count;
    std::size_t target_lat_count;
    const std::size_t coarse_lon_count = 1440;
    const std::size_t coarse_lat_count = 720;
    const std::size_t fine_inverse_cell_size = 200;
    const float fine_cell_size = 1. / static_cast<float>(fine_inverse_cell_size);
    const std::size_t fine_lon_count = 360 * fine_inverse_cell_size;
    const std::size_t fine_lat_count = 180 * fine_inverse_cell_size;
    std::string map_path;
    std::string flddph_filename;
    std::string flddph_varname;
    std::string fldfrc_filename;
    std::string fldfrc_varname;
    std::size_t inverse_target_cell_size;
    std::string projection_times_name;
    std::string return_levels_name;
    float from_lat, to_lat;
    float from_lon, to_lon;

    nvector::Vector<T, 2> coarse_to_fine(const Area& area, const nvector::View<T, 2>& coarse_flddph) const;
    template<typename Function>
    constexpr void fine_to_med_dx_dy(std::size_t area_size_x, T const* p_fine_flddph, float lat, float lon, int dx, int dy, Function&& func) const;
    template<typename Function>
    void fine_to_med(const Area& area, const nvector::Vector<T, 2>& fine_flddph, Function&& func);

  public:
    void coarse_to_fine_gpu(const Area& area, const nvector::View<T, 2, T*>& coarse_flddph, nvector::View<T, 2, T*>& result) const;
    void downscale(const nvector::View<T, 3>& timed_flddph,
                   netCDF::File& result_flddph,
                   netCDF::NcVar result_flddph_var,
                   netCDF::File& result_fldfrc,
                   netCDF::NcVar result_fldfrc_var);

    explicit Downscaling(const settings::SettingsNode& settings);
    void run(pipeline::Pipeline* p) override;

    inline pipeline::ModuleDescription describe() override {
        return pipeline::ModuleDescription{"downscaling", {projection_times_name, return_levels_name}, {}};
    }
};

}  // namespace modules
}  // namespace flood_processing

#endif
