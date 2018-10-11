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
#include <memory>
#include <tuple>
#include "FortranGrid.h"
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
  protected:
    struct Area {
        const char* name;
        struct {
            float lon;
            float lat;
        } origin;
        struct {
            std::size_t x;
            std::size_t y;
        } size;
        std::unique_ptr<FortranGrid<int>> grid;
        std::unique_ptr<FortranGrid<float>> flddif;
        std::unique_ptr<FortranGrid<float>> lonlat;
    };
    std::array<Area, 13> areas = std::array<Area, 13>{
        Area{"sa1", {-85, 15}, {11000, 15000}, nullptr, nullptr, nullptr},
        Area{"ca1", {-120, 40}, {12000, 7000}, nullptr, nullptr, nullptr},
        Area{"na1", {-130, 60}, {16000, 7000}, nullptr, nullptr, nullptr},
        Area{"af1", {5, 35}, {11000, 14000}, nullptr, nullptr, nullptr},
        Area{"eu1", {-20, 60}, {8000, 12000}, nullptr, nullptr, nullptr},
        Area{"eu2", {5, 60}, {13000, 8000}, nullptr, nullptr, nullptr},
        Area{"as1", {55, 60}, {9000, 11000}, nullptr, nullptr, nullptr},
        Area{"as2", {90, 60}, {12000, 8000}, nullptr, nullptr, nullptr},
        Area{"as3", {90, 35}, {13000, 10000}, nullptr, nullptr, nullptr},
        Area{"oc1", {110, -10}, {14000, 8000}, nullptr, nullptr, nullptr},
        Area{"as4", {std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()}, {10256, 10263}, nullptr, nullptr, nullptr},
        Area{"na2", {std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()}, {9102, 8384}, nullptr, nullptr, nullptr},
        Area{"eu3", {std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()}, {6888, 7385}, nullptr, nullptr, nullptr}};
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
    int from_lat, to_lat;
    int from_lon, to_lon;

    inline void coarse_to_fine(const Area& area, const nvector::View<T, 2>& coarse_flddph, nvector::Vector<T, 2>* fine_flddph);
    template<typename Function>
    inline void fine_to_med_dx_dy(std::size_t area_size_x, bool area_has_lonlat, T* p_fine_flddph, float* lat, float* lon, int dx, int dy, Function&& func);
    template<typename Function>
    inline void fine_to_med(Area& area, nvector::Vector<T, 2>* fine_flddph, Function&& func);
    inline void downscale(nvector::View<T, 3>& coarse_flddph,
                          netCDF::File& result_flddph,
                          netCDF::NcVar result_flddph_var,
                          netCDF::File& result_fldfrc,
                          netCDF::NcVar result_fldfrc_var);

  public:
    Downscaling(const settings::SettingsNode& settings);
    void run(pipeline::Pipeline* p) override;

    inline pipeline::ModuleDescription describe() override {
        return pipeline::ModuleDescription{"downscaling", {"projection_times", "return_levels_thresholded"}, {}};
    }
};

}  // namespace modules
}  // namespace flood_processing

#endif
