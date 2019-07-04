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

#ifndef FLOOD_PROCESSING_RASTERIZE_H
#define FLOOD_PROCESSING_RASTERIZE_H

#include <string>
#include <unordered_map>
#include <vector>
#include "nvector.h"
#include "pipeline.h"
#include "settingsnode.h"

namespace flood_processing {
namespace modules {

template<typename T>
class Rasterization : public pipeline::Module {
  protected:
    std::string resolution_mask_name;
    std::string shapefilename;
    std::vector<std::string> layernames;
    std::string fieldname;
    std::size_t max_advance;
    std::size_t adjust_scale;
    T invalid_value;
    bool adjust_max;
    std::size_t xres;
    std::size_t yres;
    void advance(nvector::View<T, 2>& result, std::size_t max_advance);
    template<typename Function>
    void rasterize(nvector::View<T, 2>& result, Function&& func);
    void combine_fine(const nvector::View<T, 2>& fine_raster, nvector::View<T, 2>& raster);

  public:
    explicit Rasterization(const settings::SettingsNode& settings);
    void run(pipeline::Pipeline* p) override;

    inline pipeline::ModuleDescription describe() override {
        if (resolution_mask_name.empty()) {
            return pipeline::ModuleDescription{"rasterization", {}, {"raster"}};
        }
        return pipeline::ModuleDescription{"rasterization", {resolution_mask_name}, {"raster"}};
    }
};

template<typename T>
class RegionIndexRasterization : public Rasterization<T> {
  protected:
    using Rasterization<T>::resolution_mask_name;
    using Rasterization<T>::rasterize;
    using Rasterization<T>::max_advance;
    using Rasterization<T>::adjust_scale;
    using Rasterization<T>::advance;
    using Rasterization<T>::combine_fine;
    using Rasterization<T>::invalid_value;
    using Rasterization<T>::xres;
    using Rasterization<T>::yres;
    std::vector<std::string> fieldnames;
    std::unordered_map<std::string, std::string> correction_map;
    std::unordered_map<std::string, std::string> iso3_to_iso2;

  public:
    explicit RegionIndexRasterization(const settings::SettingsNode& settings);
    void run(pipeline::Pipeline* p) override;

    inline pipeline::ModuleDescription describe() override {
        if (resolution_mask_name.empty()) {
            return pipeline::ModuleDescription{"region_index_rasterization", {"regions"}, {"region_index_raster"}};
        }
        return pipeline::ModuleDescription{"region_index_rasterization", {"regions", resolution_mask_name}, {"region_index_raster"}};
    }
};

}  // namespace modules
}  // namespace flood_processing

#endif
