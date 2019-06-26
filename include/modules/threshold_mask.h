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

#ifndef FLOOD_PROCESSING_THRESHOLD_MASK_H
#define FLOOD_PROCESSING_THRESHOLD_MASK_H

#include "pipeline.h"
#include "settingsnode.h"

namespace flood_processing {
namespace modules {

template<typename T>
class ThresholdMask : public pipeline::Module {
  protected:
    bool less_than;
    bool or_equal;
    bool use_raster;
    T mask_value;
    T threshold;
    std::string apply_to_name;
    std::string compare_to_name;
    std::string output_name;
    std::string raster_name;

  public:
    ThresholdMask(const settings::SettingsNode& settings) {
        use_raster = settings.has("raster");
        if (use_raster) {
            raster_name = settings["raster"].as<std::string>();
        } else {
            threshold = settings["threshold"].as<T>();
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
    void run(pipeline::Pipeline* p) override;

    inline pipeline::ModuleDescription describe() override {
        if (use_raster) {
            return pipeline::ModuleDescription{"threshold_mask", {raster_name, compare_to_name, apply_to_name}, {output_name}};
        } else {
            return pipeline::ModuleDescription{"threshold_mask", {compare_to_name, apply_to_name}, {output_name}};
        }
    }
};

}  // namespace modules
}  // namespace flood_processing

#endif
