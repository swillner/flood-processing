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

#include <sstream>
#include "pipeline.h"
#include "settingsnode.h"

namespace flood_processing {
namespace modules {

template<typename T>
class ThresholdMask : public pipeline::Module {
  protected:
    bool less_than;
    bool or_equal;
    bool use_threshold_raster;
    T mask_value;
    T threshold;
    std::string apply_to_name;
    std::string compare_to_name;
    std::string output_name;
    std::string threshold_raster_name;

  public:
    ThresholdMask(const settings::SettingsNode& settings);
    void run(pipeline::Pipeline* p) override;
    inline pipeline::ModuleDescription describe() override {
        if (use_threshold_raster) {
            return pipeline::ModuleDescription{"threshold_mask", {threshold_raster_name, compare_to_name, apply_to_name}, {output_name}};
        } else {
            return pipeline::ModuleDescription{"threshold_mask", {compare_to_name, apply_to_name}, {output_name}};
        }
    }
};

}  // namespace modules
}  // namespace flood_processing

#endif
