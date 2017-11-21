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

#ifndef FLOOD_PROCESSING_RETURN_PERIOD_THRESHOLD_H
#define FLOOD_PROCESSING_RETURN_PERIOD_THRESHOLD_H

#include "pipeline.h"
#include "settingsnode.h"

namespace flood_processing {
namespace modules {

template<typename T>
class ReturnPeriodThreshold : public pipeline::Module {
  protected:
    bool use_raster;
    T threshold;

  public:
    ReturnPeriodThreshold(const settings::SettingsNode& settings) {
        use_raster = settings["threshold"].as<std::string>() == "flopros";
        if (!use_raster) {
            threshold = settings["threshold"].as<T>();
        }
    }
    void run(pipeline::Pipeline* p) override;

    inline pipeline::ModuleDescription describe() override {
        if (use_raster) {
            return pipeline::ModuleDescription{"return_period_threshold", {"raster", "return_periods", "return_levels"}, {"return_levels_thresholded"}};
        } else {
            return pipeline::ModuleDescription{"return_period_threshold", {"return_periods", "return_levels"}, {"return_levels_thresholded"}};
        }
    }
};

}  // namespace modules
}  // namespace flood_processing

#endif
