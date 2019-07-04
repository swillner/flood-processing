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

#ifndef FLOOD_PROCESSING_RETURN_LEVEL_LOOKUP_H
#define FLOOD_PROCESSING_RETURN_LEVEL_LOOKUP_H

#include "pipeline.h"
#include "settingsnode.h"

namespace flood_processing {
namespace modules {

template<typename T>
class ReturnLevelLookup : public pipeline::Module {
  protected:
    bool interpolate;
    std::string return_periods_name;

  public:
    explicit ReturnLevelLookup(const settings::SettingsNode& settings) {
        interpolate = settings["interpolate"].as<bool>();
        return_periods_name = settings["input"]["return_periods"].as<std::string>();
    }
    void run(pipeline::Pipeline* p) override;

    inline pipeline::ModuleDescription describe() override {
        return pipeline::ModuleDescription{"return_level_lookup", {return_periods_name, "return_levels_mapping", "return_periods_mapping"}, {"return_levels"}};
    }
};

}  // namespace modules
}  // namespace flood_processing

#endif
