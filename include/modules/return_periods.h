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

#ifndef FLOOD_PROCESSING_RETURNPERIODS_H
#define FLOOD_PROCESSING_RETURNPERIODS_H

#include <cstddef>
#include <string>
#include <vector>
#include "nvector.h"
#include "pipeline.h"

namespace settings {
class SettingsNode;
}  // namespace settings

namespace flood_processing {
namespace modules {

template<typename T>
class ReturnPeriods : public pipeline::Module {
  protected:
    enum class Distribution { GEV, GUM };
    Distribution distribution;
    std::vector<std::size_t> from_vec;
    std::vector<std::size_t> to_vec;
    std::size_t length;

    template<typename Distribution>
    inline T anderson_darling_test(const std::vector<T>& data, const Distribution& d);
    template<typename Distribution>
    inline T kolmogorov_smirnov_test(const std::vector<T>& data, const Distribution& d);
    nvector::Vector<T, 3> return_periods(nvector::Vector<T, 3>& history_discharge, nvector::Vector<T, 3>& projection_discharge);

  public:
    explicit ReturnPeriods(const settings::SettingsNode& settings);
    void run(pipeline::Pipeline* p) override;

    inline pipeline::ModuleDescription describe() override {
        return pipeline::ModuleDescription{"return_periods", {"history_discharge", "projection_discharge"}, {"return_periods"}};
    }
};

}  // namespace modules
}  // namespace flood_processing

#endif
