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

#ifndef FLOOD_PROCESSING_SUM_PER_REGION_H
#define FLOOD_PROCESSING_SUM_PER_REGION_H

#include <algorithm>
#include "netcdf/File.h"
#include "nvector.h"
#include "pipeline.h"
#include "settingsnode.h"

namespace flood_processing {
namespace modules {

// TODO: make input names controllable in settings

template<typename T>
class SumPerRegion : public pipeline::Module {
  protected:
    std::string inputname;
    std::string outputname;
    std::string filename;
    std::string varname;

  public:
    explicit SumPerRegion(const settings::SettingsNode& settings) {
        outputname = settings["outputname"].as<std::string>();
        filename = settings["filename"].as<std::string>("");
        if (!filename.empty()) {
            varname = settings["variable"].as<std::string>();
        } else {
            inputname = settings["inputname"].as<std::string>();
        }
    }
    void run(pipeline::Pipeline* p) override;

    inline pipeline::ModuleDescription describe() override {
        if (filename.empty()) {
            return pipeline::ModuleDescription{"sum_per_region", {inputname, "region_index_raster", "regions"}, {outputname}};
        }
        return pipeline::ModuleDescription{"sum_per_region", {"region_index_raster", "regions"}, {outputname}};
    }
};

template<typename T>
class PerRegionWriter2d : public pipeline::Module {
  protected:
    std::string filename;
    std::string inputarrayname;
    std::string inputdim1name;
    std::string regionvarname;
    std::string varname;

  public:
    explicit PerRegionWriter2d(const settings::SettingsNode& settings) {
        filename = settings["filename"].as<std::string>();
        inputarrayname = settings["input"]["array"].as<std::string>();
        inputdim1name = settings["input"]["time"].as<std::string>();
        regionvarname = settings["regionvarname"].as<std::string>();
        varname = settings["variable"].as<std::string>();
    }
    void run(pipeline::Pipeline* p) override {
        auto regions = p->consume<std::vector<std::string>>("regions");
        auto data = p->consume<nvector::Vector<T, 2>>(inputarrayname);
        auto dim1 = p->consume<netCDF::DimVar<double>>(inputdim1name);
        netCDF::File file(filename, 'w');

        std::vector<const char*> regions_char(regions->size());
        std::transform(std::begin(*regions), std::end(*regions), std::begin(regions_char), [](const std::string& r) { return r.c_str(); });
        auto regiondim = file.addDim(regionvarname, regions->size());
        file.set<const char*>(file.var<const char*>(regionvarname, {regiondim}), regions_char);

        if (data->template size<0>() != dim1->size()) {
            throw std::runtime_error("Sizes do not match: " + inputarrayname + " and " + inputdim1name);
        }
        if (data->template size<1>() != regions->size()) {
            throw std::runtime_error("Sizes do not match: " + inputarrayname + " and regions");
        }
        netCDF::NcVar var = file.var<T>(varname, {file.dimvar(*dim1), regiondim});
        file.set<T>(var, *data);
    }
    inline pipeline::ModuleDescription describe() override { return pipeline::ModuleDescription{"per_region_writer", {inputarrayname, "regions", "time"}, {}}; }
};

}  // namespace modules
}  // namespace flood_processing

#endif
