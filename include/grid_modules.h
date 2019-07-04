/*
  Copyright (C) 2017 Sven Willner <sven.willner@gmail.com>

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as published
  by the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef GRID_MODULES_H
#define GRID_MODULES_H

#include <fstream>
#include <stdexcept>
#include <vector>
#include "netcdf/File.h"
#include "nvector.h"
#include "pipeline.h"
#include "settingsnode.h"

namespace pipeline {

inline void check_file_exists(const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile.good()) {
        throw std::runtime_error(filename + " not found");
    }
}

class ArrayReaderModule : public pipeline::Module {
  protected:
    std::string filename;
    std::string varname;
    std::string outputname;
    std::string type_name;

  public:
    explicit ArrayReaderModule(const settings::SettingsNode& settings) {
        filename = settings["filename"].as<std::string>();
        varname = settings["variable"].as<std::string>();
        outputname = settings["output"].as<std::string>();
        type_name = settings["type"].as<std::string>();
    }
    void run(pipeline::Pipeline* p) override {
        check_file_exists(filename);
        netCDF::File file(filename, 'r');
        if (type_name == "string") {
            const auto data = file.get<const char*>(varname);
            auto result = std::make_shared<std::vector<std::string>>(data.size());
            std::transform(std::begin(data), std::end(data), std::begin(*result), [](const char* r) { return r; });
            p->provide<std::vector<std::string>>(outputname, result);
        } else if (type_name == "double") {
            auto result = std::make_shared<std::vector<double>>(file.get<double>(varname));
            p->provide<std::vector<double>>(outputname, result);
        } else {
            throw std::runtime_error("Unknown type '" + type_name + "'");
        }
    }
    inline pipeline::ModuleDescription describe() override { return pipeline::ModuleDescription{"array_reader", {}, {outputname}}; }
};

template<typename T>
class GridReaderModule : public pipeline::Module {
  protected:
    std::string filename;
    std::string varname;
    std::string outputgridname;
    std::string outputtimename;

  public:
    explicit GridReaderModule(const settings::SettingsNode& settings) {
        filename = settings["filename"].as<std::string>();
        varname = settings["variable"].as<std::string>();
        outputgridname = settings["output"]["grid"].as<std::string>("");
        outputtimename = settings["output"]["time"].as<std::string>("");
    }
    void run(pipeline::Pipeline* p) override {
        check_file_exists(filename);
        netCDF::File file(filename, 'r');
        const auto var = file.var(varname);
        if (!outputgridname.empty()) {
            if (check_dimensions(var, {"time", "lat", "lon"}) || check_dimensions(var, {"time", "latitude", "longitude"})) {
                auto grid = std::make_shared<nvector::Vector<T, 3>>(file.get<T, 3>(var));
                p->provide<nvector::Vector<T, 3>>(outputgridname, grid);
            } else if (check_dimensions(var, {"lat", "lon"}) || check_dimensions(var, {"latitude", "longitude"})) {
                auto grid = file.get<T, 2>(var);
                p->provide<nvector::Vector<T, 3>>(outputgridname,
                                                  std::make_shared<nvector::Vector<T, 3>>(grid.data(), 1, grid.template size<0>(), grid.template size<1>()));
            } else {
                throw std::runtime_error(filename + " - " + varname + ": Unexpected dimensions");
            }
        }
        if (!outputtimename.empty()) {
            if (check_dimensions(var, {"time", "lat", "lon"}) || check_dimensions(var, {"time", "latitude", "longitude"})) {
                auto time = std::make_shared<netCDF::DimVar<double>>(file.dimvar<double>(var.getDim(0)));
                p->provide<netCDF::DimVar<double>>(outputtimename, time);
            } else {
                throw std::runtime_error(filename + " - " + varname + ": Unexpected dimensions");
            }
        }
    }
    inline pipeline::ModuleDescription describe() override {
        std::vector<std::string> outputs;
        if (!outputgridname.empty()) {
            outputs.push_back(outputgridname);
        }
        if (!outputtimename.empty()) {
            outputs.push_back(outputtimename);
        }
        return pipeline::ModuleDescription{"grid_reader", {}, outputs};
    }
};

template<typename T>
class TimeSliceModule : public pipeline::Module {
  protected:
    std::string inputname;
    std::string outputname;
    std::size_t i;

  public:
    explicit TimeSliceModule(const settings::SettingsNode& settings) {
        inputname = settings["inputname"].as<std::string>();
        outputname = settings["outputname"].as<std::string>();
        i = settings["index"].as<std::size_t>();
    }
    void run(pipeline::Pipeline* p) override {
        auto input = p->consume<nvector::Vector<T, 3>>(inputname);
        auto grid = input->template split<nvector::Split<false, true, true>>().at(i);
        auto output = std::make_shared<nvector::Vector<T, 3>>(0, 1, grid.template size<0>(), grid.template size<1>());
        std::copy(grid.begin(), grid.end(), output->begin());
        p->provide<nvector::Vector<T, 3>>(outputname, output);
    }
    inline pipeline::ModuleDescription describe() override { return pipeline::ModuleDescription{"time_slice", {inputname}, {outputname}}; }
};

template<typename T>
class GridWriterModule : public pipeline::Module {
  protected:
    std::string filename;
    std::string varname;
    std::string inputgridname;
    std::string inputtimename;

  public:
    explicit GridWriterModule(const settings::SettingsNode& settings) {
        filename = settings["filename"].as<std::string>();
        varname = settings["variable"].as<std::string>();
        inputgridname = settings["input"]["grid"].as<std::string>();
        inputtimename = settings["input"]["time"].as<std::string>("");
    }
    void run(pipeline::Pipeline* p) override {
        auto grid = p->consume<nvector::Vector<T, 3>>(inputgridname);
        netCDF::File file(filename, 'w');
        if (inputtimename.empty()) {
            if (grid->template size<0>() != 1) {
                throw std::runtime_error("Only 2d input grid supported");
            }
            auto var = file.var<T>(varname, {file.lat(grid->template size<1>()), file.lon(grid->template size<2>())});
            file.set<T>(var, *grid);
        } else {
            auto time = p->consume<netCDF::DimVar<double>>(inputtimename);
            if (grid->template size<0>() != time->size()) {
                throw std::runtime_error("Sizes do not match: " + inputgridname + " and " + inputtimename);
            }
            auto var = file.var<T>(varname, {file.dimvar(*time), file.lat(grid->template size<1>()), file.lon(grid->template size<2>())});
            file.set<T>(var, *grid);
        }
    }
    inline pipeline::ModuleDescription describe() override {
        if (inputtimename.empty()) {
            return pipeline::ModuleDescription{"grid_writer", {inputgridname}, {}};
        }
        return pipeline::ModuleDescription{"grid_writer", {inputgridname, inputtimename}, {}};
    }
};

template<typename T>
class RegionRasterGridWriterModule : public GridWriterModule<T> {
  protected:
    using GridWriterModule<T>::filename;
    using GridWriterModule<T>::inputgridname;
    using GridWriterModule<T>::varname;
    std::string regionvarname;
    std::string inputregionsname;

  public:
    explicit RegionRasterGridWriterModule(const settings::SettingsNode& settings) : GridWriterModule<T>(settings) {
        regionvarname = settings["region_variable"].as<std::string>();
        inputregionsname = settings["input"]["regions"].as<std::string>();
    }
    void run(pipeline::Pipeline* p) override {
        auto regions = p->consume<std::vector<std::string>>(inputregionsname);
        auto data = p->consume<nvector::Vector<T, 3>>(inputgridname);
        if (data->template size<0>() != 1) {
            throw std::runtime_error("Only 2d input grid supported");
        }
        netCDF::File file(filename, 'w');
        file.set<T>(file.var<T>(varname, {file.lat(data->template size<1>()), file.lon(data->template size<2>())}), *data);
        std::vector<const char*> regions_char(regions->size());
        std::transform(std::begin(*regions), std::end(*regions), std::begin(regions_char), [](const std::string& r) { return r.c_str(); });
        file.set<const char*>(file.var<const char*>(regionvarname, {file.addDim(regionvarname, regions->size())}), regions_char);
    }
    inline pipeline::ModuleDescription describe() override {
        return pipeline::ModuleDescription{"region_raster_grid_writer", {inputgridname, inputregionsname}, {}};
    }
};

}  // namespace pipeline

#endif
