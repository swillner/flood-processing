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

void check_file_exists(const std::string& filename) {
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
    ArrayReaderModule(const settings::SettingsNode& settings) {
        filename = settings["filename"].as<std::string>();
        varname = settings["varname"].as<std::string>();
        outputname = settings["output"].as<std::string>();
        type_name = settings["type"].as<std::string>();
    }
    void run(pipeline::Pipeline* p) override {
        check_file_exists(filename);
        netCDF::File file(filename, 'r');
        if (type_name == "string") {
            const auto data = file.get<const char*>(varname);
            auto result = std::make_shared<std::vector<std::string>>();
            result->reserve(data.size());
            for (const auto& d : data) {
                result->push_back(d);
            }
            p->provide<std::vector<std::string>>(outputname, result);
        } else if (type_name == "double") {
            const auto data = file.get<double>(varname);
            auto result = std::make_shared<std::vector<double>>();
            result->reserve(data.size());
            for (const auto& d : data) {
                result->push_back(d);
            }
            p->provide<std::vector<double>>(outputname, result);
        }
    }
    inline pipeline::ModuleDescription describe() override { return pipeline::ModuleDescription{"array_reader", {}, {outputname}}; }
};

template<typename T>
class GridReader2dModule : public pipeline::Module {
  protected:
    std::string filename;
    std::string varname;
    std::string outputname;

  public:
    GridReader2dModule(const settings::SettingsNode& settings) {
        filename = settings["filename"].as<std::string>();
        varname = settings["varname"].as<std::string>();
        outputname = settings["output"].as<std::string>();
    }
    void run(pipeline::Pipeline* p) override {
        check_file_exists(filename);
        netCDF::File file(filename, 'r');
        p->provide<nvector::Vector<T, 2>>(outputname, std::make_shared<nvector::Vector<T, 2>>(file.get<T, 2>(varname)));
    }
    inline pipeline::ModuleDescription describe() override { return pipeline::ModuleDescription{"grid_reader2d", {}, {outputname}}; }
};

template<typename T>
class GridReader3dModule : public pipeline::Module {
  protected:
    std::string filename;
    std::string varname;
    std::string outputgridsname;
    std::string outputdim1name;

  public:
    GridReader3dModule(const settings::SettingsNode& settings) {
        filename = settings["filename"].as<std::string>();
        varname = settings["varname"].as<std::string>();
        outputgridsname = settings["output"]["grids"].as<std::string>("");
        outputdim1name = settings["output"]["dim1"].as<std::string>("");
    }
    void run(pipeline::Pipeline* p) override {
        check_file_exists(filename);
        netCDF::File file(filename, 'r');
        if (!outputgridsname.empty()) {
            auto grids = std::make_shared<nvector::Vector<T, 3>>(file.get<T, 3>(varname));
            p->provide<nvector::Vector<T, 3>>(outputgridsname, grids);
        }
        if (!outputdim1name.empty()) {
            auto dim1 = std::make_shared<netCDF::DimVar<double>>(file.dimvar<double>(file.var(varname).getDim(0)));
            p->provide<netCDF::DimVar<double>>(outputdim1name, dim1);
        }
    }
    inline pipeline::ModuleDescription describe() override {
        std::vector<std::string> outputs;
        if (!outputgridsname.empty()) {
            outputs.push_back(outputgridsname);
        }
        if (!outputdim1name.empty()) {
            outputs.push_back(outputdim1name);
        }
        return pipeline::ModuleDescription{"grid_reader3d", {}, outputs};
    }
};

template<typename T>
class Grid3dSliceModule : public pipeline::Module {
  protected:
    std::string inputname;
    std::string outputname;
    std::size_t i;

  public:
    Grid3dSliceModule(const settings::SettingsNode& settings) {
        inputname = settings["inputname"].as<std::string>();
        outputname = settings["outputname"].as<std::string>();
        i = settings["index"].as<std::size_t>();
    }
    void run(pipeline::Pipeline* p) override {
        auto input = p->consume<nvector::Vector<T, 3>>(inputname);
        auto grid = input->template split<nvector::Split<false, true, true>>().at(i);
        auto output = std::make_shared<nvector::Vector<T, 2>>(0, grid.template size<0>(), grid.template size<1>());
        std::copy(grid.begin(), grid.end(), output->begin());
        p->provide<nvector::Vector<T, 2>>(outputname, output);
    }
    inline pipeline::ModuleDescription describe() override { return pipeline::ModuleDescription{"grid3d_slice", {inputname}, {outputname}}; }
};

template<typename T>
class GridWriter2dModule : public pipeline::Module {
  protected:
    std::string filename;
    std::string varname;
    std::string inputname;

  public:
    GridWriter2dModule(const settings::SettingsNode& settings) {
        filename = settings["filename"].as<std::string>();
        varname = settings["varname"].as<std::string>();
        inputname = settings["input"].as<std::string>();
    }
    void run(pipeline::Pipeline* p) override {
        auto data = p->consume<nvector::Vector<T, 2>>(inputname);
        netCDF::File file(filename, 'w');
        file.set<T>(file.var<T>(varname, {file.lat(data->template size<0>()), file.lon(data->template size<1>())}), *data);
    }
    inline pipeline::ModuleDescription describe() override { return pipeline::ModuleDescription{"grid_writer2d", {inputname}, {}}; }
};

template<typename T>
class RegionRasterGridWriter2dModule : public GridWriter2dModule<T> {
  protected:
    using GridWriter2dModule<T>::filename;
    using GridWriter2dModule<T>::inputname;
    using GridWriter2dModule<T>::varname;
    std::string regionvarname;

  public:
    RegionRasterGridWriter2dModule(const settings::SettingsNode& settings) : GridWriter2dModule<T>(settings) {
        regionvarname = settings["regionvarname"].as<std::string>();
    }
    void run(pipeline::Pipeline* p) override {
        auto regions = p->consume<std::vector<std::string>>("regions");
        auto data = p->consume<nvector::Vector<T, 2>>(inputname);
        netCDF::File file(filename, 'w');
        file.set<T>(file.var<T>(varname, {file.lat(data->template size<0>()), file.lon(data->template size<1>())}), *data);
        std::vector<const char*> regions_char;
        regions_char.reserve(regions->size());
        for (const auto& r : *regions) {
            regions_char.push_back(r.c_str());
        }
        file.set<const char*>(file.var<const char*>(regionvarname, {file.addDim(regionvarname, regions->size())}), regions_char);
    }
    inline pipeline::ModuleDescription describe() override { return pipeline::ModuleDescription{"region_raster_grid_writer2d", {inputname, "regions"}, {}}; }
};

template<typename T>
class GridWriter3dModule : public pipeline::Module {
  protected:
    std::string filename;
    std::string varname;
    std::string inputgridsname;
    std::string inputdim1name;

  public:
    GridWriter3dModule(const settings::SettingsNode& settings) {
        filename = settings["filename"].as<std::string>();
        varname = settings["varname"].as<std::string>();
        inputgridsname = settings["input"]["grids"].as<std::string>();
        inputdim1name = settings["input"]["dim1"].as<std::string>();
    }
    void run(pipeline::Pipeline* p) override {
        auto grids = p->consume<nvector::Vector<T, 3>>(inputgridsname);
        auto dim1 = p->consume<netCDF::DimVar<double>>(inputdim1name);
        netCDF::File file(filename, 'w');
        if (grids->template size<0>() != dim1->size()) {
            throw std::runtime_error("Sizes do not match: " + inputgridsname + " and " + inputdim1name);
        }
        netCDF::NcVar var = file.var<T>(varname, {file.dimvar(*dim1), file.lat(grids->template size<1>()), file.lon(grids->template size<2>())});
        file.set<T>(var, *grids);
    }
    inline pipeline::ModuleDescription describe() override { return pipeline::ModuleDescription{"grid_writer3d", {inputgridsname, inputdim1name}, {}}; }
};

}  // namespace pipeline

#endif
