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

#include <fstream>
#include <iostream>
#include <stdexcept>
#include "grid_modules.h"
#include "modules/downscaling.h"
#include "modules/rasterization.h"
#include "modules/return_level_lookup.h"
#include "modules/return_period_threshold.h"
#include "modules/return_periods.h"
#include "modules/sum_per_region.h"
#include "nvector.h"
#include "pipeline.h"
#include "settingsnode.h"
#include "version.h"

using T = float;

static void run(const settings::SettingsNode& settings) {
    pipeline::Pipeline p;
    for (const auto& m : settings["steps"].as_sequence()) {
        const auto& mod = m["module"].as<settings::hstring>();
        switch (mod) {
            case settings::hstring::hash("array_reader"):
                p.register_module(new pipeline::ArrayReaderModule(m));
                break;
            case settings::hstring::hash("grid_reader2d"):
                p.register_module(new pipeline::GridReader2dModule<T>(m));
                break;
            case settings::hstring::hash("grid_reader3d"):
                p.register_module(new pipeline::GridReader3dModule<T>(m));
                break;
            case settings::hstring::hash("region_raster_grid_writer2d"):
                p.register_module(new pipeline::RegionRasterGridWriter2dModule<T>(m));
                break;
            case settings::hstring::hash("grid3d_slice"):
                p.register_module(new pipeline::Grid3dSliceModule<T>(m));
                break;
            case settings::hstring::hash("grid_writer2d"):
                p.register_module(new pipeline::GridWriter2dModule<T>(m));
                break;
            case settings::hstring::hash("grid_writer3d"):
                p.register_module(new pipeline::GridWriter3dModule<T>(m));
                break;
            case settings::hstring::hash("downscaling"):
                p.register_module(new flood_processing::modules::Downscaling<T>(m));
                break;
            case settings::hstring::hash("return_periods"):
                p.register_module(new flood_processing::modules::ReturnPeriods<T>(m));
                break;
            case settings::hstring::hash("return_level_lookup"):
                p.register_module(new flood_processing::modules::ReturnLevelLookup<T>(m));
                break;
            case settings::hstring::hash("return_period_threshold"):
                p.register_module(new flood_processing::modules::ReturnPeriodThreshold<T>(m));
                break;
            case settings::hstring::hash("rasterization"):
                p.register_module(new flood_processing::modules::Rasterization<T>(m));
                break;
            case settings::hstring::hash("region_index_rasterization"):
                p.register_module(new flood_processing::modules::RegionIndexRasterization<T>(m));
                break;
            case settings::hstring::hash("sum_per_region"):
                p.register_module(new flood_processing::modules::SumPerRegion<T>(m));
                break;
            case settings::hstring::hash("per_region_writer2d"):
                p.register_module(new flood_processing::modules::PerRegionWriter2d<T>(m));
                break;
            default:
                throw std::runtime_error("unknown module '" + std::string(mod) + "'");
        }
    }
    p.run(settings["verbose"].as<bool>(false));
}

static void print_usage(const char* program_name) {
    std::cerr << "Flooding data processing\n"
                 "Version: " FLOOD_PROCESSING_VERSION
                 "\n\n"
                 "Author:  Sven Willner <sven.willner@pik-potsdam.de>\n"
                 "\n"
                 "Usage:   "
              << program_name
              << " (<option> | <settingsfile>)\n"
                 "Options:\n"
                 "  -h, --help     Print this help text\n"
                 "  -v, --version  Print version"
              << std::endl;
}

int main(int argc, char* argv[]) {
#ifndef DEBUG
    try {
#endif
        if (argc != 2) {
            print_usage(argv[0]);
            return 1;
        }
        const std::string arg = argv[1];
        if (arg.length() > 1 && arg[0] == '-') {
            if (arg == "--version" || arg == "-v") {
                std::cout << FLOOD_PROCESSING_VERSION << std::endl;
            } else if (arg == "--help" || arg == "-h") {
                print_usage(argv[0]);
            } else {
                print_usage(argv[0]);
                return 1;
            }
        } else {
            if (arg == "-") {
                std::cin >> std::noskipws;
                run(settings::SettingsNode(std::cin));
            } else {
                std::ifstream settings_file(arg);
                run(settings::SettingsNode(settings_file));
            }
            return 0;
        }
#ifndef DEBUG
    } catch (std::exception& ex) {
        std::cerr << ex.what() << std::endl;
        return 255;
    }
#endif
    return 0;
}
