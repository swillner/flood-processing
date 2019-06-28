# Copyright (C) 2017 Sven Willner <sven.willner@pik-potsdam.de>
#
# This file is part of flood-processing.
#
# flood-processing is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# flood-processing is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with flood-processing.  If not, see <http://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
include(cmake/helpers.cmake)
include(cmake/libraries.cmake)

set_default_build_type(Release)

project(flood_processing)

file(GLOB MODULE_SOURCES src/modules/*.cpp)
add_executable(flood_processing src/main.cpp ${MODULE_SOURCES})

target_include_directories(flood_processing PRIVATE include lib/lmoments/include lib/cpp-library)
target_compile_options(flood_processing PRIVATE -std=c++11)

set_advanced_cpp_warnings(flood_processing)
set_build_type_specifics(flood_processing)
add_git_version(flood_processing
  FALLBACK_VERSION 1.0.0)

option(PARALLELIZATION "" ON)
if(PARALLELIZATION)
  find_package(OpenMP REQUIRED)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

include_netcdf_cxx4(flood_processing ON v4.3.0)

option(WITH_GDAL "" OFF)
if(WITH_GDAL)
  #include GDAL library
  find_package(GDAL REQUIRED)
  message(STATUS "GDAL include directory: ${GDAL_INCLUDE_DIR}")
  message(STATUS "GDAL library: ${GDAL_LIBRARY}")
  mark_as_advanced(GDAL_CONFIG GDAL_INCLUDE_DIR GDAL_LIBRARY)
  target_link_libraries(flood_processing PRIVATE gdal)
  target_include_directories(flood_processing PRIVATE ${GDAL_INCLUDE_DIR})
  target_compile_definitions(flood_processing PRIVATE FLOOD_PROCESSING_WITH_GDAL)
endif()

include(lib/settingsnode/settingsnode.cmake)
include_settingsnode(flood_processing)
include_yaml_cpp(flood_processing ON "yaml-cpp-0.6.2")

add_cpp_tools(flood_processing)
