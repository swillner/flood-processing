#  Copyright (C) 2017 Sven Willner <sven.willner@gmail.com>
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as published
#  by the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

set(SETTINGSNODE_MODULES_PATH ${CMAKE_CURRENT_LIST_DIR})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${SETTINGSNODE_MODULES_PATH}/cmake)
include(${SETTINGSNODE_MODULES_PATH}/cmake/libraries.cmake)

function(include_settingsnode TARGET)
  add_library(settingsnode INTERFACE)
  set_property(TARGET settingsnode PROPERTY INTERFACE_COMPILE_DEFINITIONS "SETTINGSNODE_WITH_YAML")
  set_property(TARGET settingsnode PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${SETTINGSNODE_MODULES_PATH}/include)
  target_link_libraries(${TARGET} settingsnode)
  include_yaml_cpp(${TARGET} ON "yaml-cpp-0.6.1")
endfunction()
