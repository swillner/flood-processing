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

#ifndef PIPELINE_H
#define PIPELINE_H

#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace pipeline {

class Pipeline;

using DataDescription = std::string;

struct ModuleDescription {
    std::string name;
    std::vector<DataDescription> inputs;
    std::vector<DataDescription> outputs;
};

class Module {
  public:
    virtual void run(Pipeline* p) = 0;
    virtual ModuleDescription describe() = 0;
};

class Pipeline {
  protected:
    class Data {
      protected:
        std::size_t want_count = 0;
        std::shared_ptr<void> data;

      public:
        void want() { ++want_count; }
        bool consume() {
            --want_count;
            return want_count == 0;
        }
        bool is_set() const { return static_cast<bool>(data); }
        template<typename T>
        std::shared_ptr<T> get() {
            return std::static_pointer_cast<T>(data);
        }
        template<typename T>
        void set(std::shared_ptr<T> d) {
            data = d;
        }
    };
    std::unordered_map<DataDescription, Data> warehouse;
    std::vector<std::tuple<std::unique_ptr<Module>, ModuleDescription>> modules;

  public:
    template<typename T>
    std::shared_ptr<T> consume(DataDescription desc) {
        auto data = warehouse.find(desc);
        if (data == warehouse.end()) {
            throw std::runtime_error("input " + desc + " has not been registered");
        }
        if (!data->second.is_set()) {
            throw std::runtime_error(desc + " not set");
        }
        auto result = data->second.get<T>();
        if (data->second.consume()) {
            warehouse.erase(data);
        }
        return result;
    }
    bool is_wanted(DataDescription desc) const { return warehouse.find(desc) != warehouse.end(); }
    template<typename T>
    void provide(DataDescription desc, std::shared_ptr<T> data) {
        auto d = warehouse.find(desc);
        if (d != warehouse.end()) {
            d->second.set<T>(data);
        }
    }
    void register_module(Module* m) {
        ModuleDescription desc = m->describe();
        for (const auto& input : desc.inputs) {
            warehouse[input].want();
        }
        modules.emplace_back(std::make_tuple(std::unique_ptr<Module>(m), desc));
    }
    void run(bool verbose = false) {
        std::size_t i = 0;
        for (const auto& m : modules) {
            if (verbose) {
                std::cout << "[" << (i + 1) << "/" << modules.size() << "] Running module " << std::get<1>(m).name << "..." << std::endl;
            }
            std::get<0>(m)->run(this);
            ++i;
        }
    }
};

}  // namespace pipeline

#endif
