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

#ifndef NETCDFGRIDFILE_H
#define NETCDFGRIDFILE_H

#include <assert.h>
#include <ncDim.h>
#include <ncFile.h>
#include <ncGroupAtt.h>
#include <ncType.h>
#include <ncVar.h>
#include <iostream>
#include <netcdf>
#include <string>
#include <vector>

#include "nvector.h"

namespace netCDF {

inline void debugout() { std::cout << std::endl; }
template<typename Arg, typename... Args>
inline void debugout(Arg arg, Args... args) {
    std::cout << arg << " ";
    debugout(args...);
}

template<typename T>
struct NetCDFType {};
template<>
struct NetCDFType<double> {
    static const netCDF::NcType::ncType type = netCDF::NcType::nc_DOUBLE;
};
template<>
struct NetCDFType<float> {
    static const netCDF::NcType::ncType type = netCDF::NcType::nc_FLOAT;
};
template<>
struct NetCDFType<std::int16_t> {
    static const netCDF::NcType::ncType type = netCDF::NcType::nc_USHORT;
};
template<>
struct NetCDFType<const char*> {
    static const netCDF::NcType::ncType type = netCDF::NcType::nc_STRING;
};

template<typename T>
class DimVar : public std::vector<T> {
  protected:
    const std::string name_m;
    std::vector<std::tuple<std::string, NcType, std::size_t, void*>> attributes;

  public:
    using std::vector<T>::size;

    DimVar(NcDim dim, NcVar var) : name_m(dim.getName()), std::vector<T>(dim.getSize()) {
        var.getVar(&(*this)[0]);
        for (const auto& att : var.getAtts()) {
            std::size_t size = att.second.getAttLength();
            void* value = malloc(size);
            att.second.getValues(value);
            attributes.emplace_back(std::make_tuple(att.first, att.second.getType(), size, value));
        }
    }
    DimVar(DimVar&) = delete;
    DimVar(DimVar&&) = default;
    ~DimVar() {
        for (auto& att : attributes) {
            free(std::get<3>(att));
        }
    }
    NcDim write_to(NcFile& file) const {
        netCDF::NcDim result = file.addDim(name_m, size());
        netCDF::NcVar var = file.addVar(name_m, NetCDFType<T>::type, {result});
        var.putVar(&(*this)[0]);
        for (const auto& att : attributes) {
            var.putAtt(std::get<0>(att), std::get<1>(att), std::get<2>(att), std::get<3>(att));
        }
        return result;
    }
    const std::string& name() const { return name_m; }
};

class File : public netCDF::NcFile {
  protected:
    template<typename T, std::size_t c, std::size_t dim, std::size_t... Ns>
    struct get_helper {
        static inline nvector::Vector<T, dim> reserve(const netCDF::NcVar& var_p) { return get_helper<T, c + 1, dim, Ns..., c>::reserve(var_p); }
    };

    template<typename T, std::size_t dim, std::size_t... Ns>
    struct get_helper<T, dim, dim, Ns...> {
        static inline nvector::Vector<T, dim> reserve(const netCDF::NcVar& var_p) { return nvector::Vector<T, dim>(0, var_p.getDim(Ns).getSize()...); }
    };

    template<typename T, std::size_t c, std::size_t dim, std::size_t... Ns>
    struct set_helper {
        static inline void set(netCDF::NcVar& var, const nvector::Vector<T, dim>& data) { set_helper<T, c + 1, dim, Ns..., c>::set(var, data); }
    };

    template<typename T, std::size_t dim, std::size_t... Ns>
    struct set_helper<T, dim, dim, Ns...> {
        static inline void set(netCDF::NcVar& var, const nvector::Vector<T, dim>& data) { var.putVar({(0 * Ns)...}, {data.size(Ns)...}, &data.data()[0]); }
    };

    template<typename T, std::size_t c, std::size_t dim, std::size_t... Ns>
    struct var_set_helper {
        template<typename... Indices>
        static inline void set(netCDF::NcVar& var, const nvector::Vector<T, dim>& data, Indices&&... indices) {
            var_set_helper<T, c + 1, dim, Ns..., c>::set(var, data, std::forward<Indices>(indices)...);
        }
    };

    template<typename T, std::size_t dim, std::size_t... Ns>
    struct var_set_helper<T, dim, dim, Ns...> {
        template<typename... Indices>
        static inline void set(netCDF::NcVar& var, const nvector::Vector<T, dim>& data, Indices&&... indices) {
            var.putVar({std::forward<Indices>(indices)..., (0 * Ns)...}, {(0 * sizeof(Indices) + 1)..., data.size(Ns)...}, &data.data()[0]);
        }
    };

  public:
    File(const std::string& filename, char mode) {
        switch (mode) {
            case 'r':
                open(filename, netCDF::NcFile::read);
                break;
            case 'w':
                open(filename, netCDF::NcFile::replace, netCDF::NcFile::nc4);
                break;
            default:
                throw std::runtime_error("unknown file mode");
        }
    }

    netCDF::NcVar var(const std::string& varname) { return getVar(varname); }

    template<typename T>
    netCDF::NcVar var(const std::string& varname, const std::vector<netCDF::NcDim>& dims) {
        netCDF::NcVar res = addVar(varname, NetCDFType<T>::type, dims);
        res.setCompression(false, true, 7);
        res.setFill<T>(true, std::numeric_limits<T>::quiet_NaN());
        return res;
    }

    template<typename T>
    std::vector<T> get(const std::string& varname) {
        auto var = getVar(varname);
        std::vector<T> result(var.getDim(0).getSize());
        var.getVar(&result[0]);
        return result;
    }

    template<typename T, std::size_t dim>
    nvector::Vector<T, dim> get(const std::string& varname) {
        auto var = getVar(varname);
        auto result = get_helper<T, 0, dim>::reserve(var);
        var.getVar(&result.data()[0]);
        return result;
    }

    template<typename T, std::size_t dim>
    nvector::Vector<T, dim> get(const netCDF::NcVar& var) {
        auto result = get_helper<T, 0, dim>::get(var);
        var.getVar(&result.data()[0]);
        return result;
    }

    template<typename T>
    DimVar<T> dimvar(const netCDF::NcDim& dim) {
        return DimVar<T>(dim, getVar(dim.getName()));
    }

    template<typename T>
    netCDF::NcDim dimvar(const DimVar<T>& dimvar) {
        return dimvar.write_to(*this);
    }

    template<typename T, std::size_t dim>
    nvector::Vector<T, dim> reserve(const std::string& varname) {
        return get_helper<T, 0, dim>::reserve(getVar(varname));
    }

    template<typename T, std::size_t dim>
    nvector::Vector<T, dim> reserve(const netCDF::NcVar& var) {
        return get_helper<T, 0, dim>::reserve(var);
    }

    template<unsigned char dim>
    inline std::size_t size(const netCDF::NcVar& var) const {
        return var.getDim(dim).getSize();
    }

    netCDF::NcDim time(std::size_t steps_count) {  // TODO
        netCDF::NcDim time_dim = addDim("time", steps_count);
        netCDF::NcVar time_var = addVar("time", NetCDFType<double>::type, {time_dim});  // TODO double
        for (std::size_t t = 0; t < steps_count; ++t) {
            time_var.putVar({t}, (double)t);
        }
        return time_dim;
    }

    netCDF::NcDim lat(std::size_t lat_count) {
        netCDF::NcDim lat_dim = addDim("lat", lat_count);
        netCDF::NcVar lat_var = addVar("lat", NetCDFType<double>::type, {lat_dim});
        lat_var.putAtt("standard_name", "latitude");
        lat_var.putAtt("long_name", "latitude");
        lat_var.putAtt("units", "degrees_north");
        lat_var.putAtt("axis", "Y");
        for (unsigned int lat = 0; lat < lat_count; ++lat) {
            lat_var.putVar({lat}, 90.0 - 180.0 * (lat + 0.5) / lat_count);
        }
        return lat_dim;
    }

    netCDF::NcDim lon(std::size_t lon_count) {
        netCDF::NcDim lon_dim = addDim("lon", lon_count);
        netCDF::NcVar lon_var = addVar("lon", NetCDFType<double>::type, {lon_dim});
        lon_var.putAtt("standard_name", "longitude");
        lon_var.putAtt("long_name", "longitude");
        lon_var.putAtt("units", "degrees_east");
        lon_var.putAtt("axis", "X");
        for (unsigned int lon = 0; lon < lon_count; ++lon) {
            lon_var.putVar({lon}, -180.0 + 360.0 * (lon + 0.5) / lon_count);
        }
        return lon_dim;
    }

    template<typename T>
    void set(netCDF::NcVar var, const std::vector<T>& data) {
        var.putVar(&data[0]);
    }

    template<typename T, std::size_t dim>
    void set(netCDF::NcVar var, const nvector::Vector<T, dim>& data) {
        set_helper<T, 0, dim>::set(var, data);
    }

    template<typename T, std::size_t dim, typename... Indices>
    void set(netCDF::NcVar var, const nvector::Vector<T, dim>& data, Indices&&... indices) {
        var_set_helper<T, 0, dim>::set(var, data, std::forward<Indices>(indices)...);
    }

    netCDF::NcDim dim(const netCDF::NcDim& dim) {
        if (dim.isUnlimited()) {
            return addDim(dim.getName());
        } else {
            return addDim(dim.getName(), dim.getSize());
        }
    }

    netCDF::NcVar var(const netCDF::NcVar& var_p) {
        std::vector<netCDF::NcDim> dims;
        std::vector<std::size_t> sizes;
        std::vector<std::size_t> indices;
        std::size_t size = var_p.getType().getSize();
        for (const auto& dim : var_p.getDims()) {
            dims.emplace_back(getDim(dim.getName()));
            sizes.emplace_back(dim.getSize());
            indices.emplace_back(0);
            size *= dim.getSize();
        }
        netCDF::NcVar res = addVar(var_p.getName(), var_p.getType().getTypeClass(), dims);
        res.setCompression(false, true, 7);
        for (const auto& att : var_p.getAtts()) {
            std::size_t size = att.second.getAttLength();
            void* value = malloc(size);
            att.second.getValues(value);
            res.putAtt(att.first, att.second.getType(), size, value);
            free(value);
        }
        void* fill_value = malloc(var_p.getType().getSize());
        void* value = malloc(size);
        bool fill_mode;
        var_p.getFillModeParameters(fill_mode, fill_value);
        if (fill_mode) {
            res.setFill(fill_mode, fill_value);
        }
        var_p.getVar(value);
        res.putVar(indices, sizes, value);
        free(fill_value);
        free(value);
        return res;
    }
};

template<>
inline netCDF::NcVar File::var<const char*>(const std::string& varname, const std::vector<netCDF::NcDim>& dims) {
    netCDF::NcVar res = addVar(varname, NetCDFType<const char*>::type, dims);
    res.setCompression(false, true, 7);
    return res;
}

}  // namespace netCDF

#endif
